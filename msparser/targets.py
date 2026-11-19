"""Target ion generation from user-entered components.

Rules:
- Formulas are neutral. The user lists the allowed charge states (+1, +2, -1, -2).
- mode "proton": charge by gain/loss of H+.
- mode "intrinsic": charge by electron removal/addition only; the component
  also acts as a charge carrier for adducts with proton-mode molecules.
- Adducts use entered species only: a molecule with entered intrinsic ions, and
  pairs of two different components (two molecules sharing H+, or a cation
  with an anion). No default adducts, neutral losses, radical ions or
  combinations of a component with itself.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from molmass import ELECTRON, ELEMENTS, Formula

ELECTRON_MASS = ELECTRON.mass


@dataclass(frozen=True)
class Component:
    name: str
    formula: str
    charges: tuple[int, ...]
    mode: str = "proton"

    def composition(self) -> Counter:
        f = Formula(self.formula)
        if f.charge != 0:
            raise ValueError(f"{self.name}: enter a neutral formula, charge is set separately")
        comp = Counter()
        for sym, item in f.composition().asdict().items():
            if sym not in ELEMENTS:
                raise ValueError(f"{self.name}: unsupported symbol '{sym}' in formula")
            comp[sym] = int(item[0])
        return comp


def normalize_mode(mode) -> str:
    """'molecule' and 'ion' are accepted as synonyms of 'proton' and 'intrinsic'."""
    m = str(mode).strip().lower()
    return {"molecule": "proton", "ion": "intrinsic"}.get(m, m)


def load_components(path: str | Path) -> list[Component]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    comps, seen = [], set()
    for d in data["components"]:
        charges = tuple(int(c) for c in d.get("charges", []))
        if any(c not in (-2, -1, 1, 2) for c in charges):
            raise ValueError(f"{d['name']}: charges must be among -2, -1, +1, +2")
        mode = normalize_mode(d.get("mode", "proton"))
        if mode not in ("proton", "intrinsic"):
            raise ValueError(f"{d['name']}: mode must be 'molecule' (or 'proton') or 'ion' (or 'intrinsic')")
        if d["name"] in seen:
            raise ValueError(f"duplicate component name {d['name']}")
        seen.add(d["name"])
        comps.append(Component(str(d["name"]), str(d["formula"]), charges, mode))
    return comps


def _light_isotope(sym: str):
    return max(ELEMENTS[sym].isotopes.values(), key=lambda i: i.abundance)


def mono_mass(comp: Counter) -> float:
    return sum(n * _light_isotope(s).mass for s, n in comp.items())


def ion_mz(comp: Counter, z: int) -> float:
    return (mono_mass(comp) - z * ELECTRON_MASS) / abs(z)


def best_isotopologue(comp: Counter, z: int) -> tuple[float, float]:
    """Most abundant single-substitution isotopologue: (m/z shift, ratio to monoisotopic)."""
    best = (np.nan, np.nan)
    for sym, n in comp.items():
        if n <= 0:
            continue
        light = _light_isotope(sym)
        for iso in ELEMENTS[sym].isotopes.values():
            if iso.massnumber == light.massnumber or iso.abundance <= 0:
                continue
            ratio = n * iso.abundance / light.abundance
            if not best[1] >= ratio:
                best = ((iso.mass - light.mass) / abs(z), ratio)
    return best


def _charge_suffix(z: int) -> str:
    return ("" if abs(z) == 1 else str(abs(z))) + ("+" if z > 0 else "-")


def _label(parts: Counter, z: int) -> str:
    s = "[M"
    for (sign, name), n in parts.items():
        s += sign + (str(n) if n > 1 else "") + name
    return s + "]" + _charge_suffix(z)


def _add_pairs(components: list[Component], add) -> None:
    """Ions made of two different components."""
    for i, a in enumerate(components):
        for b in components[i + 1:]:
            name, req = f"{a.name} + {b.name}", {a.name, b.name}
            both = a.composition() + b.composition()
            if a.mode == "proton" and b.mode == "proton":
                # two molecules sharing H+: [M1+M2+H]+, [M1+M2-H]-, and the doubly charged forms
                for z in (1, 2, -1, -2):
                    s = 1 if z > 0 else -1
                    ok = z in a.charges or z in b.charges or (abs(z) == 2 and s in a.charges and s in b.charges)
                    comp = Counter(both)
                    comp["H"] += z
                    if not ok or comp["H"] < 0:
                        continue
                    h = ("+" if z > 0 else "-") + ("2" if abs(z) == 2 else "") + "H"
                    add(name, f"[M1+M2{h}]" + _charge_suffix(z), z, comp, req, "pair")
            elif a.mode == "intrinsic" and b.mode == "intrinsic":
                # cation with anion, when the charges do not cancel: e.g. [Ca+Cl]+
                for qa in a.charges:
                    for qb in b.charges:
                        z = qa + qb
                        if qa * qb < 0 and z != 0 and abs(z) <= 2:
                            add(name, "[M1+M2]" + _charge_suffix(z), z, both, req, "pair")


def generate_targets(components: list[Component], pairs: bool = True) -> pd.DataFrame:
    intrinsic = [c for c in components if c.mode == "intrinsic"]
    rows = []

    def add(comp_name, ion, z, comp, requires, kind="ion"):
        shift, ratio = best_isotopologue(comp, z)
        mz = ion_mz(comp, z)
        rows.append(
            dict(
                target_id=f"{comp_name} {ion}", component=comp_name, ion=ion, charge=z, type=kind,
                mz=mz, iso_mz=mz + shift, iso_ratio=ratio,
                requires=";".join(sorted(requires)), use=1,
            )
        )

    for c in components:
        base = c.composition()
        if c.mode == "intrinsic":
            for z in c.charges:
                add(c.name, "[M]" + _charge_suffix(z), z, base, {c.name})
            continue
        for z in c.charges:
            s = 1 if z > 0 else -1
            # carriers: (label sign, label name, composition delta, charge, required component)
            carriers = [("+" if s > 0 else "-", "H", Counter({"H": s}), s, None)]
            for x in intrinsic:
                for q in x.charges:
                    if q * s > 0 and abs(q) <= abs(z):
                        carriers.append(("+", x.formula, x.composition(), q, x.name))
            for k in range(1, abs(z) + 1):
                for combo in combinations_with_replacement(range(len(carriers)), k):
                    if sum(carriers[i][3] for i in combo) != z:
                        continue
                    comp, parts, req = Counter(base), Counter(), {c.name}
                    for i in combo:
                        sign, name, delta, _, need = carriers[i]
                        for sym, n in delta.items():
                            comp[sym] += n
                        parts[(sign, name)] += 1
                        if need:
                            req.add(need)
                    if min(comp.values()) < 0:
                        continue  # not enough H to remove
                    add(c.name, _label(parts, z), z, comp, req, "adduct" if len(req) > 1 else "ion")
    if pairs:
        _add_pairs(components, add)
    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError("no targets generated: check charges of the components")
    df = df.drop_duplicates("target_id")
    order = df["type"].map({"ion": 0, "adduct": 1, "pair": 2})
    return df.assign(_o=order).sort_values(["_o", "component", "charge", "mz"]).drop(columns="_o").reset_index(drop=True)


def targets_for_file(targets: pd.DataFrame, injected: set[str] | None, polarity: int) -> pd.DataFrame:
    """Targets applicable to one acquisition: confirmed, matching polarity, all
    required components injected."""
    t = targets[(targets["use"] == 1) & (np.sign(targets["charge"]) == polarity)]
    if injected is not None:
        ok = t["requires"].map(lambda r: set(r.split(";")) <= injected)
        t = t[ok]
    return t.reset_index(drop=True)

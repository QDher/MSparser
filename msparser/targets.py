"""Target ion generation from user-entered components.

Rules:
- Formulas are neutral. The user lists the allowed charge states (-3 to +3).
- mode "proton": charge by gain/loss of H+.
- mode "intrinsic": charge by electron removal/addition only; the component
  also acts as a charge carrier for adducts with proton-mode molecules.
- Adducts use entered species only. One molecule, or two different molecules,
  gain or lose up to 3 protons and attach up to 3 entered intrinsic ions
  (including exchanges such as [M-2H+Na]-). Protons and attached ions do not
  count as molecules; an ion never holds more than two molecules, nor the
  same molecule twice. A cation with an anion is also proposed. No default
  adducts, neutral losses or radical ions.
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
        if any(c not in (-3, -2, -1, 1, 2, 3) for c in charges):
            raise ValueError(f"{d['name']}: charges must be between -3 and +3 (not 0)")
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


def _h_part(h: int) -> str:
    return "" if h == 0 else ("+" if h > 0 else "-") + (str(abs(h)) if abs(h) > 1 else "") + "H"


def _label(core: str, h: int, ions: Counter, z: int) -> str:
    """[M+H]+, [M-2H+Na]-, [M+H+Na]2+, [M1+M2-H]-, [M1+M2+Na]+ ..."""
    s = "[" + core + _h_part(h)
    for name, n in ions.items():
        s += "+" + (str(n) if n > 1 else "") + name
    return s + "]" + _charge_suffix(z)


# Limits of one ion. Molecules are the components that gain or lose H+; the
# protons themselves and the entered ions do not count as molecules.
MAX_MOLECULES = 2     # molecules in one ion (never the same one twice)
MAX_H = 3             # protons gained or lost
MAX_ADDUCT_IONS = 3   # entered ions attached


def generate_targets(components: list[Component], pairs: bool = True) -> pd.DataFrame:
    molecules = [c for c in components if c.mode == "proton"]
    intrinsic = [c for c in components if c.mode == "intrinsic"]
    carriers = [(x, q) for x in intrinsic for q in x.charges]  # entered ions, with each of their charges
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

    def with_protons_and_ions(name, core, base, charges, req0, alone):
        """One or two molecules reach each allowed charge z by gaining or losing
        h protons and attaching entered ions, cations or anions (not both at once):
        [M+H]+, [M+Na]+, [M+H+Na]2+, [M-2H+Na]-, [M-H+2Na]+, [M+Cl]-, [M1+M2+Na]+ ..."""
        for z in charges:
            for k in range(MAX_ADDUCT_IONS + 1):
                for combo in combinations_with_replacement(range(len(carriers)), k):
                    qs = [carriers[i][1] for i in combo]
                    if any(q > 0 for q in qs) and any(q < 0 for q in qs):
                        continue
                    h = z - sum(qs)
                    if abs(h) > MAX_H or (k == 0 and h == 0):
                        continue
                    comp, ions, req = Counter(base), Counter(), set(req0)
                    comp["H"] += h
                    for i in combo:
                        x = carriers[i][0]
                        comp.update(x.composition())
                        ions[x.formula] += 1
                        req.add(x.name)
                    if comp["H"] < 0:
                        continue  # not enough H to remove
                    kind = ("adduct" if k else "ion") if alone else "pair"
                    add(name, _label(core, h, ions, z), z, comp, req, kind)

    for c in intrinsic:
        for z in c.charges:
            add(c.name, "[M]" + _charge_suffix(z), z, c.composition(), {c.name})
    for c in molecules:
        with_protons_and_ions(c.name, "M", c.composition(), c.charges, {c.name}, alone=True)
    if pairs and MAX_MOLECULES >= 2:
        for i, a in enumerate(molecules):
            for b in molecules[i + 1:]:
                # charges of the pair: those of either molecule, or the sum of one charge of each
                charges = set(a.charges) | set(b.charges) | {
                    za + zb for za in a.charges for zb in b.charges if za * zb > 0 and abs(za + zb) <= 3}
                with_protons_and_ions(f"{a.name} + {b.name}", "M1+M2", a.composition() + b.composition(),
                                      sorted(charges), {a.name, b.name}, alone=False)
        for i, a in enumerate(intrinsic):
            for b in intrinsic[i + 1:]:
                # cation with anion, when the charges do not cancel: e.g. [Ca+Cl]+
                for qa in a.charges:
                    for qb in b.charges:
                        z = qa + qb
                        if qa * qb < 0 and z != 0 and abs(z) <= 3:
                            add(f"{a.name} + {b.name}", "[M1+M2]" + _charge_suffix(z), z,
                                a.composition() + b.composition(), {a.name, b.name}, "pair")
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

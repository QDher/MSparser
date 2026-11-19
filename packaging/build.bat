@echo off
rem Build the Windows application. Run from the repository root:  packaging\build.bat
rem Requires Python 3.10+ on the PATH. Result: dist\MSparser\MSparser.exe
python -m venv .build-venv || exit /b 1
call .build-venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install ".[ui]" pyinstaller || exit /b 1
pyinstaller --noconfirm packaging\msparser.spec || exit /b 1
echo.
echo Built: dist\MSparser\MSparser.exe
echo Distribute the whole dist\MSparser folder (zip it).

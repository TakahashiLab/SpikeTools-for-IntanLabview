@echo off
setlocal
if defined PREPHY_PYTHON goto run
if exist "%~dp0.venv\Scripts\python.exe" (
  set "PREPHY_PYTHON=%~dp0.venv\Scripts\python.exe"
  goto run
)
if exist "%USERPROFILE%\.conda\envs\phy2\python.exe" (
  set "PREPHY_PYTHON=%USERPROFILE%\.conda\envs\phy2\python.exe"
  goto run
)
echo Set PREPHY_PYTHON to your Python executable, or create .venv in this folder.
echo See PREPHY_MIGRATION.md for setup instructions.
pause
exit /b 1
:run
if not exist "%PREPHY_PYTHON%" (
  echo Python was not found: "%PREPHY_PYTHON%"
  pause
  exit /b 1
)
"%PREPHY_PYTHON%" "%~dp0prephy_gui.py" %*
if errorlevel 1 (
  pause
  exit /b 1
)

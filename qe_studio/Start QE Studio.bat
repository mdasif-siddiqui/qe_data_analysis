@echo off
rem Double-click to start QE Studio (runs in WSL Ubuntu) and open it in the browser.
cd /d "%~dp0"

rem Check that the default WSL system can run bash (fails if WSL/Ubuntu is missing or another distro is default)
wsl -e bash -c "exit 0" >nul 2>&1
if errorlevel 1 (
  echo ERROR: could not start Ubuntu in WSL.
  echo Check in PowerShell:   wsl -l -v
  echo Ubuntu must be listed and marked with * as the default. If not:   wsl --set-default Ubuntu
  echo If WSL is not installed:   wsl --install -d Ubuntu
  pause
  exit /b 1
)

start "" http://localhost:8050
rem Remove Windows line endings from the scripts (they break bash), then start
wsl -e bash -c "d=\"$(wslpath '%~dp0')\"; sed -i 's/\r$//' \"$d/start.sh\" \"$d/setup_wsl.sh\" \"$d/build_gipaw.sh\"; bash \"$d/start.sh\""
pause

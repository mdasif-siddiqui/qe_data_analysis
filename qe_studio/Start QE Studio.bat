@echo off
rem Double-click to start QE Studio (runs in WSL Ubuntu) and open it in the browser.
cd /d "%~dp0"
start "" http://localhost:8050
wsl -e bash -c "bash \"$(wslpath '%~dp0')start.sh\""
pause

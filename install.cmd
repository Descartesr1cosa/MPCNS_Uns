@echo off
rem Windows (Developer Command Prompt / MSYS2 + MS-MPI): build and create the run folders of ONE model.
rem   usage: install.cmd C:\path\to\work [Mercury|Mars]
if "%~1"=="" ( echo usage: %0 ^<work-dir^> [Mercury^|Mars] & exit /b 2 )
set MODEL=%~2
if "%MODEL%"=="" set MODEL=Mercury
cmake -S "%~dp0." -B "%~dp0build" -DCMAKE_BUILD_TYPE=Release -DUNS_BUILD_TESTS=OFF -DUNS_CASE_MODEL=%MODEL% || exit /b 1
cmake --build "%~dp0build" --config Release -j || exit /b 1
cmake --install "%~dp0build" --config Release --prefix "%~1"

@echo off
setlocal

if exist "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars64.bat" (
	call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars64.bat" > nul
	goto run
)

if exist "C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\VC\Auxiliary\Build\vcvars64.bat" (
	call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\VC\Auxiliary\Build\vcvars64.bat" > nul
	goto run
)

echo Cannot setup environment
exit /b 1

:run
if not exist ..\bin (
  mkdir ..\bin
	copy plot.* ..\bin
)
set CommonCompilerFlags=-nologo -W4 -Od -Zi -wd4996
cl %CommonCompilerFlags% -Fe..\bin\main.exe main.cpp

@echo off
rem throw the first parameter away
set filedir=%~dp0
set params=%1
:loop
shift
if [%1]==[] goto afterloop
set params=%params% %1
goto loop
:afterloop

Rscript.exe "%filedir%\landscapes_and_ML.R" %params%
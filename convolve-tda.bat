@echo off
rem throw the first parameter away
set filedir=%~dp0
shift
set params=%0
:loop
shift
if [%1]==[] goto afterloop
set params=%params% %1
goto loop
:afterloop

Rscript "%filedir%\convolve_images.R" %params%
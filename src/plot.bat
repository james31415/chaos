@echo off
main > out.csv
if errorlevel 1 (
  echo Plotting failed 
	goto end
)
gnuplot -c plot.gp
:end

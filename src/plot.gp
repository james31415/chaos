set datafile separator ","
splot "out.csv" using 2:3:4 with lines lc rgb "#ff0000"
pause -1

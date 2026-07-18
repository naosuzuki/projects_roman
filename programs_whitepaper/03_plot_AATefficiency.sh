cat -> gnubody << eof
set term png
set output 'AATefficiency.png'
set xr[20.0:26.0]
set yr[0:1.2]
set xlabel 'Mag'
set ylabel 'AAT Efficiency'
plot '../data/AATmaglimit_v2.dat' pt 4,1.0+0.871884*x-0.140627*x**2+0.00751911*x**3-0.000133365*x**4 lt 7 lw 4
eof

gnuplot gnubody
rm -f gnubody


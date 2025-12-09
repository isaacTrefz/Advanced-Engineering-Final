set term gif animate delay 10 size 1200,400
set output "fluid-field.gif"
set xlabel "x (m)"
set ylabel "y (m)"
set xrange [0:4]
set yrange [0:1]
#set size ratio 0.4
set size {1.0, 0.25}
set key off
#set object 1 rect from 0.3,0.3 to 0.4,0.7 fc rgb "blue" fillstyle solid front

set object 1 polygon from 0.625,0.375 to 1,0.75 to 1.125,0.625 to 0.75,0.25 to 0.625,0.375
set object 1 fillstyle solid 1.0 fillcolor "purple"


Time = 0.1
SKIP = 1
NT = 100
dt = Time/NT

do for [k=0:NT/SKIP] {
	t = SKIP*(k+1)*dt
	plot "Final-fluid-field.dat" index k u 1:2:3:4 every 2:2 w vectors head filled
  set title sprintf("2D Navier-Stokes, nu = 0.3, t = %.3f s",t)
}

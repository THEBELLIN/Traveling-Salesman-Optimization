set terminal wxt size 1920, 1080
set term wxt font ",1"
set style line 1 \
linecolor rgb '#FFaa00' \
linetype 1 linewidth 2 \
pointtype 7 pointsize 1
set title 'plot points'
plot "../data/data_points.dat" \
using 1:2:(sprintf("%s",strcol(3))) \
with labels offset 0.5,0.5 \
point pt 7 ps 1 lc rgb '#FFaa00'


set style line 1 \
linecolor rgb '#FFaa00' \
linetype 1 linewidth 2 \
pointtype 7 pointsize 1 
replot "../data/data_lines.dat" with linespoints linestyle 1
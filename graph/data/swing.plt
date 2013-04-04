#set style line 1 lc rgb '#0060ad' lt 1 lw 1 pt 1 ps 1   # --- blue
plot 'swing.dat' using 1:2 with lines lt 1 title 'P1', \
'swing.dat' using 1:3 with lines lt 2 title 'P2', \
'swing.dat' using 1:4 with lines lt 3 title 'a1', \
'swing.dat' using 1:5 with lines lt 4 title 'a2'\


set key left box
set samples 150
plot [0:2*3.1415926] cos(x), 1.05714 -1.09135*x + 0.173693*x*x, 1+0.424413*x-1.14831*x*x+0.344016*x*x*x-0.027376*x*x*x*x
#plot [0:2*3.1415926] cos(x), 1.00001-0.036657*x-0.390086*x*x-0.122929*x*x*x+0.108838*x*x*x*x-0.0189183*x*x*x*x*x+0.00100365*x*x*x*x*x*x
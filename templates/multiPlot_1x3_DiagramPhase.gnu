#!/usr/bin/gnuplot
#
# Demonstration of a common use scenario of the multiplot environment.
#
# AUTHOR: Hagen Wierstorf
#

reset

# wxt terminal
#set terminal wxt size 900,600 enhanced font 'Verdana,10' persist
# png
set terminal pngcairo size 900,300 enhanced font 'Verdana,10'
set output 'layout.6.png'
# svg
#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif' \
#fsize '10'
#set output 'multiplot3.svg'

# color definitions
set style line 1 lc rgb 'black' pt 7   # circle
set style line 2 lc rgb '#E8991E'
set style line 3 lc rgb '#070BE8'
set style line 4 lc rgb '#8313E8'
unset key

# Enable the use of macros
set macros

set tics scale 1.0
set ytics 1
set xrange [0:1]
set yrange [0:4.5]

# MACROS
# x- and ytics for each row resp. column
NOXTICS = "set xtics ('' 0, '' 0.2,'' 0.4,'' 0.6, '' 0.8, '' 1); \
          unset xlabel"
XTICS = "set xtics ('0' 0, '0.2' 0.2, '0.4' 0.4,'0.6' 0.6, '0.8' 0.8, '1' 1);\
          set xlabel 'J_2 / J_1'"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%.0f'; set ylabel 'h / J_1'"
# Placement of the a,b,c,d labels in the graphs
POS = "at graph 0.03,0.95 font ',10'"
POS_J3 = "at graph 0.75,0.94 font ',10'"
POS_PM = "at graph 0.49,0.7 font ',12'"
POS_AF = "at graph 0.175,0.3 font ',12'"
POS_SAF = "at graph 0.76,0.3 font ',12'"
#POS_SD = "at graph 0.48,0.13 font ',12'"

# Margins for each row resp. column
TMARGIN = "set tmargin at screen .96; set bmargin at screen 0.20"
BMARGIN = "set tmargin at screen 0.54; set bmargin at screen 0.14"
LMARGIN = "set lmargin at screen 0.05; set rmargin at screen 0.35"
CMARGIN = "set lmargin at screen 0.37; set rmargin at screen 0.67"
RMARGIN = "set lmargin at screen 0.69; set rmargin at screen 0.99"


### Start multiplot (1x2 layout)
set multiplot layout 1,3 rowsfirst

# --- GRAPH a
AF = "./DiagramasEstadoFundamental/J3(0.00)/AF_mZero-Gamma-J2.dat"
SAF = "./DiagramasEstadoFundamental/J3(0.00)/SAF_mZero-Gamma-J2.dat"
AF_SAF = "./DiagramasEstadoFundamental/J3(0.00)/Transition--AF-SAF.dat"

@XTICS; @YTICS; @TMARGIN; ; @LMARGIN
set label 1 'a)' @POS
set label 2 'J_{3} = 0.0' @POS_J3
set label 3 "PM" @POS_PM 
set label 4 "AF" @POS_AF 
set label 5 "SAF" @POS_SAF 
set label 6 at 0.475, 1.95
set label 6 "0.56" tc rgb '#8313E8'

plot AF using ($3<=0.51?$3:1/0):2 w l ls 3, \
    SAF using ($3>0.56?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.56 && $3>=0.505?$3:1/0):2 w l ls 3 dt 2, \
    "<echo '.56 1.7'" with points ls 1, \
    AF_SAF using 1:($2<=1.5?$2:1/0) w l ls 3 dt 2



# --- GRAPH b
AF = "./DiagramasEstadoFundamental/J3(-0.05)/AF_mZero-Gamma-J2.dat"
SAF = "./DiagramasEstadoFundamental/J3(-0.05)/SAF_mZero-Gamma-J2.dat"
AF_SAF = "./DiagramasEstadoFundamental/J3(-0.05)/Transition--AF-SAF.dat"

@XTICS; @NOYTICS; @TMARGIN; @CMARGIN 
set label 1 'b)' @POS
set label 2 'J_{3} = -0.05' @POS_J3
set label 6 at 0.46, 2.1
set label 6 "0.52" tc rgb '#8313E8'

plot AF using ($3<=0.51?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.52?$3:1/0):2 w l ls 3, \
    AF_SAF using 1:($2<=1.7?$2:1/0) w l ls 3 dt 2, \
    "<echo '.52 1.725'" with points ls 1, \


# --- GRAPH c
AF = "./DiagramasEstadoFundamental/J3(-1.00)/AF_mZero-Gamma-J2.dat"
SAF = "./DiagramasEstadoFundamental/J3(-1.00)/SAF_mZero-Gamma-J2.dat"
AF_SAF = "./DiagramasEstadoFundamental/J3(-1.00)/Transition--AF-SAF.dat"

@TMARGIN; @RMARGIN; @XTICS; @NOYTICS
set label 1 'c)' @POS
set label 2 'J_{3} = -0.10' @POS_J3
set label 6 at 0.46, 2.1
set label 6 "" tc rgb '#8313E8'

plot AF using ($3<=0.51?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.51?$3:1/0):2 w l ls 3, \
    AF_SAF using 1:($2<=2.?$2:1/0) w l ls 3 dt 2, \

unset multiplot
### End multiplot

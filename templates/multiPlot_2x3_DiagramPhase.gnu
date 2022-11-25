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
set terminal pngcairo size 900,600 enhanced font 'Verdana,10'
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
set yrange [0:5.5]

# MACROS
# x- and ytics for each row resp. column
NOXTICS = "set xtics ('' 0, '' 0.2,'' 0.4,'' 0.6, '' 0.8, '' 1); \
          unset xlabel"
XTICS = "set xtics ('0' 0, '0.2' 0.2, '0.4' 0.4,'0.6' 0.6, '0.8' 0.8, '1' 1);\
          set xlabel 'J_2 / J_1'"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%.0f'; set ylabel 'T / J_1'"
# Placement of the a,b,c,d labels in the graphs
POS = "at graph 0.03,0.95 font ',10'"
Pos_h = "at graph 0.8,0.94 font ',10'"
POS_PM = "at graph 0.49,0.7 font ',12'"
POS_AF = "at graph 0.175,0.2 font ',12'"
POS_SAF = "at graph 0.76,0.2 font ',12'"
POS_SD = "at graph 0.48,0.13 font ',12'"

# Margins for each row resp. column
TMARGIN = "set tmargin at screen .96; set bmargin at screen 0.56"
BMARGIN = "set tmargin at screen 0.54; set bmargin at screen 0.14"
LMARGIN = "set lmargin at screen 0.05; set rmargin at screen 0.35"
CMARGIN = "set lmargin at screen 0.37; set rmargin at screen 0.67"
RMARGIN = "set lmargin at screen 0.69; set rmargin at screen 0.99"





### Start multiplot (2x3 layout)
set multiplot layout 2,3 rowsfirst

# --- GRAPH a
AF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma1/J3-0.05/AF_mZero-Temp-J2.dat"
SAF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma1/J3-0.05/SAF_mZero-Temp-J2.dat"
AF_SAF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma1/J3-0.05/Transition--AF-SAF.dat"

@NOXTICS; @YTICS; @TMARGIN; ; @LMARGIN
set label 1 'a)' @POS
set label 2 'h = 1' @Pos_h
set label 3 "PM" @POS_PM 
set label 4 "AF" @POS_AF 
set label 5 "SAF" @POS_SAF 
set label 6 at 0.535, 2.15
set label 6 "0.62" tc rgb '#8313E8'

plot AF using ($3<=0.51?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.62?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.62 && $3>=0.51?$3:1/0):2 w l ls 3 dt 2, \
    "<echo '.62 1.87'" with points ls 1, \
    AF_SAF using 1:($2<=1.3?$2:1/0) w l ls 3 dt 2



# --- GRAPH b
AF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma2/J3-0.05/AF_mZero-Temp-J2.dat"
SAF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma2/J3-0.05/SAF_mZero-Temp-J2.dat"

@NOXTICS; @NOYTICS; @TMARGIN; @CMARGIN 
set label 1 'b)' @POS
set label 2 'h = 2' @Pos_h
set label 6 at 0.325, .3
set label 6 "0.44" tc rgb '#8313E8'
set label 7 at 0.58, .3
set label 7 "0.57" tc rgb '#8313E8'
#set label 7 'SD' @POS_SD

plot AF using ($3<=0.5?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.5?$3:1/0):2 w l ls 3, \
    #AF_SD w l ls 3 dt 2, \
    #SD_PM using ($2<=0.6 && $2>=0.415?$2:1/0):1 w l ls 3 dt 4, \
    #SD_SAF using 2:($1<=1.1?$1:1/0) w l ls 3 dt 2,\
    #SD_SAF using 2:($1>=1.1 && $1<=1.40?$1:1/0) w l ls 3 dt 2


# --- GRAPH c
AF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma3/J3-0.05/AF_mZero-Temp-J2.dat"
SAF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma3/J3-0.05/SAF_mZero-Temp-J2.dat"
#AF_SD = "./DiagramaDeFase(J3-0.2)/Transition(J3-0_2)--AF-SD.dat"
#SD_SAF = "./DiagramaDeFase(J3-0.2)/Transition(J3-0_2)--SD-SAF.dat"
#SD_PM = "./DiagramaDeFase(J3-0.2)/Transition(J3-0_2)--SD-PM.dat"

@TMARGIN; @RMARGIN; @NOXTICS; @NOYTICS
set label 1 'c)' @POS
set label 2 'h = 3' @Pos_h
set label 4 at 0.07, 1
set label 5 at 0.85, 1
set label 6 at 0.1, .3
set label 6 "0.16" tc rgb '#8313E8'
set label 7 at 0.8, .3
set label 7 "0.84" tc rgb '#8313E8'

plot AF using ($3<=0.32?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.78?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.78 && $3>=0.72?$3:1/0):2 w l ls 3 dt 4, \
    #"<echo '.78 1.63'" with points ls 1, \
    #AF_SD w l ls 3 dt 4, \
    #SD_PM using ($1<=0.7 && $1>=0.3?$1:1/0):2 w l ls 3 dt 4, \
    #SD_SAF using 2:($1<=1.3?$1:1/0) w l ls 3 dt 2,\
    #SD_SAF using 2:($1>=1.1 && $1<=1.40?$1:1/0) w l ls 3 dt 4


# --- GRAPH d
AF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma1/J3-0.05/AF_mZero-Temp-J2.dat"
SAF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma1/J3-0.05/SAF_mZero-Temp-J2.dat"
AF_SAF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma1/J3-0.05/Transition--AF-SAF.dat"

@BMARGIN; @LMARGIN; @XTICS; @YTICS
set label 1 'd)' @POS
set label 2 'h = 1' @Pos_h
#set label 7 '' @POS_SD
set label 6 at 0.5, 2.35
set label 6 "0.59" tc rgb '#8313E8'

plot AF using ($3<=0.52?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.59?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.59 && $3>=0.52?$3:1/0):2 w l ls 3 dt 4, \
    "<echo '.59 2.09057'" with points ls 1, \
    AF_SAF w l ls 3 dt 4

# --- GRAPH e
AF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma2/J3-0.05/AF_mZero-Temp-J2.dat"
SAF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma2/J3-0.05/SAF_mZero-Temp-J2.dat"

@BMARGIN; @CMARGIN; @XTICS; @NOYTICS
set label 1 'e)' @POS
set label 2 'h = 2' @Pos_h
set label 6 at 0.465, 2.65
set label 6 "0.53" tc rgb '#8313E8'

plot AF using ($3<=0.518?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.53?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.53 && $3>=0.52?$3:1/0):2 w l ls 3 dt 4, \
    "<echo '.53 2.26435'" with points ls 1, \
    #AF_SAF w l ls 3 dt 4

# --- GRAPH f
AF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma3/J3-0.05/AF_mZero-Temp-J2.dat"
SAF = "./DiagramasQuanticos/J3_Ferromagnetico/Gamma3/J3-0.05/SAF_mZero-Temp-J2.dat"

@BMARGIN; @RMARGIN; @XTICS; @NOYTICS
set label 1 'f)' @POS
set label 2 'h = 3' @Pos_h
set label 6 ""

plot AF using ($3<=0.518?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.51?$3:1/0):2 w l ls 3, \
    AF_SAF using 1:($2<=2.9?$2:1/0) w l ls 3 dt 4

unset multiplot
### End multiplot

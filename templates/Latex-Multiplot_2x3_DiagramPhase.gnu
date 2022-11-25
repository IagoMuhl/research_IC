reset
set terminal lua tikz linewidth 2 standalone
set output 'dash.tex' 
set encoding utf8

#Set color line
set style line 1 lc rgb 'black' pt 7   # circle
set style line 2 lc rgb '#E8991E'
set style line 3 lc rgb '#070BE8'
set style line 5 lc rgb '#070BE8' dt 4
set style line 4 lc rgb '#8313E8'
unset key

# Enable the use of macros
set macros

set tics scale 0.5
#Determina o tamanho dos traços de marcação dentro do gráfico, quanto menor valor, menor o traço, em 0 não aparece nada, o gráfico fica "liso" por dentro.

set ytics 1
set xrange [0:1]
set yrange [0:5.5]

# MACROS
# x- and ytics for each row resp. column
NOXTICS = "set xtics ('' 0, '' 0.2,'' 0.4,'' 0.6, '' 0.8, '' 1); \
          unset xlabel "
XTICS = "set xtics ('$0$' 0, '$0.2$' 0.2, '$0.4$' 0.4,'$0.6$' 0.6, '$0.8$' 0.8, '$1$' 1) font ',10';\
          set xlabel '$J_{2} / J_{1} $' font ',10'"
NOYTICS = "set format y '' ; unset ylabel"
YTICS = "set ytics font ',9' offset 1,0; set format y '%.0f'; set ylabel '$T / J_{1}$' font ',10'"

#offset 1,0 Determina a distância do rótulo em relação ao eixo em questão
#font ',9' Exigimos uma fonte tamanho 9

# Placement of the a,b,c,d labels in the graphs
POS = "at graph 0.02,0.94 font ',8'"
POS_J3 = "at graph 0.58,0.94 font ',8'"
POS_PM = "at graph 0.41,0.7 font ',10'"
POS_AF = "at graph 0.106,0.15 font ',10'"
POS_SAF = "at graph 0.68,0.15 font ',10'"
POS_SD = "at graph 0.41,0.08 font ',10'"

# Margins for each row resp. column
TMARGIN = "set tmargin at screen .96; set bmargin at screen 0.56"
BMARGIN = "set tmargin at screen 0.54; set bmargin at screen 0.14"
LMARGIN = "set lmargin at screen 0.05; set rmargin at screen 0.35"
CMARGIN = "set lmargin at screen 0.37; set rmargin at screen 0.67"
RMARGIN = "set lmargin at screen 0.69; set rmargin at screen 0.99"





### Start multiplot (2x2 layout)
set multiplot layout 2,3 rowsfirst

# --- GRAPH a
AF = "./DiagramasDeFase(J3-0.0)/AF(J3-0_0)_mZero-Temp-J2-Z.dat"
SAF = "./DiagramasDeFase(J3-0.0)/SAF(J3-0_0)_mZero-Temp-J2-Z.dat"
AF_SAF = "./DiagramasDeFase(J3-0.0)/Transition(J3-0_0)--AF-SAF.dat"

@NOXTICS; @YTICS; @TMARGIN; ; @LMARGIN
set label 1 '$a)$' @POS
set label 2 '$\phantom{-}J_{3} = 0.0$' @POS_J3
set label 3 "$PM$" @POS_PM
set label 4 "$AF$" @POS_AF
set label 5 "$SAF$" @POS_SAF
set label 6 at 0.49, 2.3
set label 6 "$0.66$" font ',9' tc rgb '#8313E8'

plot AF using ($3<=0.51?$3:1/0):2 w l ls 3, \
    SAF using ($3>0.66?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.66 && $3>=0.5?$3:1/0):2 w l ls 5, \
    AF_SAF using 2:($1<=0.7?$1:1/0) w l ls 5, \
    "<echo '.66 1.96'" with points ls 1
    



# --- GRAPH b
AF = "./DiagramaDeFase(J3-0.1)/AF_mZero-Temp-J2-Z.dat"
SAF = "./DiagramaDeFase(J3-0.1)/SAF_mZero-Temp-J2-Z.dat"
AF_SD = "./DiagramaDeFase(J3-0.1)/Transition--AF-SD.dat"
SD_SAF = "./DiagramaDeFase(J3-0.1)/Transition--SD-SAF.dat"
SD_PM = "./DiagramaDeFase(J3-0.1)/Transition--SD-PM.dat"

@NOXTICS; @NOYTICS; @TMARGIN; @CMARGIN 
set label 1 '$b)$' @POS
set label 2 '$\phantom{-}J_{3} = 0.1$' @POS_J3
set label 6 at 0.55, 2.15
set label 6 "$0.72$" tc rgb '#8313E8'
set label 7 '$SD$' @POS_SD

plot AF using ($3<=0.423?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.72?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.72 && $3>=0.635?$3:1/0):2 w l ls 5, \
    AF_SD w l ls 5, \
    SD_PM using ($2<=0.6 && $2>=0.415?$2:1/0):1 w l ls 5, \
    SD_SAF using 2:($1<=1.1?$1:1/0) w l ls 5,\
    SD_SAF using 2:($1>=1.1 && $1<=1.40?$1:1/0) w l ls 5, \
    "<echo '.72 1.82'" with points ls 1


# --- GRAPH c
AF = "./DiagramaDeFase(J3-0.2)/AF_mZero-Temp-J2-Z.dat"
SAF = "./DiagramaDeFase(J3-0.2)/SAF_mZero-Temp-J2-Z.dat"
AF_SD = "./DiagramaDeFase(J3-0.2)/Transition(J3-0_2)--AF-SD.dat"
SD_SAF = "./DiagramaDeFase(J3-0.2)/Transition(J3-0_2)--SD-SAF.dat"
SD_PM = "./DiagramaDeFase(J3-0.2)/Transition(J3-0_2)--SD-PM.dat"

@TMARGIN; @RMARGIN; @NOXTICS; @NOYTICS
set label 1 '$c)$' @POS
set label 2 '$\phantom{-}J_{3} = 0.2$' @POS_J3
set label 6 at 0.605, 1.95
set label 6 "$0.78$" tc rgb '#8313E8'

plot AF using ($3<=0.32?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.78?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.78 && $3>=0.72?$3:1/0):2 w l ls 5, \
    AF_SD w l ls 5, \
    SD_PM using ($1<=0.7 && $1>=0.3?$1:1/0):2 w l ls 5, \
    SD_SAF using 2:($1<=1.3?$1:1/0) w l ls 5,\
    "<echo '.78 1.63'" with points ls 1


# --- GRAPH d
AF = "./DiagramaDeFase(J3--0.1)/AF_mZero-Temp-J2-Z.dat"
SAF = "./DiagramaDeFase(J3--0.1)/SAF_mZero-Temp-J2-Z.dat"
AF_SAF = "./DiagramaDeFase(J3--0.1)/Transition--AF-SAF.dat"

@BMARGIN; @LMARGIN; @XTICS; @YTICS
set label 1 '$d)$' @POS
set label 2 '$J_{3} = -0.1$' @POS_J3
set label 7 '' @POS_SD
set label 6 at 0.43, 2.4
set label 6 "$0.59$" tc rgb '#8313E8'

plot AF using ($3<=0.52?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.59?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.59 && $3>=0.52?$3:1/0):2 w l ls 5, \
    AF_SAF w l ls 5, \
    "<echo '.59 2.09057'" with points ls 1



# --- GRAPH e
AF = "./DiagramaDeFase(J3--0.2)/AF_mZero-Temp-J2-Z.dat"
SAF = "./DiagramaDeFase(J3--0.2)/SAF_mZero-Temp-J2-Z.dat"
AF_SAF = "./DiagramaDeFase(J3--0.2)/Transition--AF-SAF.dat"

@BMARGIN; @CMARGIN; @XTICS; @NOYTICS
set label 1 '$e)$' @POS
set label 2 '$J_{3} = -0.2$' @POS_J3
set label 6 at 0.405, 2.8
set label 6 "$0.53$" tc rgb '#8313E8'

plot AF using ($3<=0.518?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.53?$3:1/0):2 w l ls 3, \
    SAF using ($3<=0.53 && $3>=0.52?$3:1/0):2 w l ls 5, \
    AF_SAF w l ls 5, \
    "<echo '.53 2.26435'" with points ls 1

# --- GRAPH f
AF = "./DiagramaDeFase(J3--0.3)/AF_mZero-Temp-J2-Z.dat"
SAF = "./DiagramaDeFase(J3--0.3)/SAF_mZero-Temp-J2-Z.dat"
AF_SAF = "./DiagramaDeFase(J3--0.3)/Transition--AF-SAF.dat"

@BMARGIN; @RMARGIN; @XTICS; @NOYTICS
set label 1 '$f)$' @POS
set label 2 '$J_{3} = -0.3$' @POS_J3
set label 6 ""

plot AF using ($3<=0.518?$3:1/0):2 w l ls 3, \
    SAF using ($3>=0.51?$3:1/0):2 w l ls 3, \
    AF_SAF using 1:($2<=2.9?$2:1/0) w l ls 5

unset multiplot
### End multiplot


### Transformando em EPS a partir do PDF
#pdftops -eps dash.pdf

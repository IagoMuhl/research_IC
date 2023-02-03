reset
set terminal lua tikz linewidth 3 standalone
set output 'dash.tex' 
set encoding utf8

#Set color line
set style line 1 lc rgb 'black' pt 7   # circle
set style line 2 lc rgb '#E8991E'
set style line 3 lc rgb '#1E90FF'
set style line 5 lc rgb '#1E90FF' dt 4
set style line 4 lc rgb '#0000CD'
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
AF = "J3=0.20/AF_J2-T.dat"
SAF = "J3=0.20/SAF_J2-T.dat"
AF_SD = "J3=0.20/AF-SD(J2-T).dat"
SAF_SD = "J3=0.20/SAF-SD(J2-T).dat"
SD = "J3=0.20/SD_J2-T.dat"

@NOXTICS; @YTICS; @TMARGIN; ; @LMARGIN
set label 1 '$a)$' @POS
set label 2 '$\phantom{-}J_{3} = 0.2$' @POS_J3
set label 3 "$PM$" @POS_PM
set label 3 "$SD$" @POS_SD
set label 4 "$AF$" @POS_AF
set label 5 "$SAF$" @POS_SAF
set label 6 at 0.6, 2.
set label 6 "$0.78$" font ',9' tc rgb '#1C1C1C'

plot AF using ($1<=0.32?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.78?$1:2/0):2 w l ls 3, \
    SAF using ($1<=0.78 && $1>=0.6?$1:2/0):2 w l ls 3 dt 4, \
    "<echo '.78 1.61763'" with points ls 1, \
    SD using ($1>=0.32?$1:2/0):2 w l ls 3 dt 4,\
    AF_SD using 1:2 w l ls 3 dt 4,\
    SAF_SD using ($1<=0.72?$1:2/0):2 w l ls 3 dt 4,\
    



# --- GRAPH b
AF = "J3=0.10/AF_J2-T.dat"
SAF = "J3=0.10/SAF_J2-T.dat"
SD = "J3=0.10/SD_J2-T.dat"
AF_SD = "J3=0.10/AF-SD(J2-T).dat"
SAF_SD = "J3=0.10/SAF-SD(J2-T).dat"

@NOXTICS; @NOYTICS; @TMARGIN; @CMARGIN 
set label 1 '$b)$' @POS
set label 2 '$\phantom{-}J_{3} = 0.1$' @POS_J3
set label 6 at 0.55, 2.15
set label 6 "$0.72$" tc rgb '#1C1C1C'


plot AF using ($1<=0.43?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.72?$1:2/0):2 w l ls 3, \
    SAF using ($1<=0.72 && $1>=0.6?$1:2/0):2 w l ls 3 dt 4, \
    AF_SD using 1:2 w l ls 3 dt 4, \
    SD using ($1<=0.6 && $1>=0.415?$1:2/0):2 w l ls 3 dt 4, \
    SAF_SD using 1:2 w l ls 3 dt 4,\
    "<echo '.72 1.79470'" with points ls 1, \


# --- GRAPH c
AF = "J3=0.00/AF_J2-T.dat"
SAF = "J3=0.00/SAF_J2-T.dat"
AF_SAF = "J3=0.00/J2_fix-T.dat"

@TMARGIN; @RMARGIN; @NOXTICS; @NOYTICS
set label 1 '$c)$' @POS
set label 2 '$\phantom{-}J_{3} = 0.0$' @POS_J3
set label 3 "$$" @POS_SD
set label 6 at 0.48, 2.35
set label 6 "$0.66$" tc rgb '#1C1C1C'


plot AF using ($1<=0.52?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.66?$1:2/0):2 w l ls 3, \
    SAF using ($1<=0.78 && $1>=0.72?$1:2/0):2 w l ls 3 dt 4, \
    AF_SAF using 1:2 w l ls 3 dt 4, \
    "<echo '.655 1.94347'" with points ls 1, \


# --- GRAPH d
AF = "J3=-0.10/AF_J2-T.dat"
SAF = "J3=-0.10/SAF_J2-T.dat"
AF_SAF = "J3=-0.10/J2_fix-T.dat"

@BMARGIN; @LMARGIN; @XTICS; @YTICS
set label 1 '$d)$' @POS
set label 2 '$J_{3} = -0.1$' @POS_J3
set label 7 '' @POS_SD
set label 6 at 0.4, 2.45
set label 6 "$0.59$" tc rgb '#1C1C1C'

plot AF using ($1<=0.52?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.59?$1:2/0):2 w l ls 3, \
    SAF using ($1<=0.59 && $1>=0.52?$1:2/0):2 w l ls 3 dt 4, \
    AF_SAF w l ls 3 dt 4, \
    "<echo '.595 2.11885'" with points ls 1



# --- GRAPH e
AF = "J3=-0.20/AF_J2-T.dat"
SAF = "J3=-0.20/SAF_J2-T.dat"
AF_SAF = "J3=-0.20/J2_fix-T.dat"

@BMARGIN; @CMARGIN; @XTICS; @NOYTICS
set label 1 '$e)$' @POS
set label 2 '$J_{3} = -0.2$' @POS_J3
set label 6 at 0.405, 2.8
set label 6 "$0.53$" tc rgb '#1C1C1C'

plot AF using ($1<=0.518?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.53?$1:2/0):2 w l ls 3, \
    SAF using ($1<=0.53 && $1>=0.52?$1:2/0):2 w l ls 3 dt 4, \
    AF_SAF using 1:2 w l ls 3 dt 4, \
    "<echo '.53 2.26435'" with points ls 1
    

# --- GRAPH f
AF = "J3=-0.30/AF_J2-T.dat"
SAF = "J3=-0.30/SAF_J2-T.dat"
AF_SAF = "J3=-0.30/J2_fix-T.dat"

@BMARGIN; @RMARGIN; @XTICS; @NOYTICS
set label 1 '$f)$' @POS
set label 2 '$J_{3} = -0.3$' @POS_J3
set label 6 ""

plot AF using ($1<=0.518?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.51?$1:2/0):2 w l ls 3, \
    AF_SAF using 1:($2<=2.9?$2:1/0) w l ls 3 dt 4

unset multiplot
### End multiplot


### Transformando em EPS a partir do PDF
#pdftops -eps dash.pdf

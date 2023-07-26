reset
set terminal lua tikz linewidth 2 standalone
set output 'dash.tex' 
set encoding utf8

#Set color line
set style line 1 lc rgb 'black' pt 7   # circle
set style line 2 lc rgb '#E8991E'
set style line 3 lc rgb 'black'
set style line 5 lc rgb '#1E90FF' dt 4
set style line 4 lc rgb '#0000CD'
unset key

# Enable the use of macros
set macros

set tics scale 0.5
#Determina o tamanho dos traços de marcação dentro do gráfico, quanto menor valor, menor o traço, em 0 não aparece nada, o gráfico fica "liso" por dentro.

set ytics 1
set xrange [0:1]
set yrange [0:4.5]

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
POS_J3 = "at graph 0.64,0.92 font ',8'"
POS_PM = "at graph 0.39,0.65 font ',10'"
POS_AF = "at graph 0.06,0.15 font ',10'"
POS_SAF = "at graph 0.70,0.15 font ',10'"
POS_SD = "at graph 0.4,0.06 font ',10'"

# Margins for each row resp. column
TMARGIN = "set tmargin at screen .96; set bmargin at screen 0.56"
BMARGIN = "set tmargin at screen 0.54; set bmargin at screen 0.14"
LMARGIN = "set lmargin at screen 0.05; set rmargin at screen 0.35"
CMARGIN = "set lmargin at screen 0.37; set rmargin at screen 0.67"
RMARGIN = "set lmargin at screen 0.69; set rmargin at screen 0.99"





### Start multiplot (2x2 layout)
set multiplot layout 2,3 rowsfirst

# --- GRAPH a
AF = "Quantico-Math/J3=0.2|Gamma/Gamma = 1/AF-PM(J2-T).dat"
SAF = "Quantico-Math/J3=0.2|Gamma/Gamma = 1/SAF-PM(J2-T).dat"
AF_SD = "Quantico-Math/J3=0.2|Gamma/Gamma = 1/AF-SD(J2-T).dat"
SAF_SD = "Quantico-Math/J3=0.2|Gamma/Gamma = 1/SD-SAF(J2-T).dat"
SD = "Quantico-Math/J3=0.2|Gamma/Gamma = 1/SD-PM(J2-T).dat"

@NOXTICS; @YTICS; @TMARGIN; ; @LMARGIN
set label 1 '$a)$' @POS
set label 2 '$\phantom{-}h = 1$' @POS_J3
set label 7 "$PM$" @POS_PM
set label 3 "$SD$" @POS_SD
set label 4 "$AF$" @POS_AF
set label 5 "$SAF$" @POS_SAF
set label 6 at 0.57, 1.6
set label 6 "$0.77$" font ',9' tc rgb '#1C1C1C'

plot AF using ($1<=0.32?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.77?$1:2/0):2 w l ls 3, \
    SAF using ($1<=0.77 && $1>=0.6?$1:2/0):2 w l ls 3 dt 4, \
    "<echo '0.77    1.37376 '" with points ls 1, \
    SD using 1:2 w l ls 3 dt 4,\
    AF_SD using 1:2 w l ls 3 dt 4,\
    SAF_SD using ($1<=0.72?$1:2/0):2 w l ls 3 dt 4,\
    



# --- GRAPH b
AF = "Quantico-Math/J3=0.2|Gamma/Gamma = 2/AF-PM(J2-T).dat"
SAF = "Quantico-Math/J3=0.2|Gamma/Gamma = 2/SAF-PM(J2-T).dat"

POS_AF = "at graph 0.18,0.25 font ',10'"
POS_SAF = "at graph 0.53,0.25 font ',10'"

@NOXTICS; @NOYTICS; @TMARGIN; @CMARGIN 
set label 1 '$b)$' @POS
set label 2 '$\phantom{-}h = 2$' @POS_J3
set label 4 "$AF$" @POS_AF
set label 5 "$SAF$" @POS_SAF
set label 3 "$$" @POS_SD
set label 6 at 0.55, 2.15
set label 6 "$$" tc rgb '#1C1C1C'
set label 7 "$PM$" @POS_PM


plot AF using 1:2 w l ls 3, \
    SAF using 1:2 w l ls 3, \
    #"<echo '.72 1.79470'" with points ls 1, \


# --- GRAPH c
AF = "Quantico-Math/J3=0.2|Gamma/Gamma = 3/zero.dat"
SAF = "Quantico-Math/J3=0.2|Gamma/Gamma = 3/zero.dat"

POS_AF = "at graph 0.065,0.23 font ',10'"
POS_SAF = "at graph 0.65,0.23 font ',10'"

@TMARGIN; @RMARGIN; @NOXTICS; @NOYTICS
set label 1 '$c)$' @POS
set label 2 '$\phantom{-}h = 3$' @POS_J3
set label 3 "$$" @POS_SD
set label 4 "$$" @POS_AF
set label 5 "$$" @POS_SAF
set label 7 "$PM$" @POS_PM



plot AF using ($1<=0.52?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.66?$1:2/0):2 w l ls 3, \
    SAF using ($1<=0.78 && $1>=0.72?$1:2/0):2 w l ls 3 dt 4, \



# --- GRAPH d
AF = "Quantico-Math/J3=-0.2|Gamma/Gamma = 1/AF-PM(J2-T).dat"
SAF = "Quantico-Math/J3=-0.2|Gamma/Gamma = 1/SAF-PM(J2-T).dat"
AF_SAF = "Quantico-Math/J3=-0.2|Gamma/Gamma = 1/AF-SAF(J2-T).dat"


POS_AF = "at graph 0.1,0.23 font ',10'"
POS_SAF = "at graph 0.66,0.23 font ',10'"


@BMARGIN; @LMARGIN; @XTICS; @YTICS
set label 1 '$d)$' @POS
set label 2 '$\phantom{-}h = 1$' @POS_J3
set label 4 "$AF$" @POS_AF
set label 5 "$SAF$" @POS_SAF
set label 7 '' @POS_SD
set label 6 at 0.52, 1.95
set label 6 "$0.52$" tc rgb '#1C1C1C'
set label 7 "$PM$" @POS_PM

plot AF using ($1<=0.52?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.52?$1:2/0):2 w l ls 3, \
    SAF using ($1<=0.52 && $1>=0.52?$1:2/0):2 w l ls 3 dt 4, \
    AF_SAF w l ls 3 dt 4, \
    "<echo '0.52    2.07211'" with points ls 1



# --- GRAPH e
AF = "Quantico-Math/J3=-0.2|Gamma/Gamma = 2/AF-PM(J2-T).dat"
SAF = "Quantico-Math/J3=-0.2|Gamma/Gamma = 2/SAF-PM(J2-T).dat"
AF_SAF = "Quantico-Math/J3=-0.2|Gamma/Gamma = 2/AF-SAF(J2-T).dat"

@BMARGIN; @CMARGIN; @XTICS; @NOYTICS
set label 1 '$e)$' @POS
set label 2 '$\phantom{-}h = 2$' @POS_J3
set label 6 "$$" tc rgb '#1C1C1C'
set label 7 "$PM$" @POS_PM

plot AF using ($1<=0.51?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.51?$1:2/0):2 w l ls 3, \
    AF_SAF using 1:2 w l ls 3 dt 4

    

# --- GRAPH f
AF = "Quantico-Math/J3=-0.2|Gamma/Gamma = 3/AF-PM(J2-T).dat"
SAF = "Quantico-Math/J3=-0.2|Gamma/Gamma = 3/SAF-PM(J2-T).dat"


@BMARGIN; @RMARGIN; @XTICS; @NOYTICS
set label 1 '$f)$' @POS
set label 2 '$\phantom{-}h = 3$' @POS_J3
set label 4 "$AF$" @POS_AF
set label 5 "$SAF$" @POS_SAF
set label 6 ""
set label 7 "$PM$" @POS_PM

plot AF using 1:2 w l ls 3, \
    SAF using 1:2 w l ls 3, \

unset multiplot
### End multiplot


### Transformando em EPS a partir do PDF
#pdftops -eps dash.pdf

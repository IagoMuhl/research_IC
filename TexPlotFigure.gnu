reset
set terminal lua tikz linewidth 3 standalone
set output 'NomeDoArquivo.tex' 
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

#set ytics scale 0.2
#set xtics scale 0.1
#Determina o tamanho dos traços de marcação dentro do gráfico, quanto menor valor, menor o traço, em 0 não aparece nada, o gráfico fica "liso" por dentro.

#set ytics 1

set xrange [0.3:0.8]
set yrange [0:1]

# MACROS
# x- and ytics for each row resp. column

#offset 1,0 Determina a distância do rótulo em relação ao eixo em questão
#font ',9' Exigimos uma fonte tamanho 9

# Placement of the a,b,c,d labels in the graphs
# Posição da escrita na figura
#POS = "at graph 0.02,0.94 font ',8'"
#POS_J3 = "at graph 0.58,0.94 font ',8'"
#POS_PM = "at graph 0.41,0.7 font ',10'"
#POS_AF = "at graph 0.106,0.15 font ',10'"

POS_SAF = "at graph 0.77,0.9 font ',10'"
POS_SD = "at graph 0.78,0.83 font ',10'"
POS_Negative02 = "at graph 0.29, 0.1 font ',10'"
POS_Negative01 = "at graph 0.5, 0.5 font ',10'"
POS_None = "at graph 0.625, 0.5 font ',10'"
POS_Positive01 = "at graph 0.752, 0.5 font ',10'"
POS_Positive02 = "at graph 0.88, 0.5 font ',10'"


### Start multiplot (2x2 layout)
# 
#set multiplot layout 2,3 rowsfirst

# --- GRAPH 
# Leitura dos dados
Negative02 = '-0.2_Jump_J2-T-m.dat'
Negative01 = '-0.1_Jump_J2-T-m.dat'
None = '0.0_Jump_J2-T-m.dat'
Positive01 = '0.1_Jump_J2-T-m.dat'
Positive02 = '0.2_Jump_J2-T-m.dat'




set label 1 "$SAF$" @POS_SAF
set label 2 "$SD$" @POS_SD 
set label 3 "$J_{3}=-0.2$" @POS_Negative02
set label 4 "$J_{3}=-0.1$" @POS_Negative01 rotate by -60
set label 5 "$J_{3}=0.0$" @POS_None rotate by -60
set label 6 "$J_{3}=0.1$" @POS_Positive01 rotate by -60
set label 7 "$J_{3}=0.2$" @POS_Positive02 rotate by -60
#set label 1 at 0.7, 0.9
#set label 1 "SAF" tc rgb 'black'

#set label 2 at 0.708, 0.8
#set label 2 "SD" tc rgb 'black'




plot Positive02 using ($1>=0.7?$1:3/0):3 w l ls 1, \
    Positive02 using ($1<=0.7?$1:3/0):3 w l ls 3, \
    Positive01 using ($1>=0.6?$1:3/0):3 w l ls 1, \
    None using ($1>=0.52?$1:3/0):3 w l ls 1, \
    Negative01 using ($1>=0.52?$1:3/0):3 w l ls 1, \
    Negative02 using ($1>=0.52?$1:3/0):3 w l ls 1, \
    "<echo '0.75, 0.9'" with points ls 1, \
    "<echo '0.75, 0.83'" with points ls 3 \
    

#unset multiplot
### End multiplot


### Transformando em EPS a partir do PDF
#pdftops -eps NomeDoArquivo.pdf

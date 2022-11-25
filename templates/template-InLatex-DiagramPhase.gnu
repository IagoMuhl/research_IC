reset
set terminal lua tikz linewidth 3 standalone
set output 'dash.tex' 
set encoding utf8
#set output 'DiagramaDeFaseJ3-0_0.png'
set xlabel "$J_{2} / |J_{1}|$"
set ylabel "$T / |J_{1}|$"
AF = "./AF(J3-0_0)_mZero-Temp-J2-Z.dat"
SAF = "./SAF(J3-0_0)_mZero-Temp-J2-Z.dat"
#AF_SD = "./Transition--AF-SD.dat"
#SD_SAF = "./Transition--SD-SAF.dat"
#SD_PM = "./Transition--SD-PM.dat"
#set title '$J_{3} = +0.0$'
set yrange [0:4]
set xrange [0:1]

#set xtics (0.5, 0.66)

#Style
set style line 1 lc rgb 'black' pt 7   # circle
set style line 2 lc rgb '#E8991E'
set style line 3 lc rgb '#070BE8'
set style line 4 lc rgb '#8313E8'

#Label
set label 1 at 0.48, 2.5
set label 1 "PM" tc rgb '#E8991E'

set label 2 at 0.2, 1.
set label 2 "AF" tc rgb '#E8991E'

set label 3 at 0.78, 1.
set label 3 "SAF" tc rgb '#E8991E'

#set label 5 at 0.48, .5
#set label 5 "SD" tc lt 4

set label 4 at 0.6, 2.15
set label 4 "$0.66$" tc rgb '#8313E8'

#Legenda
set key off

set label 5 at 0.0, 3.78
set label 5 "a)" tc rgb 'black'

set label 6 at 0.83, 3.78
set label 6 "$J_{3} = 0.0$" tc rgb 'black'

#Plot
plot AF using ($3<=0.51?$3:1/0):2 w l ls 3, \
    SAF using ($3>0.66?$3:1/0):2 w l ls 4, \
    SAF using ($3<0.66 && $3>=0.0?$3:1/0):2 w l ls 4 dt 4, \
    "<echo '.66 1.95'" with points ls 1, \
    #AF_SD w l lt 2 dt 2, \
    #SD_PM using ($1<0.7 && $1>=0.2?$1:1/0):2 w l lt 2 dt 2, \
    #SD_SAF using 2:($1<1.6?$1:1/0) w l lt 2 dt 2
    
#Lógica if-else
#$3<=0.5?  
#Se o elem. da coluna 3 for menor ou igual a 0.5
#$3:1/0
#Plote a coluna 3, senão não plote nada (NaN=1/0)

#Operadores
#O Operador .or. no gnuplot é ||
#O Operador .and. no gnuplot é &&
#

reset
set terminal pngcairo dashed enhanced size 600,400 font 'arial,12' fontscale 1.0
set encoding utf8
set output 'DiagramaDimmer(J3=0.2).png'
set xlabel "J_{2} / J_{1}"
set ylabel "J_{3} / J_{1}"

#AF = "./Gamma1/J3(0.05)/AF_mZero-Temp-J2_J3(0.05).dat"
#SAF = "./Gamma1/J3(0.05)/SAF_mZero-Temp-J2_J3(0.05).dat"
#SD_PM = "./Gamma1/J3(0.05)/Transition--SD-PM_J3(0.05).dat"
#AF_SD = "./Gamma1/J3(0.05)/Transition--AF-SD_J3(0.05).dat"
#SD_SAF = "./Gamma1/J3(0.05)/Transition--SD-SAF_J3(0.05).dat"
#set title '$J_{3} = +0.0$'
AF = "./AF_J2-T.dat"
SAF = "./SAF_J2-T.dat"
SD = "./SD_J2-T.dat"
AF_SD = "AF-SD(J2-T).dat"
SAF_SD = "SAF-SD(J2-T).dat"

set yrange [0:3]
set xrange [0:1]

#set xtics (0.5, 0.66)

#Style
set style line 1 lc rgb 'black' pt 7   # circle
set style line 2 lc rgb '#E8991E'
set style line 3 lc rgb '#070BE8'
set style line 4 lc rgb '#8313E8'

#Label
set label 1 at 0.48, 2
set label 1 "PM" tc rgb '#8313E8'

set label 2 at 0.2, 1
set label 2 "AF" tc rgb '#8313E8'

set label 3 at 0.78, 1
set label 3 "SAF" tc rgb '#8313E8'

set label 5 at 0.484, .4
set label 5 "SD" tc rgb '#8313E8'

#Ponto Crítico
set label 4 at 0.7, 1.76
set label 4 "0.78" tc rgb '#8313E8'

#set label 9 at 0.65, .2
#set label 9 "0.64" tc rgb '#8313E8'

#Legenda
set key off

#set label 6 at 0.0, 4.78
#set label 6 "b)" tc rgb 'black'

set label 7 at 0.82, 3.2
set label 7 "J_{3} = 0.20" tc rgb 'black'

#set label 8 at 0.7, 4
#set label 8 "h = 1" tc rgb 'black'

#Plot
#plot AF_SAF using 1:2 w l ls 3, \
    #AF_SAF using 3:4 w l ls 3, \
   # AF using ($1<=0.518?$1:2/0):2 w l ls 3, \
    #SAF using ($1>=0.51?$1:2/0):2 w l ls 3, \
    #SAF using ($1<=0.5 && $1>=0.5?$1:2/0):2 w l ls 3 dt 3, \
    #AF_SAF using 1:($2<=2.7?$2:1/0) w l ls 3 dt 2, \
    #"<echo '.525 2.23617'" with points ls 1, \
    
plot AF using ($1<=0.32?$1:2/0):2 w l ls 3, \
    SAF using ($1>=0.78?$1:2/0):2 w l ls 3, \
    SD using ($1>=0.32?$1:2/0):2 w l ls 3 dt 2, \
    AF_SD using ($1<=0.43?$1:2/0):2 w l ls 3 dt 2, \
    SAF_SD using 1:2 w l ls 3 dt 2, \
    SAF using ($1<=0.78 && $1>=0.6?$1:2/0):2 w l ls 3 dt 2, \
    "<echo '.78 1.61763'" with points ls 1, \
    
    #SD_PM using ($1>=0.48 && $1<=0.55?$1:1/0):2 w l ls 3 dt 2, \
    #SD_SAF using 1:($2<=0.859?$2:1/0) w l ls 3 dt 2, \
    #AF_SD using 1:($2<=.751?$2:1/0) w l ls 3 dt 2, \
    
    
#Lógica if-else
#$3<=0.5?  
#Se o elem. da coluna 3 for menor ou igual a 0.5
#$3:1/0
#Plote a coluna 3, senão não plote nada (NaN=1/0)

#Operadores
#O Operador .or. no gnuplot é ||
#O Operador .and. no gnuplot é &&
#

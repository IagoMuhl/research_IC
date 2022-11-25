reset
set terminal pngcairo dashed enhanced size 600,400 font 'arial,12' fontscale 1.0
set encoding utf8
set output 'DiagramaClassico-J3(0.00).png'
set xlabel "J_{2} / J_{1}"
set ylabel "T / J_{1}"

#AF = "./Gamma1/J3(0.05)/AF_mZero-Temp-J2_J3(0.05).dat"
#SAF = "./Gamma1/J3(0.05)/SAF_mZero-Temp-J2_J3(0.05).dat"
#SD_PM = "./Gamma1/J3(0.05)/Transition--SD-PM_J3(0.05).dat"
#AF_SD = "./Gamma1/J3(0.05)/Transition--AF-SD_J3(0.05).dat"
#SD_SAF = "./Gamma1/J3(0.05)/Transition--SD-SAF_J3(0.05).dat"
#set title '$J_{3} = +0.0$'
AF = "./AF_T-J2.dat"
SAF = "./SAF_T-J2.dat"
AF_SAF = "./J2_fix-T.dat"

set yrange [0:5]
set xrange [0:1]

#set xtics (0.5, 0.66)

#Style
set style line 1 lc rgb 'black' pt 7   # circle
set style line 2 lc rgb '#E8991E'
set style line 3 lc rgb '#070BE8'
set style line 4 lc rgb '#8313E8'

#Label
set label 1 at 0.48, 4
set label 1 "PM" tc rgb '#8313E8'

set label 2 at 0.05, 1.
set label 2 "AF" tc rgb '#8313E8'

set label 3 at 0.9, 1.
set label 3 "SAF" tc rgb '#8313E8'

#set label 5 at 0.484, .4
#set label 5 "SD" tc rgb '#8313E8'

#Ponto Crítico
#set label 4 at 0.625, 1.85
#set label 4 "0.67" tc rgb '#8313E8'

#set label 9 at 0.65, .2
#set label 9 "0.64" tc rgb '#8313E8'

#Legenda
set key off

#set label 6 at 0.0, 4.78
#set label 6 "b)" tc rgb 'black'

set label 7 at 0.2, 4
set label 7 "J_{3} = 0.00" tc rgb 'black'

#set label 8 at 0.7, 4
#set label 8 "h = 1" tc rgb 'black'

#Plot
plot AF using ($2<=0.52?$2:1/0):1 w l ls 3, \
    SAF using ($2>=0.57?$2:1/0):1 w l ls 3, \
    AF_SAF using 1:($2<=1.5?$2:1/0) w l ls 3 dt 2, \
    "<echo '.675 1.63572'" with points ls 1#, \
    #SAF using ($3<=0.675 && $3>=0.53?$3:1/0):2 w l ls 3 dt 2, \
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

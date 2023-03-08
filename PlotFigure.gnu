reset
set terminal pngcairo dashed enhanced size 600,400 linewidth 2 font 'arial,12' fontscale 1.0
set encoding utf8
set output 'MagJump.png'
set xlabel "J_{2} / J_{1}"
set ylabel "{/Symbol D}m" rotate by 0


#AF = "./Gamma1/J3(0.05)/AF_mZero-Temp-J2_J3(0.05).dat"
#SAF = "./Gamma1/J3(0.05)/SAF_mZero-Temp-J2_J3(0.05).dat"
#SD_PM = "./Gamma1/J3(0.05)/Transition--SD-PM_J3(0.05).dat"
#AF_SD = "./Gamma1/J3(0.05)/Transition--AF-SD_J3(0.05).dat"
#SD_SAF = "./Gamma1/J3(0.05)/Transition--SD-SAF_J3(0.05).dat"
#set title '$J_{3} = +0.0$'

Negative02 = '-0.2_Jump_J2-T-m.dat'
Negative01 = '-0.1_Jump_J2-T-m.dat'
None = '0.0_Jump_J2-T-m.dat'
Positive01 = '0.1_Jump_J2-T-m.dat'
Positive02 = '0.2_Jump_J2-T-m.dat'

set yrange [0:1]
set xrange [0.3:0.8]

#set xtics (0.5, 0.66)

#Style
set style line 1 lc rgb 'black' pt 7   # circle
set style line 2 lc rgb '#E8991E'
set style line 3 lc rgb '#428af5'
set style line 4 lc rgb '#8313E8'

set label 1 at 0.7, 0.9
set label 1 "SAF" tc rgb 'black'

set label 2 at 0.708, 0.8
set label 2 "SD" tc rgb 'black'

set label 3 at 0.475, 0.1
set label 3 "J_3=-0.2" 
set label 3 font 'reboto,10'

set label 4 at 0.552, 0.5 
set label 4 "J_3=-0.1" rotate by -60
set label 4 font 'reboto,10'

set label 5 at 0.62, 0.5 
set label 5 "J_3=0.0" rotate by -60
set label 5 font 'reboto,10'

set label 6 at 0.683, 0.5 
set label 6 "J_3=0.1" rotate by -60
set label 6 font 'reboto,10'

set label 7 at 0.746, 0.5 
set label 7 "J_3=0.2" rotate by -60
set label 7 font 'reboto,10'

#Label
#set label 1 at 0.6, 2.8
#set label 1 "PM" tc rgb '#070BE8'

#set label 2 at 0.2, 1
#set label 2 "AF" tc rgb '#070BE8'

#set label 3 at 0.78, 1
#set label 3 "SAF" tc rgb '#070BE8'

#set label 5 at 0.484, .4
#set label 5 "SD" tc rgb '#8313E8'

#Ponto Crítico
#set label 4 at 0.6, 1.9
#set label 4 "0.59" tc rgb '#8313E8'

#set label 9 at 0.65, .2
#set label 9 "0.64" tc rgb '#8313E8'

#Legenda
set key off

#set label 6 at 0.0, 4.78
#set label 6 "b)" tc rgb 'black'

#set label 7 at 0.88, 3.2
#set label 7 "J_{3} = -0.10" tc rgb 'black'

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
    

plot Positive02 using ($1>=0.7?$1:3/0):3 w l ls 1, \
    Positive02 using ($1<=0.7?$1:3/0):3 w l ls 3, \
    Positive01 using ($1>=0.6?$1:3/0):3 w l ls 1, \
    None using ($1>=0.52?$1:3/0):3 w l ls 1, \
    Negative01 using ($1>=0.52?$1:3/0):3 w l ls 1, \
    Negative02 using ($1>=0.52?$1:3/0):3 w l ls 1, \
    "<echo '0.75, 0.9'" with points ls 1, \
    "<echo '0.75, 0.8'" with points ls 3 \
    #M_SAF using ($1<=0.59?$1:2/0):2 w l ls 3 dt 4, \
    #SAF using ($1<=0.59 && $1>=0.7?$1:2/0):2 w l ls 2 dt 2, \
    #SAF_SD using 1:2 w l ls 2 dt 2, \
    #"<echo '.597 2.12'" with points ls 1, \
    #SAF using ($1<=0.66 && $1>=0.45?$1:2/0):2 w l ls 3 dt 2, \
#plot AF using ($1<=0.32?$1:2/0):2 w l ls 3, \
   # SD using ($1>=0.32?$1:2/0):2 w l ls 3 dt 2, \
    #AF_SD using ($1<=0.43?$1:2/0):2 w l ls 3 dt 2, \
   # SAF_SD using 1:2 w l ls 3 dt 2, \
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

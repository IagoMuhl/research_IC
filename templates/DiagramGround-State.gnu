reset
#set terminal pngcairo dashed enhanced size 600,400 font 'arial,12' fontscale 1.0
#set output 'DiagramaDeEstados.png'

set terminal lua tikz linewidth 5 standalone
set output 'diagramStates.tex' 

#Setup
set encoding utf8
set xlabel "$J_{2} / J_{1}$" font ',14'
set ylabel "$J_{3} / J_{1}$" font ',14'
set yrange [-0.3:0.5]
set xrange [0:1]
set xtics font ',14'
set ytics font ',14'

#Style
#ps-> point size
#pt-> point type
#lc-> line color
set style line 1 lc rgb 'black' pt 7 ps 2  # circle
set style line 5 lc rgb '#E8991E' pt 7 ps 2  # circle
set style line 2 lc rgb '#E8991E'
set style line 3 lc rgb '#070BE8'
set style line 4 lc rgb '#8313E8'

#Label
set label 2 at 0.158, 0.1
set label 2 "$AF$" tc rgb '#8313E8' font ',16'

set label 3 at 0.742, 0.1
set label 3 "$SAF$" tc rgb '#8313E8' font ',16'

set label 5 at 0.455, .45
set label 5 "$SD$" tc rgb '#8313E8' font ',16'

#Legenda
set key off

#Plot
AF = "GroundState/1:2(H)-3:4(V).dat"

plot AF using ($1<=0.5 ?$1:1/0):2 w l ls 3, \
     AF using ($1>=0.5 ?$1:1/0):2 w l ls 3, \
     AF using 3:4 w l ls 3, \
    "<echo '.125 -.15'" with points ls 1, \
    "<echo '.175 -.15'" with points ls 5, \
    "<echo '.225 -.15'" with points ls 1, \
    "<echo '.275 -.15'" with points ls 5, \
    "<echo '.125 -.10'" with points ls 5, \
    "<echo '.175 -.10'" with points ls 1, \
    "<echo '.225 -.10'" with points ls 5, \
    "<echo '.275 -.10'" with points ls 1, \
    "<echo '.125 -0.05'" with points ls 1, \
    "<echo '.175 -.05'" with points ls 5, \
    "<echo '.225 -0.05'" with points ls 1, \
    "<echo '.275 -.05'" with points ls 5, \
    "<echo '.125 -.0'" with points ls 5, \
    "<echo '.175 -0.0'" with points ls 1, \
    "<echo '.225 -.0'" with points ls 5, \
    "<echo '.275 -0.0'" with points ls 1, \
    "<echo '.725 -.15'" with points ls 1, \
    "<echo '.775 -.15'" with points ls 5, \
    "<echo '.825 -.15'" with points ls 1, \
    "<echo '.875 -.15'" with points ls 5, \
    "<echo '.725 -.10'" with points ls 1, \
    "<echo '.775 -.10'" with points ls 5, \
    "<echo '.825 -.10'" with points ls 1, \
    "<echo '.875 -.10'" with points ls 5, \
    "<echo '.725 -0.05'" with points ls 1, \
    "<echo '.775 -.05'" with points ls 5, \
    "<echo '.825 -0.05'" with points ls 1, \
    "<echo '.875 -.05'" with points ls 5, \
    "<echo '.725 -.0'" with points ls 1, \
    "<echo '.775 -0.0'" with points ls 5, \
    "<echo '.825 -.0'" with points ls 1, \
    "<echo '.875 -0.0'" with points ls 5, \
    "<echo '.425 .35'" with points ls 1, \
    "<echo '.475 .35'" with points ls 1, \
    "<echo '.525 .35'" with points ls 5, \
    "<echo '.575 .35'" with points ls 5, \
    "<echo '.425 .3'" with points ls 5, \
    "<echo '.475 .3'" with points ls 5, \
    "<echo '.525 .3'" with points ls 1, \
    "<echo '.575 .3'" with points ls 1, \
    "<echo '.425 .25'" with points ls 1, \
    "<echo '.475 .25'" with points ls 1, \
    "<echo '.525 .25'" with points ls 5, \
    "<echo '.575 .25'" with points ls 5, \
    "<echo '.425 .2'" with points ls 5, \
    "<echo '.475 .2'" with points ls 5, \
    "<echo '.525 .2'" with points ls 1, \
    "<echo '.575 .2'" with points ls 1

#pdftops -eps diagramStates.pdf
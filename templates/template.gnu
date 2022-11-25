reset
set terminal pngcairo dashed enhanced size 600,400 font 'arial,12' fontscale 1.0
set encoding utf8
set output 'figSAF_J3(-0_1)_J2(0_49).png'
set xlabel "T"
set ylabel "m_1"
#m = "./mIsZero.dat"
#n = "./plot.dat"
#o = "./TporM.dat"
set title '4 sítios (SAF): J_1 = 1, J_2 = .49, J_3 = -0.1'

#set key off

plot [-10:10] sin(x)
    #o using 1:2 w l title "J_3 = 0.1"


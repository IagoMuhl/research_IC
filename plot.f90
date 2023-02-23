program plot
    implicit none

    call system("gnuplot -p templates/DiagramPhase.gnu")
    !call system("gnuplot -p templates/DiagramGround-State.gnu")
    !call system("gnuplot -p templates/multiPlot_2x3_DiagramPhase.gnu")
    
end program plot
program data
    implicit none
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: values

    
    ! using keyword arguments
    call date_and_time(date,time,zone,values)
    call date_and_time(DATE=date,ZONE=zone)
    call date_and_time(TIME=time)
    call date_and_time(VALUES=values)
    ! print '(a,2x,a,2x,a)', date, time, zone
    ! print '(8i5)', values
    WRITE (*, '(I2,A,I2,A,I2,A,I2)') values(5),':',values(6),' do dia ',values(3),'/',values(2)

end program data
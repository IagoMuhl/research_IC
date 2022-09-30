program normal
    implicit none
    integer, parameter:: minConfig=2
    integer, parameter:: maxConfig= minConfig**4
    integer:: k
    integer:: s1,s2,s3,s4,m1,m2,m3,m4
    integer, dimension(16,4):: s(maxConfig,4)
    real*8::H,J1,J2,J3
 
    J3=0.05d0; J2=0.3d0; J1=1.d0; k=1
    !!!!!!!!!!    CONTRUÇÃO DA BASE     !!!!!!!!!!!!
    do s1 =-1,1,2
       do s2 =-1,1,2
          do s3 =-1,1,2
             do s4 =-1,1,2
                     
                s(k,1)=s1
                s(k,2)=s2
                s(k,3)=s3
                s(k,4)=s4
                
 
                   
 
                m1=s1
                m2=s2
                m3=s3
                m4=s4
 
                print*, s1,s2,s3,s4
                print *, s(15,4)
                
 
                !!!!!!!!!!!!!!!!! HAMILTONIANO !!!!!!!!!!!!!!
 
                H = -J1*(((s1*m2)-(m1*m2)/2.d0) &
                &+(s1*m3)-(m1*m3)/2.d0) &
                &-(3*J2)*((s1*m4)-(m1*m4)/2.d0)&
                &-(4*J3)*((s1*m1)-((m1)**2.d0)/2.d0)
 
 
                H = H +J1*(((s2*m1)-(m2*m1)/2.d0)&
                &+(s2*m4)-(m2*m4)/2.d0)&
                &+(3*J2)*((s2*m3)-(m2*m3/2.d0))&
                &+(4*J3)*((s2*m2)-((m2)**2.d0)/2.d0)
 
 
                H = H -J1*(((s3*m1)-(m3*m1)/2.d0) &
                &+(s3*m4)-(m3*m4)/2.d0) &
                &-(3*J2)*((s3*m2)-(m3*m2)/2.d0)&
                &-(4*J3)*((s3*m3)-((m3)**2.d0)/2.d0)
 
 
                H = H +J1*(((s4*m2)-(m4*m2)/2.d0)&
                &+(s4*m3)-(m4*m3)/2.d0)&
                &+(3*J2)*((s4*m1)-(m4*m1/2.d0))&
                &+(4*J3)*((s4*m4)-((m4)**2.d0)/2.d0)
 
                
 
                k=k+1
 
             end do
          end do
       end do
    end do
 
 
 
 
 end program normal
 
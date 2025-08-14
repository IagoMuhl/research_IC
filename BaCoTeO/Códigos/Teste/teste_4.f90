program teste_hexa
    implicit none
    real*8:: J1, J2, J3, T, H, H_final, step, tol, a, b, erro1(6), erro2(6), Z1, Z2
    real*8:: m_1(6), m_2(6), m1_erro, m2_erro, F1, F2, F
    real*8:: H_intra(64), m1(6), m2(6), H_inter1(64), H_inter2(64), s_z(64), H1(64), H2(64)
    integer:: k, s(64,6), s1, s2, s3, s4, s5, s6, i, p, j

    J1 = 1.d0; J2 = 0.8d0; J3 = -0.111d0
    T = 0.05d0
    H = 0.9d0
    H_final = 0.d0

    s_z = 0.d0
    step = 10.d0**(-3)
    tol = 10.d0**(-8)
    a = 0.d0
    b = 0.d0
    Z1 = 0.d0; Z2 = 0.d0

    k = 1

    do s1 = -1, 1, 2
        do s2 = -1, 1, 2
            do s3 = -1, 1, 2
                do s4 = -1, 1, 2
                    do s5 = -1, 1, 2
                        do s6 = -1, 1, 2

                            s(k,1) = s1
                            s(k,2) = s2
                            s(k,3) = s3
                            s(k,4) = s4
                            s(k,5) = s5
                            s(k,6) = s6

                            k = k + 1

                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    do i = 1, 64
        H_intra(i) = J1*(s(i,1)*s(i,2) + s(i,2)*s(i,3) + s(i,3)*s(i,4) &
        & + s(i,4)*s(i,5) + s(i,5)*s(i,6) + s(i,6)*s(i,1)) &
        & + J2*(s(i,1)*s(i,3) + s(i,3)*s(i,5) + s(i,5)*s(i,1) &
        & + s(i,2)*s(i,4) + s(i,4)*s(i,6) + s(i,6)*s(i,2)) &
        & + J3*(s(i,1)*s(i,4) + s(i,2)*s(i,5) + s(i,3)*s(i,6))
    end do

    do p = 1, 6
        do i = 1, 64
            s_z(i) = s_z(i) + s(i,p)
        enddo
    enddo

    if (H<H_final) then
        step  = step
    else
        step = -step
    endif

    m1 = [-1.d0, 1.d0, 1.d0, -1.d0, 1.d0, 1.d0]
    m2 = [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0]

    open(10)
    open(11)
    open(12)

    do while (H/=H_final)

        erro1 = 1.d0; erro2 = 1.d0

        do while(abs((m1_erro+m2_erro)/2)>= tol)

            H_inter1 = 0.d0; H_inter2 = 0.d0

            do i = 1, 64

                H_inter1(i) = J1*(m2(1)*(s(i,1)+s(i,4)-m1(1)) + m1(2)*(s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m1(2))) &
                & +  J2*(m2(1)*(2*s(i,1)+2*s(i,4)-2*m1(1)+s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m1(2)) &
                & + m1(1)*(s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m1(2)) &
                & + 2*m1(2)*(s(i,1)+s(i,4)-m1(1)+s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m1(2))) &
                & + J3*(2*m1(1)*(s(i,1)+s(i,4)-m1(1)) + (m2(1) + m1(2))*(s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m1(2)))

                ! H_inter1(i) = J1*((s(i,1)-m1(1)/2)*(m2(4)) + (s(i,2)-m1(2)/2)*(m1(6)) + (s(i,3)-m1(3)/2)*(m1(5)) &
                ! & + (s(i,4)-m1(4)/2)*(m2(1)) + (s(i,5)-m1(5)/2)*(m1(3)) + (s(i,6)-m1(6)/2)*(m1(2))) &
                ! & + J2*((s(i,1)-m1(1)/2)*(m1(2)+m2(5)+m2(3)+m1(6)) + (s(i,2)-m1(2)/2)*(m2(4)+m1(1)+2*m1(5)) &
                ! & + (s(i,3)-m1(3)/2)*(m1(4)+m2(1)+2*m1(6)) + (s(i,4)-m1(4)/2)*(m1(5)+m2(2)+m2(6)+m1(3)) &
                ! & + (s(i,5)-m1(5)/2)*(m1(4)+m2(1)+2*m1(2)) + (s(i,6)-m1(6)/2)*(m1(1)+m2(4)+2*m1(3))) &
                ! & + J3*((s(i,1)-m1(1)/2)*(2*m1(1)) + (s(i,2)-m1(2)/2)*(m1(6)+m2(3)) + (s(i,3)-m1(3)/2)*(m1(5)+m2(2)) &
                ! & + (s(i,4)-m1(4)/2)*(2*m1(4)) + (s(i,5)-m1(5)/2)*(m1(3)+m2(6)) + (s(i,6)-m1(6)/2)*(m1(2)+m2(5)))

                H_inter2(i) = (s(i,1)+s(i,2)+s(i,3)+s(i,4)+s(i,5)+s(i,6)-3*m2(1))*(J1*m1(1) + J2*2*(m1(1)+m1(2)) + J3*2*m1(2))

            enddo

            H1 = H_intra + H_inter1 - H*s_z
            H2 = H_intra + H_inter2 - H*s_z


               if (T<=10.d0**(-1)) then

                  a = minval(H1)

                  H1 = H1 - a

                  b = minval(H2)

                  H2 = H2 - b

               endif

            do i = 1, 64
                Z1 = Z1 + dexp(-H1(i)/T)
                Z2 = Z2 + dexp(-H2(i)/T)
            enddo

            m_1 = m1
            m_2 = m2

            m1 = 0.d0; m2 = 0.d0

            do i = 1, 6

                do j = 1, 64

                    m1(i) = m1(i) + s(j,i)*dexp(-H1(j)/T)
                    m2(i) = m2(i) + s(j,i)*dexp(-H2(j)/T)

                enddo

                m1(i) = m1(i)/Z1
                m2(i) = m2(i)/Z2

            enddo
            

            erro1 = abs(m1 - m_1)
            erro2 = abs(m2 - m_2)

            m1_erro = maxval(erro1)
            m2_erro = maxval(erro2)

        enddo

        F1 = -T*dlog(Z1)
        F2 = -T*dlog(Z2)

        F1 = F1 + a
        F2 = F2 + b

        F = (3*F1 + F2)/4

        write(10,*) H, F, F1, F2
        write(11,*) H, m1
        write(12,*) H, m2

        print*, H

        H = H + step

        if ((max(H,H_final))-(min(H,H_final))<=abs(step)) then
               H = H_final
        endif

    enddo

    close(10)
    close(11)
    close(12)

    
end program teste_hexa
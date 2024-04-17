!        This is file : main
! Author = ufokis
! Started at: 07.04.2024
! Last Modified: qui 11 abr 2024 23:46:19
!
Program main

    use GPDsigmaPhotonProton
    implicit none
    real(wp):: Amplitude, Y, W, zeta, SNN, dsigma, bv, dsigma2
    integer:: IY, i, ier1, ier2
    SNN = 5020.0_wp
    IY = -60
    Y = 0.1_wp*float(IY)
    W = M_V*SNN*EXP(Y)
    open (unit = 2, file='dados3.dat', status='unknown')
    do i = 0, 2000
        W = float(i) + 1.d-5
        zeta = (M_V**2)/(W**2)
        XI = zeta/(2.0_wp-zeta)
        bv = 4.9_wp+4*(0.06)*log(W/90.0_wp)

        dsigma2 = (1/bv) *  calculo_derivative_cross_section()

        call AmpPhotonProton(Amplitude, ier1, ier2)
        !if (ier1 == 1 .AND. ier2 == 1) then
        dsigma = (1/(W**4))*(4*(pi**2)* &
                             (1.0_wp/137.0_wp)* &
                             ((e_q**2))/(9*(xi**2)))* &
                 (1.32_wp/(m_c**3))*Amplitude/bv
        print *, "Amplitude", W, dsigma*0.3894D3, dsigma2*0.3894D3
        write (2, *) W, dsigma*0.3894D3, dsigma2*0.3894D3
        !end if
    end do
end program main

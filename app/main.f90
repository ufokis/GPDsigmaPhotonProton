!	This is file : main
! Author= ufokis
! Started at: 07.04.2024
! Last Modified: Sun Apr  7 12:32:43 2024
!
Program  main
   use GPDsigmaPhotonProton
   implicit none
   real(wp)::Amplitude, Y, W, zeta, SNN, dsigma, bv
   integer::IY, i
   SNN = 5020.0_wp
   IY = -60
   Y = 0.1_wp*float(IY) 
   W = M_V*SNN*EXP(Y)
   open(unit = 2, file = 'dados2.dat', status = 'unknown')
   do i = 0, 2000
        W =  float(i)
        zeta = (M_V**2)/(W**2)
        XI = zeta/(2.0_wp - zeta)
        bv =  4.9_wp +  4*(0.06)*log(W/90.0_wp)
        call AmpPhotonProton(Amplitude)
        dsigma =  (1/(W**4)) *( 4*(pi**2) * &
            (1.0_wp/137.0_wp) *((e_q**2))/(9*(xi**2)))*&
            (1.32_wp/(m_c**3)) * Amplitude/bv 
       print *, "Amplitude", W, dsigma*0.3894_wp 
       write(2, * ) W, dsigma*0.3894_wp 
   end do
end program main

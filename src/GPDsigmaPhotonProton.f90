module GPDsigmaPhotonProton
    !! Módulo $ \sigma(\gamma+p --> J/\Psi P) $
   use Teste_GPDs
   use FunctionsMath
   use quadpack, wp => quadpack_RK
   ! Constantes
   implicit none
   public
   real(8), parameter:: pi = 3.1415926d0
   real(8), parameter:: m_c = 1.55d0
   real(8):: alpha_s
   real(8), parameter:: n_f = 4.d0
   real(8), parameter:: M_V = 2.d0*m_c
   real(8), parameter:: mu_f = m_c
   real(8), parameter:: mu_r = mu_f
   real(8), parameter:: N_c = 3.d0
   real(8), parameter:: QA = mu_f
   real(8), parameter:: C_f = 4/3.d0
   real(8), parameter:: e_q = 2/3.d0
   real(8):: xi
   real(8), parameter:: c_1 = C_f
   real(8), parameter:: c_2 = 1/(2.d0*N_c)
   real(8), parameter:: beta_0 = (11*N_c/3.D0) - (2*n_f/3.D0)
   character(len = 64), parameter:: name = 'nCTEQ15WZSIH_208_82'

contains
   Subroutine AmpPhotonProton(Amplitude, ier1, ier2)
      !integer, intent(out):: ier1, ier2
      !real(wp), parameter:: epsabs = 0.0_wp
      !real(wp), parameter:: epsrel = 1.5d-3
      !10**(log10(epsilon(1.0_wp))/2.0_wp+1)
      !real(wp):: err
      real(wp):: xmin, xmax, Amplitude
      real(wp):: AmpIntImg, AmpIntReal
      !print *, "tmin", tmin, epsrel      
      real(wp), parameter:: epsabs = 0.0_wp
      real(wp), parameter:: epsrel = 10**(log10(epsilon(1.0_wp))/2.0_wp+1)
      real(wp), parameter:: a = 1.d-7
      real(wp), parameter:: b = 0.99_wp
      integer, parameter:: key = 4
      integer, parameter:: limit = 100
      integer, parameter:: lenw = limit*4
      real(wp):: abserr, result, work(lenw)
      integer:: ier1, ier2, iwork(limit), last, neval
      
      call InitPDFsetByName(name)
      call InitPDF(0)
      call GetXmin(0, xmin)
      call GetXmax(0, xmax)

      call dqag(IntImg, a, b, epsabs, epsrel, key, AmpIntImg, &
            abserr, neval, ier1, limit, lenw, last, &
            iwork, work)

      call dqag(IntReal, a, b, epsabs, epsrel, key, AmpIntReal, &
            abserr, neval, ier2, limit, lenw, last, &
            iwork, work)
      !call dgauss8(IntImg, 1.d-7, xmax, &
      !epsrel, AmpIntImg, ier1, err)
      !call dgauss8(IntReal, 1.d-7, xmax, &
      !epsrel, AmpIntReal, ier2, err)
      !if (ier1 == 1 .AND. ier2 == 1) then
         Amplitude = (AmpIntReal**2 + AmpIntImg**2)
      !end if
      !AmpIntImg  = SGS3(tmin, 0.99d0, 1.d-5, IntImg)
      !AmpIntReal = SGS3(tmin, 0.99d0, 1.d-5, IntReal)
   End Subroutine AmpPhotonProton

   function IntImg(x) result(ResultIntImg)
      real(wp), intent(in):: x
      real(wp)  :: ResultIntImg
      ResultIntImg = aimag(IntComplex(x))
   end function IntImg

   function IntReal(x) result(ResultIntReal)
      real(wp), intent(in):: x
      real(wp):: ResultIntReal
      ResultIntReal = real(IntComplex(x))
   end function IntReal

   function IntComplex(x) result(ResultIntComplex)
      real(wp), intent(in)  :: x
      complex(wp):: ResultIntComplex
      real(wp)   :: PDF(-6:6), I1, I2
      call InitPDFsetByName(name)
      call InitPDF(0)

      alpha_s = alphasPDF(QA)
      ! write(*,*) ' aqui ', alpha_s
      call evolvePDF(x, QA, PDF)

      I1 = 2*PDF(0)
      I2 = (sum(PDF(1:4)) + sum(PDF(-1:-4)))/x

      ResultIntComplex = I1*Tg(x) + I2*Tq(x)

   end function IntComplex

   function Tq(x) result(ResultTq)
      real(wp)              :: x
      complex(wp)           :: ResultTq
      complex(wp):: epsilon
      epsilon = xi*dcmplx(0.d0, 1.d-8)
      !print*, "epsilon", epsilon, 1.d-8, xi
      ResultTq = (C_f*alpha_s**2/(2*pi))*fq((x-xi+epsilon)/(2*xi))
   end function Tq

   function fq(y) result(Resultfq)
      complex(wp):: y, fq_1, fq_2, fq_3, fq_4, &
                    fq_5, fq_6, fq_7, fq_8, fq_9, fq_10, Resultfq

      fq_1 = (log(4*(m_c**2)/(mu_f**2)))*(1+2*y)
      fq_2 = (log(-y)/(1+y) - log(1+y)/y)
      fq_3 = (PI**2)*(13*(1+2*y))/(48*y*(1+y))
      fq_4 = 2*log(2.d0)/(1+2*y)
      fq_5 = (log(-y) + log(1+y))/(1+2*y)
      fq_6 = (1+2*y)*(((log(-y))**2)/(1+y) &
                        - (log(1+y))**2/(y))
      fq_7 = ((3-4*y+16*y*(1+y))/(4*y*(1+y)))
      fq_8 = cdli2(1+2*y)
      fq_9 = ((7+4*y+16*y*(1+y))/(4*y*(1+y)))
      fq_10 = cdli2(-1-2*y)

      Resultfq = fq_1*fq_2-fq_3+fq_4+fq_5+fq_6 + &
                 fq_7*fq_8-fq_9*fq_10
   end function fq

   function Tg(x) result(ResultTg)
      real(wp)              :: x
      complex(wp)           :: ResultTg, ResultTgLO, ResultTgNLO

      complex(wp):: epsilon
      !print *, ' alhpas = ', alpha_s
      epsilon = xi*dcmplx(0.d0, 1.d-8)
      ResultTgLO = (xi/((x-xi+epsilon)* &
                        (x+xi-epsilon)))*alpha_s
      epsilon = xi*dcmplx(0.d0, 1.d-5)
      ResultTgNLO = (xi/((x-xi+epsilon)* &
                         (x+xi-epsilon)))*((alpha_s**2)/(4*pi)) &
                    *fg((x-xi+epsilon)/(2*xi))
      ResultTg = ResultTgLO+ResultTgNLO
   end function Tg
   function fg(y) result(Resultfg)

      complex(wp)   :: Resultfg
      complex(wp):: y, fg_1, fg_2, fg_3, fg_4, fg_5, &
                    fg_6, fg_7, fg_8, fg_9, fg_10
      complex(wp):: fg_11, fg_12, fg_13, fg_14, &
                    fg_15, fg_16, fg_17, fg_18, fg_19
      complex(wp):: fg_20, fg_21, fg_22, fg_23, &
                    fg_24, fg_25, fg_26, fg_27, fgt_1
      complex(wp):: fgt_2, fgt_3, fgt_4, fgt_5

      fg_1 = 4*(c_1-c_2)*(1+2*y*(1+y))
      fg_2 = (log(-y)/(1+y) - log(1+y)/y)
      fg_3 = (log(4*(m_c**2)/(mu_f**2)))
      fg_4 = beta_0*log((mu_r**2)/(mu_f**2))
      fg_5 = fg_1
      fg_6 = (((log(-y)**2)/(1+y)) - ((log(1+y)**2)/y))
      fg_7 = 8*c_1
      fg_8 = (2+y*(1+y)*(25+88*y*(1+y))) &
             *c_1/(48*(y**2)*((1+y)**2))
      fg_9 = (10+y*(1+y)*(7-52*y*(1+y))) &
             *c_2/(24*(y**2)*((1+y)**2))
      fg_10 = c_1*(1+6*y*(1+y)*(1+2*y*(1+y))) &
              /(y*(1+y)*((1+2*y)**2))
      fg_11 = c_2*((1+2*y)**2)/(y*(1+y))
      fg_12 = PI*(SQRT(-y*(1+y))/(y*(1+y)))* &
              ((7/2.d0)*c_1-3*c_2)
      fg_13 = 2*c_2*SQRT(-y*(1+y))/(y*(1+y))
      fg_14 = ((1+4*y)/(1+y))*ATAN(SQRT(-y/(1+y)))
      fg_15 = ((3+4*y)/y)*ATAN(SQRT((1+y)/(-y)))
      fg_16 = ((ATAN(SQRT(-y/(1+y)))**2)/(2*y*(1+y)))
      fg_17 = (7+4*y)*c_1
      fg_18 = 2*(1+2*y-2*(y**2))*c_2/(1+y)
      fg_19 = (ATAN(SQRT((1+y)/(-y)))**2)/(2*y*(1+y))
      fg_20 = (3-4*y)*c_1
      fg_21 = 2*(3+6*y+2*(y**2))*c_2/y
      fg_22 = 2.d0*a_1(y)*log(-y)
      fg_23 = 2*log(1+y)*a_1(-y-1)
      fg_24 = 2.d0*a_2(y)
      fg_25 = cdli2(1+2*y)
      fg_26 = 2.d0*a_2(-1-y)
      fg_27 = cdli2(-1-2*y)

      fgt_1 = fg_1*fg_2*fg_3+fg_4+fg_5*fg_6-fg_7
      fgt_2 = -(PI**2)*(fg_8+fg_9) - log(2.d0)*(fg_10+fg_11)
      fgt_3 = fg_12+fg_13*(fg_14+fg_15) - fg_16*(fg_17-fg_18)
      fgt_4 = -fg_19*(fg_20-fg_21) + fg_22+fg_23
      fgt_5 = fg_24*fg_25+fg_26*fg_27

      Resultfg = fgt_1+fgt_2+fgt_3+fgt_4+fgt_5
   end function fg

   function a_1(y) result(Resulta_1)
      complex(wp):: y, Resulta_1

      Resulta_1 = 5.0D0 + (16.0D0*y) - &
                  (6.0D0/(1.D0+y))
      Resulta_1 = Resulta_1 + (1.0D0/(1.0D0+2.0D0*y)**2)
      Resulta_1 = (c_1/4.0D0)*(Resulta_1 - &
                               (5.0D0/(1.0D0+2.0D0*y)))
      Resulta_1 = Resulta_1 - (c_2/2.0D0)* &
                  (2.0D0 + (3.0D0/y))
      Resulta_1 = Resulta_1 - (c_2/2.0D0)* &
                  ((8.0D0*y) - (1.0D0/(1.0D0+y)))
   end function a_1

   function a_2(y) result(Resulta_2)
      complex(wp):: y, Resulta_2

      Resulta_2 = 12.0D0 + (9.0D0/y) + (64.0D0*y) &
                  - (2.0D0/(1.0D0+y)**2)
      Resulta_2 = (c_1/8.0D0)*(Resulta_2 + &
                               (21.0D0/(1.0D0+y)) - (4.0D0/(1.0D0+2.0D0*y)))
      Resulta_2 = Resulta_2 - (c_2/4.0D0)* &
                  (8.0D0 + (3.0D0/y**2) + (11.0D0/y))
      Resulta_2 = Resulta_2 - (c_2/4.0D0)* &
                  ((32.0D0*y) - (2.0D0/(1.0D0+y)**2) &
                   + (9.0D0/(1.0D0+y)))
   end function a_2
end module GPDsigmaPhotonProton


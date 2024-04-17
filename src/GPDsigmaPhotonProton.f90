module GPDsigmaPhotonProton
    !! Módulo que calcula \(\sigma^{\gamma+p \rightarrow J/\Psi+p} \) e \( \frac{id\sigma^{\gamma+p \rightarrow j/\Psi+p}}{dt}\Big|_{t = 0}\)
    use FunctionsMath
    use Teste_GPDs
    use quadpack
    use iso_fortran_env, only: wp => real64  ! double precision
    implicit none
    ! Argumentos do  modulos
    real(kind = wp), parameter:: pi = 3.1415926_wp         !! Número de \(\pi\)
    real(kind = wp), parameter:: m_c = 1.55_wp             !! Massa do quark charm
    real(kind = wp), parameter:: n_f = 4.0_wp              !! Número de sabores ativos de sabores
    real(kind = wp), parameter:: M_V = 2.0_wp*m_c          !! Massa do meson Vetorial
    real(kind = wp), parameter:: mu_f = M_V                !! escala de fatorização
    real(kind = wp), parameter:: mu_r = mu_f               !! Fator de renormalização
    real(kind = wp), parameter:: N_c = 3.0_wp              !! Numero de cor
    real(kind = wp), parameter:: QA = mu_f                 !! momento transferido
    real(kind = wp), parameter:: C_f = 4/3.0_wp            !! Fator de cor
    real(kind = wp), parameter:: e_q = 2/3.0_wp            !! Carga do quark
    real(kind = wp), parameter:: c_1 = C_f                 !! Estrutura de cor \(c_1\)
    real(kind = wp), parameter:: c_2 = -1/(2.0_wp*N_c)     !! Estrutura de cor \(c_2\)
    real(kind = wp), parameter:: beta_0 = (11.0_wp* &
                                         N_c/3.0_wp) - (2*n_f/3.0_wp) !! Parametro \(b_0\)
    character(len = 64), parameter:: name = &
                                   'nCTEQ15WZSIH_208_82'   !! Nome da PDF utilizada

    real(kind = wp)           :: alpha_s                   !! Constante de aclopamento da QCD
    real(kind = wp)           :: xi                        !! Parametro de assimetria

contains

    function calculo_derivative_cross_section() &
        result(derivative_cross_section)
        !! Função que calcula \( \frac{id\sigma^{\gamma+p \rightarrow j/\Psi+p}}{dt}\Big|_{t = 0} \)
        real(kind = wp), parameter:: Gamma_Vll = 5.55*1.d-6      !! Tamanho do decaimento  do vetor Meson em \(l^+l^-\)
        real(kind = wp), parameter:: alpha_QED = 1.0_wp/137.0_wp !! Constante de aclopamento da QED
        real(kind = wp)           :: derivative_cross_section    !! Retorno de \(\frac{d\sigma^{\gamma+p \rightarrow j/\Psi+p}}{dt}\Big|_{t = 0} \)
        real(kind = wp)           :: PDF(-6:6), F_g              !! Vetor que recebe as  PDF \(s, c, b, g, \bar{b}, \bar{c}, \bar{s} \)

        call InitPDFsetByName(name)  ! Inicializa a pdf escolhida
        call InitPDF(0)              ! Inicializa a versão
        call evolvePDF(xi, QA, PDF)  ! adiciona valor as PDF's
        
        alpha_s = alphasPDF(QA)      ! Inicializa a constante de acoplamento da QCD
        F_g     = PDF(0)             ! Calcula a GPD \(F_g\)
        
        derivative_cross_section = (16.0_wp*pi**3/3.d0) &
                                   *(Gamma_Vll/M_V**5)* &
                                   ((alpha_s**2)*(xi*F_g)**2) &
                                   /(alpha_QED)
    end function calculo_derivative_cross_section

    Subroutine AmpPhotonProton(Amplitude, ier1, ier2)
        !! Subrotina que retorna a Amplitude (\(A^g_\text{LO}, A^g_\text{NLO}, A^q_\text{NLO}\))
        !integer, intent(out):: ier1, ier2
        !real(wp), parameter:: epsabs = 0.0_wp
        !real(wp), parameter:: epsrel = 1.5d-3
        !10**(log10(epsilon(1.0_wp))/2.0_wp+1)
        !real(wp):: err
        real(kind = wp):: xmin, xmax, Amplitude
        real(kind = wp):: AmpIntImg, AmpIntReal
        !print *, "tmin", tmin, epsrel
        real(kind = wp), parameter:: epsabs = 0
        real(kind = wp), parameter:: epsrel = 10**(log10(epsilon(1.0_wp))/2.0_wp+1)
        real(kind = wp), parameter:: a = 1.d-7
        real(kind = wp), parameter:: b = 0.99_wp
        integer, parameter:: key = 4
        integer, parameter:: limit = 100
        integer, parameter:: lenw = limit*4
        real(kind = wp):: abserr, result, work(lenw)
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
        real(kind = wp), intent(in):: x
        real(kind = wp)  :: ResultIntImg
        ResultIntImg = aimag(IntComplex(x))
    end function IntImg

    function IntReal(x) result(ResultIntReal)
        real(kind = wp), intent(in):: x
        real(kind = wp)  :: ResultIntReal
        ResultIntReal = real(IntComplex(x))
    end function IntReal

    function IntComplex(x) result(ResultIntComplex)
        real(kind = wp), intent(in):: x
        complex(kind = wp)  :: ResultIntComplex
        real(kind = wp)     :: PDF(-6:6), I1, I2
        
        call InitPDFsetByName(name)
        call InitPDF(0)
        call evolvePDF(x, QA, PDF)
        
        alpha_s = alphasPDF(QA)
        I1      = 2*PDF(0)
        I2      = (sum(PDF(1:4)) + sum(PDF(-1:-4)))/x

        ResultIntComplex = I1*Tg(x)! + I2*Tq(x)

    end function IntComplex

    function Tq(x) result(ResultTq)
        real(kind = wp)              :: x
        complex(kind = wp)           :: ResultTq
        complex(kind = wp):: epsilon
        epsilon = xi*dcmplx(0.d0, 1.d-8)
        !print*, "epsilon", epsilon, 1.d-8, xi
        ResultTq = (C_f*alpha_s**2/(2*pi))*fq((x-xi+epsilon)/(2*xi))
    end function Tq

    function fq(y) result(Resultfq)
        complex(kind = wp):: y, fq_1, fq_2, fq_3, fq_4, &
                      fq_5, fq_6, fq_7, fq_8, fq_9, &
                      fq_10, Resultfq

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

    function Tg_LO(x) result(ResultTg)
        real(kind = wp)              :: x
        complex(kind = wp)           :: ResultTg, ResultTgLO, ResultTgNLO

        complex(kind = wp):: epsilon
        !print *, ' alhpas = ', alpha_s
        epsilon = xi*dcmplx(0.d0, 1.d-8)
        ResultTgLO = (xi/((x-xi+epsilon)* &
                          (x+xi-epsilon)))*alpha_s
        ResultTg = ResultTgLO!+ResultTgNLO
    end function Tg_LO
    
    Function  Tg_NLO(x) result (ResultTg_NLO)
        real(kind = wp)     :: x
        complex (kind = wp):: ResultTg_NLO 
        
        epsilon      = xi*dcmplx(0.d0, 1.d-5)
        ResultTg_NLO = (xi/((x-xi+epsilon)*(x+xi-epsilon)))*&
                       ((alpha_s**2)/(4*pi))*&
                       fg((x-xi+epsilon)/(2*xi))
    End Function Tg_NLO 
    
    function fg(y) result(Resultfg)

        complex(wp)   ::  y, fg_1, fg_2, fg_3, fg_4, fg_5, &
                          fg_6, fg_7, fg_8, fg_9, fg_10, &
                          fg_11, fg_12, fg_13, fg_14, &
                          fg_15, fg_16, fg_17, fg_18, fg_19, &
                          fg_20, fg_21, fg_22, fg_23, &
                          fg_24, fg_25, fg_26, fg_27
        complex(wp)   ::  fgt_2, fgt_3, fgt_4, fgt_5
        complex(wp)   :: Resultfg

        fg_1 = 4 * (c_1-c_2) * (1+2 * y * (1+y))
        fg_2 = (log(- y)/(1+y) - log(1+y)/y)
        fg_3 = (log(4 * (m_c**2)/(mu_f**2)))
        fg_4 = beta_0*log((mu_r**2)/(mu_f**2))
        fg_5 = fg_1
        fg_6 = (((log( - y)**2)/(1+y)) - ((log(1+y)**2)/y))
        fg_7 = 8*c_1
        fg_8 = (2+y * (1+y) * (25+88*y * (1+y))) &
               * c_1/(48 * (y**2) * ((1+y)**2))
        fg_9 = (10+y * (1+y) * (7-52*y * (1+y))) &
               * c_2/(24 * (y**2)*((1+y)**2))
        fg_10 = c_1 * (1+6 * y * (1+y) * (1+2 * y * (1+y))) &
                /(y * (1+y) * ((1+2 * y)**2))
        fg_11 = c_2 * ((1+2 * y)**2)/(y * (1+y))
        fg_12 = PI * (SQRT(-y * (1+y))/(y * (1+y))) * &
                ((7/2.d0) * c_1-3 * c_2)
        fg_13 = 2*c_2*SQRT(- y * (1+y))/(y * (1+y))
        fg_14 = ((1+4 * y)/(1+y)) * ATAN(SQRT(- y/(1+y)))
        fg_15 = ((3+4 * y)/y)*ATAN(SQRT((1+y)/(- y)))
        fg_16 = ((ATAN(SQRT(- y/(1+y)))**2)/(2*y *(1+y)))
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
        fgt_2 = - (PI**2) * (fg_8+fg_9) - log(2.d0) * (fg_10+fg_11)
        fgt_3 = fg_12+fg_13 * (fg_14+fg_15) - fg_16 * (fg_17-fg_18)
        fgt_4 = - fg_19 * (fg_20-fg_21) + fg_22+fg_23
        fgt_5 = fg_24*fg_25+fg_26*fg_27

        Resultfg = fgt_1+fgt_2+fgt_3+fgt_4+fgt_5
    end function fg

    function a_1(y) result(Resulta_1)
        !! Função auxiliar \(a_2\)
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
        !! Função auxiliar \(a_2\)
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


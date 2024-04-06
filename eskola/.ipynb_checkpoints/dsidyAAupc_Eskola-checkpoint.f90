!
! This is file : dsidyAAupc_Eskola
! Author = eliton
! Started at: 23.01.2024
! Last Modified: seg 18 mar 2024 23:37:13
!
Program  dsidyAAupc_Eskola
    use sigmaPhotonProton
    Implicit None
		

 
    real(8):: sig
    real(8):: rapidez, W, zeta, omega, GAM, TEMPO_MIN
    real(8):: sigFA, DSIGDYAA, ESPECTRO, formFactor
    integer:: IY, ICONT
    real(8):: SNN, IDP
    integer:: IDFLUXO
    real(8), external   :: formFactorF
    real(8):: constDecay

    COMMON/IDP/ IDP
	IDP = 1
	 
	IDFLUXO = 1
    SNN = 5020.D0

    open(unit = 8, file='dados2.dat', status='unknown')

    do IY = -60, 60
        rapidez = 0.1D0*float(IY)
        DSIGDYAA = 0.D0
        ICONT = 1 
        do while (ICONT < 3)
            OMEGA = 0.5D0*M_V*EXP(rapidez)
            W = M_V*SNN*EXP(rapidez)
            zeta = (M_V**2)/W
            XI = zeta / (2.D0-zeta)

            ! Fator de Forma
            GAM = SNN / (2.D0*0.938D0)
            TEMPO_MIN = ((M_V**2.D0)/(2.D0*OMEGA*GAM))**2.D0
            !write(6, *) 'esse', TEMPO_MIN
            
            call AmpPhotonProton(sig)
            constDecay = (5.55d-3)*((3.d0/2.d0)**3)*(137.d0/(1.3*pi))/((1-8*alpha_s/(3.d0*PI))**2)
            write(*,*)'esse ',constDecay
            call fluxoFotonsNucleo(SNN, OMEGA, ESPECTRO)
            
            
            formFactor = SGS1(TEMPO_MIN, 1.D0, 1.D-3, formFactorF)
            
            ! Cálculo da seção de choque
            sigFA = (1 / (W**2)) * ((4 * (PI**2) * (1/137.0D0) * (e_q**2)) / (9 * (XI**2))) *constDecay*sig*formFactor
            DSIGDYAA = DSIGDYAA+OMEGA*ESPECTRO*sigFA

            ICONT = ICONT+1
            rapidez = -rapidez
            end do
            
            DSIGDYAA = DSIGDYAA*0.3894D0

        ! Escreve no arquivo
        write(*, *) rapidez, DSIGDYAA     
        ! Escreve no console
        write(8, *) rapidez, DSIGDYAA
        end do

    close(unit = 8)
    
end program dsidyAAupc_Eskola
!------------------------------------------------------------------------------

!==============================================================================
! FATOR DE FORMA
!------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION formfactorF(t)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN):: t
    DOUBLE PRECISION:: PI, q, formfactor
    common/q/ q
    ! Parâmetros e constantes
    PARAMETER (PI = 3.141592d0)

    ! Cálculo de q
    q = SQRT(t)

    ! Chamada para a função WoodsSaxon_Yukawa
    CALL WoodsSaxon_Yukawa(formfactor)  ! Certifique-se de ajustar tmin em WoodsSaxon_Yukawa

    ! Cálculo do quadrado do form factor
    formfactorF = formfactor**2.d0

    RETURN
END FUNCTION formfactorF

!==============================================================================
! FLUXO DE FOTONS
!------------------------------------------------------------------------------    
DOUBLE PRECISION FUNCTION fluxoFotonsF(BB)
    
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN):: BB
    DOUBLE PRECISION:: PI, flux
    DOUBLE PRECISION:: B1, B
    
    INTEGER:: IDFLUXO

    ! Variáveis comuns
    COMMON/IDFLUXO/IDFLUXO
    COMMON/B1/B1
    COMMON/B/ B


    ! Parâmetros e constantes
    PI = 3.141592d0

    ! Seleção do fluxo com base em IDFLUXO
    IF (IDFLUXO == 1) THEN
        B1 = BB
        CALL fluxoN0(flux)
        fluxoFotonsF = 2.D0*PI*B1*flux
    ELSE
        B = BB
        CALL fluxoN2(flux)
        fluxoFotonsF = 2.D0*PI*B * flux
    ENDIF

    RETURN
END FUNCTION fluxoFotonsF


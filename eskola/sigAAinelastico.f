
c	sigAAinelastico.f
c	TAB.f
c	rho.f
c	sgs.f
c	dadmul.f

C	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
C	PARAMETER (PI=3.141592)
	
c	************** ARQUIVOS DE SAÍDA ******************
C	OPEN (UNIT=10,FILE='sigAAinel.dat',STATUS='unknown')	!WOODS-SAXON (WS)
c	***************************************************

C	BMIN = 6.04D0*5.07D0
C	BMAX = 6.98D0*5.07D0
	
C	CALL sigAAinel(bmin,bmax,sigAAinelastico)
	
C	WRITE(10,*) sigAAinelastico
C	WRITE(*,*) sigAAinelastico

C	END


	SUBROUTINE sigAAinel(bmin,bmax,fc_sigAAin)	

        IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	PARAMETER (PI=3.141592)	

	EXTERNAL sigAAinelF

	EPS = 1.D-3

	fc_sigAAin = SGS1(bmin,bmax,EPS,sigAAinelF)

	RETURN
	END

C	***************************************************
C			 WOODS - SAXON 
C	***************************************************

	DOUBLE PRECISION FUNCTION sigAAinelF(BB)

        IMPLICIT DOUBLE PRECISION (A-H,K-Z)
	PARAMETER (PI=3.141592)

	COMMON/B/ B

	B = BB
	
	call TAB(overlapN)

	SIGNN = 67.6D0/0.3894D0			!valor para 5.02 TeV --> convertido mb --> GeV⁻²

	sigAAinelF = 2.D0*PI*B*(1.D0-EXP(-SIGNN*overlapN))

	RETURN
	END





c	TA.f
c	rho.f
c	sgs.f

c	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
c	PARAMETER (PI=3.141592)
c	
c	COMMON/bglau/ bglau
c
cc	************** ARQUIVOS DE SAÍDA ******************
c	OPEN (UNIT=10,FILE='TA208.dat',STATUS='unknown')	!WOODS-SAXON (WS)
cc	***************************************************
c
c	DO Ibglaufm = 1, 100, 1
c
c	bglaufm = 0.1d0*dfloat(Ibglaufm)
c	bglau = bglaufm*5.07d0
c
c	CALL TA(espessura)
c	
c	WRITE(10,*) bglaufm, espessura*(5.07d0)**2.d0
c	WRITE(*,*) bglaufm, espessura*(5.07d0)**2.d0
c
c	ENDDO
c	END


	SUBROUTINE TA(espessura)	

        IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	PARAMETER (PI=3.141592)	

	EXTERNAL TAF

	ZMIN = -40.d0
	ZMAX = 40.d0
	EPS = 1.D-3

	espessura  = SGS1(ZMIN,ZMAX,EPS,TAF)

	RETURN
	END

C	***************************************************
C			 WOODS - SAXON 
C	***************************************************

	DOUBLE PRECISION FUNCTION TAF(Z)

        IMPLICIT DOUBLE PRECISION (A-H,K-Z)

	COMMON/bglau/ bglau
	
	r = DSQRT(bglau**2.d0+Z**2.d0)

	call woodsSaxon(r, densidade)

	TAF = densidade  	

	RETURN
	END




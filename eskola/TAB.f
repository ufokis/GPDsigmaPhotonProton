
C	ARQUIVOS NECESSÁRIOS

c	TAB.f
c	rho.f
c	dadmul.f

c	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
c	PARAMETER (PI=3.141592)
	
c	COMMON/B/ B

c	************** ARQUIVOS DE SAÍDA ******************
c	OPEN (UNIT=10,FILE='TAB208.dat',STATUS='unknown')	!WOODS-SAXON (WS)
c	***************************************************

c	DO IBfm = 1, 150, 1

c	Bfm = 0.1d0*dfloat(IBfm)
c	B = Bfm*5.07d0
	
c	CALL TAB(overNuclear)
	
c	WRITE(10,*) Bfm, overNuclear*(5.07d0)**2.d0
c	WRITE(*,*) Bfm, overNuclear*(5.07d0)**2.d0

c	ENDDO
c	END


	SUBROUTINE TAB(overlapNuclear)	

	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	INTEGER N,NFNEVL,MINPTS/1.D3/,MAXPTS/1.D6/,IWK/11.D5/
	DIMENSION A(4), AA(4), WK(1100000)
	PARAMETER (PI=3.141592)

	EXTERNAL TABF

	S1MIN=0.D0
	S1MAX=30.D0*5.07D0		

	TMIN=0.D0
	TMAX=PI
	
	Z1MIN = -50.D0*5.07D0
	Z1MAX =  50.D0*5.07D0

	Z2MIN = -50.D0*5.07D0
	Z2MAX =  50.D0*5.07D0
	
	N=4
	A(1) = S1MIN
        A(2) = TMIN
	A(3) = Z1MIN
	A(4) = Z2MIN

	AA(1) = S1MAX
	AA(2) = TMAX
	AA(3) = Z1MAX
	AA(4) = Z2MAX
	
        EPS=1.D-3

	CALL DADMUL(TABF,N,A,AA,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,
     *        NFNEVL,IFFAIL)

	overlapNuclear = RESULT

	RETURN
	END

C	***************************************************
C			 WOODS - SAXON 
C	***************************************************

	DOUBLE PRECISION FUNCTION TABF(N,XX)

        IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	INTEGER N
	DIMENSION XX(4)
	PARAMETER (PI=3.141592)	

	COMMON/B/ B

	S1 = XX(1)
	T  = XX(2)
	Z1 = XX(3)
	Z2 = XX(4)
	
	r1 = DSQRT(S1**2.0+Z1**2.D0)
	S2 = DSQRT(B**2.D0+S1**2.D0+2.D0*B*S1*DCOS(T))
	r2 = DSQRT(S2**2.0+Z2**2.D0)

	call woodsSaxon(r1, densidade1)
	call woodsSaxon(r2, densidade2)

	TABF = 2.d0*S1*densidade1*densidade2  	

	RETURN
	END




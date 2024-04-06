
c	fluxoFotons.f
c	formFactors.f
c	rho.f
c	dadmul.f
c	sgs.f
c	besselFunction.f	

c	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c	 Variando em b mantendo omega fixo
c	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
c	PARAMETER (PI=3.141592)	
c
c	COMMON/OMEGA/ OMEGA
c	COMMON/SNN/ SNN
c	COMMON/B1/ B1
c	COMMON/B/ B
c
c	EXTERNAL fluxobF
c
c	open (unit=12,file="fluxo_b.dat",status='unknown')
c
c	SNN = 2760.D0
c	OMEGA = 1.D0
c	
c	DO IB1 = 1, 300, 1
c
c	B1FM = 0.1d0*dfloat(IB1)
c	B1 = B1FM*5.07D0
c
c	CALL fluxoN0(flux0)
c	B = B1
c	CALL fluxoN2(flux2)
c
c	write(12,*) B1FM, flux0*5.07D0**2.D0, flux2*5.07D0**2.D0			!fator de correção: GeV --> 1/(GeV fm^2)	
c	write(*,*) B1FM, flux0*5.07D0**2.D0, flux2*5.07D0**2.D0
c
c	ENDDO
c	END

c	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c	Variando em omega e integrando em b
c	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
c	PARAMETER (PI=3.141592)	
c
c	COMMON/SNN/ SNN
c	COMMON/OMEGA/ OMEGA
c
c	EXTERNAL FBOMEGA
c
c	open (unit=12,file="fluxoFotons.dat",status='unknown')
c
c	SNN = 2760.D0
c
c	DO IOMEGA = 1, 1000, 5
c
c	OMEGA = 1.D-3*DBLE(IOMEGA)
c
c	RESULT = SGS1(39.4D0*5.07D0,68.3D0*5.07D0,1.D-3,FBOMEGA)
c	
c	write(12,*) OMEGA, RESULT		!fator de correção: GeV --> 1/(GeV fm^2)	
c	write(*,*) OMEGA, RESULT
c
c	ENDDO
c	END


C	DOUBLE PRECISION FUNCTION fluxobF(WW)	
C
C	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
C	PARAMETER (PI=3.141592)	
C
C	COMMON/OMEGA/ OMEGA
C
C	OMEGA = WW
C
CC	CALL fluxoN0(flux)
C	CALL fluxoN2(flux)
C
C	fluxobF = flux 
C
C	RETURN
C	END


c	DOUBLE PRECISION FUNCTION FBOMEGA(B1B1)	
c
c	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
c	PARAMETER (PI=3.141592)	
c
c	COMMON/B/ B
c	COMMON/B1/ B1
c
c	B1 = B1B1
c
cc	CALL fluxoN0(flux)
c	B = B1
c	CALL fluxoN2(flux)
c
c	FBOMEGA = flux
c
c	RETURN
c	END



	SUBROUTINE fluxoN0(flux)	

C	Essa subrotina retorna o fluxo padrão com dependência em b. Tal fluxo, aqui, está separado em Aux e Result,
C	onde aux são apenas parâmetros inicias da expressão do fluxo e Result é o quadrado do integrando
C	que está sendo integrado em kt (momento transverso do fóton).

        IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	PARAMETER (PI=3.141592)	

	COMMON/OMEGA/ OMEGA
	COMMON/SNN/ SNN
	COMMON/B1/ B1

	EXTERNAL fluxoN0F

	Z = 82.D0
	ALPHA = 1.D0/137.D0
	AUX = Z**2.D0*ALPHA/(PI**2.D0*OMEGA)

C	xxxx Utilizado quando usar fator de forma eletromagnético do tipo "Point Like" - CASO UPC xxxx
c	MASSAP = 0.938D0
c	GAM = SNN/(2.D0*MASSAP)
c	CHI = OMEGA*B1/GAM
c	RESULT = (CHI**2.D0/B1**2.D0)*(BESSK1(CHI))**2.D0
c	go to 30

	RESULT = (SGS1(0.D0,1.D0,1.D-3,fluxoN0F))**2.D0

   30	flux = AUX*RESULT
	
	RETURN
	END


	DOUBLE PRECISION FUNCTION fluxoN0F(KTKT)

C	Parte principal do fluxo, onde temos a dependência no fator de forma eletromagnético e
C	temos uma função de Bessel com dependência em b (nessa caso, B1)
	
        IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	PARAMETER (PI=3.141592)	

	COMMON/q/ q
	COMMON/B1/ B1
	COMMON/SNN/ SNN
	COMMON/OMEGA/ OMEGA

	KT = KTKT

	MASSAP = 0.938D0
	GAM = SNN/(2.D0*MASSAP)

	q = DSQRT((OMEGA/GAM)**2.D0+KT**2.D0)

	AUX = KT**2.0*BESSEL_JN(1,B1*KT)/(q**2.D0)

	call WoodsSaxon_Yukawa(formfactor)
c	call dipole(formfactor)

	fluxoN0F = AUX*formfactor

	RETURN
	END


c	************************************************
c	            Fluxo Fótons Efetivo
c	************************************************

	SUBROUTINE fluxoN2(flux)

        IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	INTEGER N,NFNEVL,MINPTS/1.D3/,MAXPTS/1.D6/,IWK/11.D5/
	DIMENSION A(2), AA(2), C(2), CC(2), WK(1100000)
	PARAMETER (PI=3.141592)	

	COMMON/RA/ RA
	COMMON/B/ B

	EXTERNAL F1
	EXTERNAL F2

	NMASSA = 207.D0
	RA = 1.2D0*NMASSA**(1.D0/3.D0)*5.07D0 		!RAIO DO PRÓTON (Hofstadter) = 0.74 fm ; RA = 1.2A^(1/3) fm; fm --> 5.07 GeV^-1

C	xxx Espectro Full xxx

	N=2
	A(1) = 0.D0		!B2MIN
	A(2) = 0.D0		!TMIN
	
	AA(1) = RA		!B2MAX
	AA(2) = PI		!TMAX
	
	EPS=1.D-3

	CALL DADMUL(F1,N,A,AA,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,
     *        NFNEVL,IFFAIL)

	FULL = 2.D0*RESULT	!somatório de todas as contribuições do flúxo na região do alvo, considerando a região de overlap nuclear

	IF (B > 2.D0*RA) GO TO 10 !usado no caso UPC

C	xxx Espectro Overlap xxx

	N=2
	C(1) = 0.D0
	C(2) = -1.D0
	
	CC(1) = 1.D0
	CC(2) = 1.D0
	
	CALL DADMUL(F2,N,C,CC,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,
     *        NFNEVL,IFFAIL)

	OVERLAP = RESULT	!somatório dos flúxos na região overlap

	AUX1 = 2.D0*RA**2.D0*DACOS(B/(2.D0*RA)) 
	AUX2 = (B/2.D0)*DSQRT(4.D0*RA**2.D0-B**2.D0)
	AREAEFF = PI*RA**2.D0-AUX1+AUX2

	flux = (FULL-2.D0*OVERLAP)/AREAEFF	!Descarta a região overlap, ficando apenas com a região dos espectadores
	
        GO TO 20

   10   CONTINUE
	
	AREAEFF = PI*RA**2.D0	!no caso UPC, a área efetiva será toda a região do núcleo

	flux = FULL/AREAEFF
	
   20   CONTINUE

	RETURN
	END	


	DOUBLE PRECISION FUNCTION F1(N,XX)

C	Acompanhar pelo artigo, olhando o desenho

        IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	DIMENSION XX(2)
	PARAMETER (PI=3.141592)

	COMMON/B/ B
	COMMON/B1/ B1
	COMMON/RA/ RA

	B2 = XX(1)
	T  = XX(2)

	B1 = DSQRT(B**2.D0+B2**2.D0+2.D0*B*B2*DCOS(T))

	CALL fluxoN0(flux)
C	******************************

	F1 = B2*flux
	
	RETURN
	END



	DOUBLE PRECISION FUNCTION F2(N,XX)

C	Acompanhar pelo artigo

        IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	DIMENSION XX(2)
	PARAMETER (PI=3.141592)

	COMMON/B/ B
	COMMON/B1/ B1
	COMMON/RA/ RA

	BY = DSQRT(RA**2.D0-B**2.D0/4.D0)*XX(1)
	BX = (DSQRT(RA**2.D0-BY**2.D0)-B/2.D0)*XX(2)+B/2.D0

	B1 = DSQRT(BX**2.D0+BY**2.D0)

	CALL fluxoN0(flux)
C	******************************

	J1 = DSQRT(RA**2.D0-B**2.D0/4.D0)		!JACOBIANOS
	J2 = DSQRT(RA**2.D0-BY**2.D0)-B/2.D0		!JACOBIANOS

	F2 = J1*J2*flux

	RETURN
	END





	SUBROUTINE fluxoFotonsProton(SNN, OMEGA, ESPECTRO)

	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	PARAMETER (PI=3.141592)

	ALPHA = 1.D0/137.D0
	Q02 = 0.71D0
	MASSAP = 0.938D0			!MASSA DO PRÓTON EM GEV
	GAM = SNN/(2.D0*MASSAP)
	
	QMIN2 = OMEGA**2.D0/GAM**2.D0
	
	CHI = 1.D0+Q02/QMIN2
	CHI2 = CHI**2.D0
	CHI3 = CHI**3.D0

	TER1 = ALPHA/(2.D0*PI*OMEGA)

	TERA = 1.D0+(1.D0-2.D0*OMEGA/SNN)**2.D0
	TERB = DLOG(CHI)-11.D0/6.D0+3.D0/CHI-1.5D0/CHI2+1.D0/(3.D0*CHI3)

	ESPECTRO = TER1*TERA*TERB

	RETURN
	END


	SUBROUTINE fluxoFotonsNucleo(SNN, OMEGA, ESPECTRO)

	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	PARAMETER (PI=3.141592)

	EXTERNAL BESSK0
	EXTERNAL BESSK1

	NMASSA = 208.D0
	ZATOM = 82.D0
	RA = 1.2D0*NMASSA**(1.D0/3.D0)*5.07D0 		! RA = 1.2A^(1/3) fm; fm --> 5.07 GeV^-1
	RH = 0.74D0*5.07D0				! RAIO DO PRÓTON (Hofstadter) = 0.74 fm
	MASSAP = 0.938D0				!MASSA DO PRÓTON EM GEV
	GAM = SNN/(2.D0*MASSAP)
	ALPHA = 1.D0/137.D0
	ETA = OMEGA*(2.D0*RA)/GAM
C	ETA = OMEGA*(RH+RA)/GAM
	
	TERMO1 = 2.D0*ZATOM**2.D0*ALPHA/(PI*OMEGA)
	TERMO2 = ETA*BESSK0(ETA)*BESSK1(ETA)
	
	UBESS2 = (BESSK1(ETA))**2.D0-(BESSK0(ETA))**2.D0
	TERMO3 = 0.5D0*ETA**2.D0*UBESS2

	ESPECTRO = TERMO1*(TERMO2-TERMO3)

	RETURN
	END


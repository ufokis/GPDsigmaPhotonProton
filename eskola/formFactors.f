
C	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
C	PARAMETER (PI=3.141592)
C
C	EXTERNAL WSYF
C
C	COMMON/q/ q
C
C	open (unit=12,file="formfactor.dat",status='unknown')
C
C	DO Iq = 1, 500
C
C	q = 0.001d0*dfloat(Iq)
C
C	CALL WoodsSaxon_Yukawa(formfactor)
C
C	aux = 4.D0*PI/(q*207.D0)
C	formfactor2 = aux*SGS1(0.D0,60.D0,1.D-3,WSYF)
C
C	WRITE(12,*) q, formfactor, formfactor2
C	WRITE(*,*) q, formfactor, formfactor2
C
C	ENDDO	
C	END 


c	DOUBLE PRECISION FUNCTION WSYF(rr)	
c
c	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
c	PARAMETER (PI=3.141592)	
c
c	COMMON/q/ q
c
c	call woodsSaxon(rr,densidade)		!Não esquecer de mudar tmin = 0.1
c
c	WSYF = densidade*dsin(q*rr)*rr
c
c	RETURN
c	END



C	**************************************************
C        	      WOODSSAXON + YUKAWA
C	**************************************************

	SUBROUTINE WoodsSaxon_Yukawa(formfactor)

	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	PARAMETER (PI=3.141592)

	COMMON/q/ q

	NMASSA = 207.D0
	RA = 1.2D0*NMASSA**(1.D0/3.D0)*5.07D0 		!RAIO DO PRÓTON (Hofstadter) = 0.74 fm ; RA = 1.2A^(1/3) fm; fm --> 5.07 GeV^-1

	a = 0.7d0*5.07d0				!fm -> 5.07 GeV^-1 
	a2 = a**2.d0

	RHO = 0.138D0/(5.07D0**3.D0)

	TER1 = 4.D0*PI*RHO/(NMASSA*q**3.D0)
	TER2 = DSIN(q*RA)-q*RA*DCOS(q*RA)
	TER3 = 1.D0/(1.D0+a2*q**2.d0)
	
	formfactor = TER1*TER2*TER3

	RETURN
	END


C	**************************************************
C        	        DIPOLE FORM FACTOR
C	**************************************************

	SUBROUTINE dipole(formfactor)

	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
	PARAMETER (PI=3.141592)

	COMMON/q/ q

	formfactor = 0.088D0**2.D0/(0.088D0**2.D0+q**2.D0)

	RETURN
	END



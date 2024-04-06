c	IMPLICIT DOUBLE PRECISION (A-H,J-Z)
c	PARAMETER (PI=3.141592)
c
c	open (unit=12,file="rho.dat",status='unknown')
c
c	DO irfm = 0, 100
c
c	rfm = 0.1d0*dfloat(irfm)
c	r = rfm*5.07d0
c
c	CALL woodsSaxon(r, densidade)
c	
c	densidade = densidade*5.07d0**3.d0	
c
c	WRITE(12,*) rfm, 4.d0*PI*densidade
c
c	ENDDO	
c	END 


C	***************************************************
C	        Distribuição Nuclear Woods-Saxon
C	***************************************************

	SUBROUTINE woodsSaxon(r, densidade)

	IMPLICIT DOUBLE PRECISION (A-H,J-Z)

	rho0  = 0.16d0*(5.07D0)**(-3.D0)
	c     = 6.624d0*5.07D0
	a     = 0.549d0*5.07D0

c	rho0  = 0.1693d0*(5.07D0)**(-3.D0)	!ouro
c	c     = 6.38d0*5.07D0			!ouro
c	a     = 0.535d0*5.07D0			!ouro

	densidade = rho0/(1.d0+dexp((r-c)/a))

	RETURN
	END


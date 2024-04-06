DOUBLE PRECISION FUNCTION SGS1A(x0,xn,n,f, nome)    
        
    IMPLICIT NONE
    character(9):: nome
    DOUBLE PRECISION, EXTERNAL::f
    INTEGER::i,n
    DOUBLE PRECISION::x0,xn,h,s
    !OPEN(UNIT=8, FILE=nome, STATUS='unknown')
    IF (MOD(n,2)==0) THEN
        h=(xn-x0)/n
        s=f(x0)+f(xn)+4*f(x0+h)
        DO i=3,n-1,2
        s=s+(4*f(x0+(i*h)))+(2*f(x0+(i-1)*h))
        !write(8,*) (x0+(i*h)), s
        END DO
        s=(h*s)/3
        PRINT *,"Value of integral is",s
        ELSE
        PRINT *,"Number of interval is not even"
    END IF
    !CLOSE(UNIT=8)
    SGS1A = s
END FUNCTION SGS1A   
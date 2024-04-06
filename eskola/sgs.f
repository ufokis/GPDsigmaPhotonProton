        DOUBLE PRECISION FUNCTION SGS1(A,B,EPS,F)
        DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
        DOUBLE PRECISION SP,SA,F,SGS9,EPS,abb
        EXTERNAL F
        S = 0.d0
        U = A
    1   V = B
        IF ( U .LT. B ) THEN
        SF = SGS9 ( F,U,V )
    2   C  = (U+V)/2
        SL = SGS9 ( F,U,C )
        SG = SGS9 ( F,C,V )
        SP = SL+SG
        abb=abs(sf)
        if(abb.ne.0.0) go to 5
        abb=1.
    5   SA = ABS(SP-SF)/(abb*EPS)
        IF (  SA.GE.1.0 ) THEN
        V  = C
        SF = SL
        GOTO 2
        END IF
        U = V
        S = S+SP
        GOTO 1
        END IF
        SGS1=S
        RETURN
        END
C
        DOUBLE PRECISION FUNCTION SGS9 ( F,A,B )
        DOUBLE PRECISION A,B,H,S,C,X,F
        EXTERNAL F
        C = (A+B)/2
        H = (B-A)/2
        X = .96028985E0*H
        S = .10122853E0*(F(C+X)+F(C-X))
        X = .79666647E0*H
        S = S + .22238103E0*(F(C+X)+F(C-X))
        X = .52553240E0*H
        S = S + .31370664E0*(F(C+X)+F(C-X))
        X = .18343464E0*H
        S = S + .36268378E0*(F(C+X)+F(C-X))
        SGS9 = S * H
        RETURN
        END


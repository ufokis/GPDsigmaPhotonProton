program check
    use GPDsigmaPhotonProton
    implicit none
    complex(wp):: A
    integer     :: i
    real(wp)    :: x
    xi = 0.1_wp
    open(unit = 2, file = 'dados_test.dat', status = 'unknown')
    do i = 0, 999
        x = (float(i)+1.d-7)/1000.0_wp
        A = IntComplex(x)
        write(2, *) x, real(A), aimag(A) 
    end do    
end program check

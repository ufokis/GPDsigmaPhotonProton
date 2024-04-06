program calculo_pdf

    real(8)   :: x, QA
    complex(8):: ResultIntComplex
    character:: name*40
    real(8)   :: PDF(-6:6), alphasPDF, I1, I2

    name = 'nCTEQ15_208_82'
    call InitPDFsetByName(name)


    call InitPDF(0)

	call GetXmin(0,xmin)
    call GetXmax(0,xmax)
    call GetQ2min(0,q2min)
    call GetQ2max(0,q2max)
    print *,'xmin=',xmin,' xmax=',xmax,' Q2min=',q2min,' Q2max=',q2max
    call GetMinMax(0,xmin,xmax,q2min,q2max)
    print *,'xmin=',xmin,' xmax=',xmax,' Q2min=',q2min,' Q2max=',q2max
    
	QA = 1.3d0
    alpha_s = alphasPDF(2.37d0)
    ! write(*,*) 'aqui', alpha_s
    open(unit = 2, file='pdf1.dat', status='unknown')
    open(unit = 3, file='pdf2.dat', status='unknown')
    do x = 5.d-6, 0.95d0, 1.d-6
        call evolvePDF(x, QA, PDF)     
        I1 = PDF(0)
        I2 = (sum(PDF(1:4)) + sum(PDF(-1:-4)))
        write(2,*) x, I1
        write(3,*) x, I2*1.d-2

    end do        
    

end program calculo_pdf

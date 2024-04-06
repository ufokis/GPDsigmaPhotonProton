program Test_PDFs
    character(64) :: name
    double precision f(-6:6), xmin, xmax, q2min, q2max, x, Q, alphasPDF, I1, I2
    double precision :: step
    integer :: num_steps
    character(20) :: arg
    integer :: i

    ! Inicializando os argumentos
    Q = 0.0d0
    name = ''

    ! Analisando os argumentos da linha de comando
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        select case (arg)
            case ("-Q")
                call get_command_argument(i+1, arg)
                read(arg, *) Q
            case ("-name")
                call get_command_argument(i+1, name)
        end select
    end do

    ! Verificando se os argumentos necessários foram fornecidos
    if (Q == 0.0d0 .or. name == '') then
        print *, "Erro: Por favor, forneça o valor de Q e o nome do conjunto de dados."
        stop
    end if

    call InitPDFsetByName(name)

    call InitPDF(0)
    call GetXmin(0,xmin)
    call GetXmax(0,xmax)
    
    print *,'xmin=',xmin,' xmax=',xmax,' Q2min=',q2min,' Q2max=',q2max
    call GetMinMax(0,xmin,xmax,q2min,q2max)
    print *,'xmin=',xmin,' xmax=',xmax,' Q2min=',q2min,' Q2max=',q2max

    ! Adaptação para x variar de xmin a xmax com passo constante
    
    step = 0.0001d0  ! Defina aqui o tamanho do passo desejado
    num_steps = nint((xmax - 1.d-7) / step)

    open(unit = 2, file='pdf1_' // trim(name) // '.dat', status='unknown')
    open(unit = 3, file='pdf2_' // trim(name) // '.dat', status='unknown')
   
    do ix = 0, num_steps
        x = 1.d-7 + real(ix) * step

        call evolvePDF(x, Q, f)

        I1 = f(0)
        I2 = (sum(f(1:4)) + sum(f(-1:-4)))/x
        
        write(2,*) x, I1
        write(3,*) x, I2*1.d-2
    end do

end program Test_PDFs

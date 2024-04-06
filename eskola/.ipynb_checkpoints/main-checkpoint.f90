program lhapdf_example
    implicit none

    ! Declarações de variáveis
    integer :: ierr
    real(kind=8) :: x, Q, pdf_gluon, pdf_strange

    ! Defina o ponto x e a escala de energia Q
    x = 0.01d0
    Q = 10.d0

    ! Inicialize o LHAPDF
    call lhapdf_init("NNPDF31_nnlo_as_0118", ierr)

    ! Verifique se ocorreu um erro na inicialização
    if (ierr /= 0) then
        print *, "Erro ao inicializar LHAPDF"
        stop
    end if

    ! Obtenha as PDFs para o gluon e o quark s
    pdf_gluon = lhapdf_xfx(0, x, Q, ierr) / x
    pdf_strange = lhapdf_xfx(3, x, Q, ierr) / x

    ! Verifique se ocorreu um erro ao obter as PDFs
    if (ierr /= 0) then
        print *, "Erro ao obter PDFs"
        stop
    end if

    ! Imprima as PDFs obtidas
    print *, "PDF do gluon em x =", x, "e Q =", Q, ":", pdf_gluon
    print *, "PDF do quark s em x =", x, "e Q =", Q, ":", pdf_strange

    ! Libere recursos do LHAPDF
    call lhapdf_cleanup()

end program lhapdf_example
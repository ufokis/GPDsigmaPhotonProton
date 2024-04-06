module SGS1_lib

    interface
        double precision function SGS2(A, B, EPS, F)
            double precision:: A, B, EPS, F
            external F
        end function
    end interface
    
    double precision::resultado

    contains

        double precision function quadratura_gaussiana(a, b, f)
            double precision, intent(in) :: a, b
            double precision:: f
            external f
            double precision :: xg(8), wg(8)
            double precision :: integral
            integer :: i

             ! Pontos e pesos de quadratura gaussiana (8 pontos)
            xg = (/ -0.1834346424956498d0, -0.5255324099163290d0, -0.7966664774136267d0, &
                    0.0d0, 0.1834346424956498d0, 0.5255324099163290d0, 0.7966664774136267d0, &
                    0.9602898564975363d0 /)

            wg = (/ 0.3626837833783620d0, 0.3137066458778873d0, 0.2223810344533745d0, &
                    0.1012285362903763d0, 0.3626837833783620d0, 0.3137066458778873d0, &
                    0.2223810344533745d0, 0.1012285362903763d0 /)

            ! Inicialize a variável de integração
            integral = 0.0

            ! Realize a quadratura gaussiana
            do i = 1, 8
            integral = integral + wg(i) * f(0.5 * (b - a) * xg(i) + 0.5 * (a + b))
            end do

            ! Ajuste para o intervalo [a, b]
            integral = 0.5 * (b - a) * integral

            ! Retorne o resultado da integração
            quadratura_gaussiana = integral
        end function quadratura_gaussiana


        recursive double precision function sgs3(a, b, eps, f) result(resultado)
            double precision, intent(in) :: a, b, eps
            double precision:: f
            external f
            double precision :: integral1, integral2, m

            ! Calcula a integral com a quadratura gaussiana em [a, b]
            integral1 = quadratura_gaussiana(a, b, f)

            ! Calcula a integral com a quadratura gaussiana em [a, m] e [m, b]
            m = 0.5 * (a + b)
            
            integral2 = quadratura_gaussiana(a, m, f) + quadratura_gaussiana(m, b, f)

            ! Verifica a tolerância e refinamento adaptativo
            if (abs(integral1 - integral2) < eps) then
                resultado = integral2
            else
                ! Refina recursivamente os intervalos
                resultado = sgs3(a, m, eps, f) + sgs3(m, b, eps, f)
            end if
        end function sgs3

end module SGS1_lib   
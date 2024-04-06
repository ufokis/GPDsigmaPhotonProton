#!/bin/bash

# Limpar a tela
clear

# Definir variáveis para os arquivos de código-fonte
SRC_FILES="teste.f90 sgs.f fast_clog.f90  besselFunction.f dadmul.f formFactors.f TA.f rho.f sigAAinelastico.f TAB.f fluxoFotons.f"

# Nome do executável de saída
OUTPUT_EXECUTABLE="sigmaAA.out"

# Nome do arquivo de log do Valgrind
VALGRIND_LOG="valgrind_log.txt"

# Compilar usando gfortran com opções para o Valgrind
echo "Compilando o código-fonte para análise de memória ou não"

gfortran -c sgs.f
gfortran -c besselFunction.f
gfortran -c dadmul.f 
gfortran -c formFactors.f
gfortran -c TA.f
gfortran -c rho.f
gfortran -c sigAAinelastico.f
gfortran -c TAB.f
gfortran -c fluxoFotons.f
gfortran -c fast_clog.f90
gfortran -c SGS1_lib.f90
gfortran -c LI2.f90


gfortran -c sgs2.f
gfortran -c SGS1_lib.f90
#gfortran -c nCTEQ-code/Cteq6Pdf.f
#gfortran -c nCTEQ-code/setNpdf2.f
#gfortran -c nCTEQ-code/Cteq6Pdfab.f
#gfortran -c nCETEQ15_77.f90 
gfortran -c sigmaPhotonProton.f90

gfortran -c dsidyAAupc_Eskola.f90

gfortran -g -O0  fast_clog.o SGS1_lib.o LI2.o sgs2.o sgs.o  besselFunction.o dadmul.o formFactors.o TA.o rho.o sigAAinelastico.o TAB.o fluxoFotons.o sigmaPhotonProton.o dsidyAAupc_Eskola.o -o $OUTPUT_EXECUTABLE -L/home/eliton/lhapdf2/lib -lLHAPDF

# Verificar se a compilação foi bem-sucedida
if [ $? -eq 0 ]; then
    echo "Compilação concluída com sucesso."

    echo "Executando $OUTPUT_EXECUTABLE Normal"
    
    #./sigmaAA.out

    # Executar o programa compilado com Valgrind
    #echo "Executando $OUTPUT_EXECUTABLE com Valgrind..."
    valgrind --leak-check=full  --show-leak-kinds=all --log-file=$VALGRIND_LOG --track-origins=yes ./sigmaAA.out

    # Exibir o resultado do Valgrind
    #cat $VALGRIND_LOG

    #echo "Executando $OUTPUT_EXECUTABLE GDB"
    #gdb ./sigmaA.out
else
    echo "Erro durante a compilação. Verifique os arquivos de código-fonte."
fi

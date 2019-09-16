CC = g++

CFLAGS =  -c -Wall `root-config --cflags`

CLIBS = `root-config --glibs` -lboost_system -lboost_filesystem `root-config --cflags`

all: PALS_avalanche

PALS_avalanche: PALS_avalanche.cpp PALS_avalanche_lib.cpp PALS_avalanche_fileTools.cpp PALS_avalanche_fitFunctions.cpp PALS_avalanche_resultsSaver.cpp PALS_avalanche_libTools.cpp PALS_avalanche_deconvolution.o
	$(CC) PALS_avalanche.cpp PALS_avalanche_lib.cpp PALS_avalanche_fileTools.cpp PALS_avalanche_fitFunctions.cpp PALS_avalanche_resultsSaver.cpp PALS_avalanche_libTools.cpp PALS_avalanche_deconvolution.o -o PALS_avalanche_pr0.exe $(CLIBS) -std=c++11 -O2 -larmadillo -llapack -lblas

PALS_avalanche_deconvolution.o: PALS_avalanche_deconvolution.cpp PALS_avalanche_libTools.cpp
	$(CC) $(CFLAGS) PALS_avalanche_deconvolution.cpp PALS_avalanche_libTools.cpp -std=c++11 -O2 -larmadillo -llapack -lblas
	
clean: 
	rm -rf *o
	rm -rf *.root
	rm -rf Results/
	rm -rf *.exe

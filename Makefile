CC = g++

CFLAGS = -c -Wall `root-config --cflags`

CLIBS = `root-config --glibs` -lboost_system -lboost_filesystem

all: PALS_avalanche

PALS_avalanche: PALS_avalanche.o PALS_avalanche_lib.o
	$(CC) PALS_avalanche.o PALS_avalanche_lib.o -o PALS_avalanche.exe $(CLIBS) -std=c++11 -O2 -larmadillo -llapack -lblas

PALS_avalanche.o: PALS_avalanche.cpp
	$(CC) $(CFLAGS) PALS_avalanche.cpp -std=c++11 -O2 -larmadillo -llapack -lblas

PALS_avalanche_lib.o: PALS_avalanche_lib.cpp
	$(CC) $(CFLAGS) PALS_avalanche_lib.cpp -std=c++11 -O2 -larmadillo -llapack -lblas
	
			
clean: 
	rm -rf *o 
	rm -rf *.png
	rm -rf *.pdf
	rm -rf *.root

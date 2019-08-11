CC = g++

CFLAGS = -Wall `root-config --cflags`

CLIBS = `root-config --glibs` -lboost_system -lboost_filesystem `root-config --cflags`

all: PALS_avalanche

PALS_avalanche: PALS_avalanche.cpp PALS_avalanche_lib.cpp PALS_avalanche_fileTools.cpp PALS_avalanche_fitFunctions.cpp PALS_avalanche_resultsSaver.cpp PALS_avalanche_libTools.cpp
	$(CC) PALS_avalanche.cpp PALS_avalanche_lib.cpp PALS_avalanche_fileTools.cpp PALS_avalanche_fitFunctions.cpp PALS_avalanche_resultsSaver.cpp PALS_avalanche_libTools.cpp -o PALS_avalanche.exe $(CLIBS) -std=c++11

clean: 
	rm -rf *o
	rm -rf *.root
	rm -rf Results/

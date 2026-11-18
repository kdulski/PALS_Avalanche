CXX      ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra
# ROOT flags come last so that ROOT's -std=c++XX (>= c++17) wins
CPPFLAGS += -Isrc -Ithird_party
ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs)
ROOTGLIBS   = $(shell root-config --glibs)

CORE  = src/Config.cpp src/Convolution.cpp src/Distribution.cpp src/FormulaModel.cpp src/LifetimeModel.cpp \
        src/ModelSetup.cpp src/Nnls.cpp src/Spectrum.cpp src/TextReport.cpp
ROOTS = src/FitSession.cpp src/Fitter.cpp src/RootFormula.cpp src/RootIO.cpp
GUI   = src/gui/PalsGui.cpp src/gui/main_gui.cpp
EXE   = PALS_avalanche

all: $(EXE) gui

$(EXE): $(CORE:.cpp=.o) $(ROOTS:.cpp=.o) src/main.o
	$(CXX) $^ -o $@ $(ROOTLIBS)

gui: PALS_avalanche_gui

PALS_avalanche_gui: $(CORE:.cpp=.o) $(ROOTS:.cpp=.o) $(GUI:.cpp=.o)
	$(CXX) $^ -o $@ $(ROOTGLIBS)

src/gui/%.o: src/gui/%.cpp src/gui/*.h src/*.h
	$(CXX) -std=c++17 $(CXXFLAGS) $(CPPFLAGS) -Isrc/gui $(ROOTCFLAGS) -c $< -o $@

src/%.o: src/%.cpp src/*.h
	$(CXX) -std=c++17 $(CXXFLAGS) $(CPPFLAGS) $(ROOTCFLAGS) -c $< -o $@

# Tests and benchmark of the model; they do not need ROOT
tests/test_model: tests/test_model.cpp $(CORE) src/*.h tests/*.h
	$(CXX) -std=c++17 -O2 -Wall $(CPPFLAGS) tests/test_model.cpp $(CORE) -o $@

test: tests/test_model
	./tests/test_model ColonCancer_AB_20122018.dat_TDiff.dat config/tumor.json

clean:
	rm -f src/*.o src/gui/*.o $(EXE) PALS_avalanche_gui tests/test_model

.PHONY: all gui test clean

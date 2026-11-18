// PALS_avalanche - discrete fit of positron annihilation lifetime spectra
//
// Usage: PALS_avalanche <data file> [config.json]
//   data file   - time differences / histogram (text) or ROOT file
//   config.json - fit configuration (default: FitDetails.json)
// The graphical interface is PALS_avalanche_gui.

#include "Config.h"
#include "FitSession.h"
#include "Spectrum.h"

#include <TH1.h>
#include <TROOT.h>

#include <chrono>
#include <iostream>

using namespace pals;

int main(int argc, char* argv[])
{
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <data file> [config.json]\n";
        return 2;
    }
    const std::string dataPath = argv[1];
    const std::string configPath = argc > 2 ? argv[2] : "FitDetails.json";
    const auto start = std::chrono::steady_clock::now();

    try {
        gROOT->SetBatch(kTRUE);
        TH1::AddDirectory(kFALSE);

        FitSession session;
        session.setConfig(Config::load(configPath), configPath);
        session.loadData(dataPath);

        std::cout << "\n------------- Preliminary estimates -------------\n";
        printEstimates(std::cout, session.spectrum(), session.estimate());

        session.fit();
        std::cout << "\n";
        session.writeReport(std::cout);
        const std::string base = session.defaultOutputBase();
        session.save(base);
        std::cout << "\nResults: " << base << ".txt / .root\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Elapsed time: " << ms.count() << " ms\n";
    return 0;
}

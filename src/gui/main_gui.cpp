// PALS_avalanche_gui - graphical interface
//
// Usage: PALS_avalanche_gui [data file] [config.json]

#include "PalsGui.h"

#include <TApplication.h>
#include <TGClient.h>
#include <TH1.h>
#include <TStyle.h>

#include <exception>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
    TApplication app("PALS_avalanche_gui", &argc, argv);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);

    auto* gui = new pals::PalsGui(gClient->GetRoot(), 1500, 950);
    try {
        if (app.Argc() > 2) gui->openConfig(app.Argv(2));
        if (app.Argc() > 1) gui->openData(app.Argv(1));
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
    }
    app.Run();
    return 0;
}

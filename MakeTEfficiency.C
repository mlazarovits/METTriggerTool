#include <vector>
#include <iostream>
#include <filesystem>
#include "TEfficiency.h"

using std::cout;
using std::endl;
namespace fs = std::filesystem;

#include <TChain.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1D.h>

#include <iostream>
#include <string>


#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

std::vector<std::string> getEOSFiles(const std::string& eosPath)
{
    std::vector<std::string> files;

    std::string command =
        "xrdfs root://cmseos.fnal.gov ls " + eosPath;

    FILE* pipe = popen(command.c_str(), "r");

    if (!pipe) {
        std::cerr << "ERROR: Could not run xrdfs" << std::endl;
        return files;
    }

    char buffer[4096];

    while (fgets(buffer, sizeof(buffer), pipe)) {

        std::string file(buffer);

        // Remove newline
        if (!file.empty() && file.back() == '\n')
            file.pop_back();

        // Only keep ROOT files
        if (file.size() >= 5 &&
            file.substr(file.size() - 5) == ".root") {

            files.push_back(file);
        }
    }

    pclose(pipe);

    return files;
}


// ------------------------------------------------------------
// Configuration
// ------------------------------------------------------------

const int NBINS = 20;
const double MET_MIN = 0.0;
const double MET_MAX = 150.0;

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------

void makeTriggerEfficiency(TChain* chain, string triggername)
{
    if (!chain) {
        std::cerr << "ERROR: Null TChain!" << std::endl;
        return;
    }

    // --------------------------------------------------------
    // Set up branches
    // --------------------------------------------------------

    float MET;
    bool trigger;

    chain->SetBranchAddress("Met_pt", &MET);
    chain->SetBranchAddress(triggername.c_str(), &trigger);

    // --------------------------------------------------------
    // Create TEfficiency
    // --------------------------------------------------------

    TEfficiency* eff = new TEfficiency(
        "triggerEfficiency",
        "Trigger Efficiency;MET [GeV];Efficiency",
        NBINS,
        MET_MIN,
        MET_MAX
    );

    // Optional: choose how the uncertainty intervals are calculated
    eff->SetStatisticOption(TEfficiency::kFCP);

    // --------------------------------------------------------
    // Loop over events
    // --------------------------------------------------------

    Long64_t nEntries = 100000;//chain->GetEntries();

    std::cout << "Processing " << nEntries << " events..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {

        chain->GetEntry(i);

        // ----------------------------------------------------
        // Denominator selection
        // ----------------------------------------------------

        // Put your denominator/preselection here.
        //
        // Example:
        //
        // if (MET < 50)
        //     continue;

        bool denominator = true;

        if (!denominator)
            continue;

        // ----------------------------------------------------
        // Fill efficiency
        // ----------------------------------------------------

        eff->Fill(trigger, MET);
    }

    // --------------------------------------------------------
    // Save to ROOT file
    // --------------------------------------------------------

    string fout = "triggerEfficiency.root";
    TFile* output = TFile::Open(fout.c_str(), "RECREATE");
    if (!output || output->IsZombie()) {
        std::cerr << "ERROR: Could not create output file!"
                  << std::endl;
        return;
    }

    eff->Write();

    output->Close();

    std::cout << "Saved efficiency to " << fout
              << std::endl;

    delete eff;
}


void MakeTEfficiency(){

	
	TChain* chain = new TChain("tree/llpgtree");

	std::string path =
	    "/store/user/lpcsusylep/jaking/KUCMSNtuple/";
	string ntupledir = "kucmsntuple_EGamma_R18_InvMetPho30_noSV_v31/EGamma/kucmsntuple_EGamma_R18_InvMetPho30_noSV_v31_EGamma_MINIAOD_Run2018C-12Nov2019_UL2018-v2/260220_193726/0000/";
	
	std::vector<std::string> files = getEOSFiles(path+ntupledir);
	cout << "Looping over " << files.size() << " files" << endl;	
	for (const auto& file : files) {
	
	    std::string xrootdFile =
	        "root://cmseos.fnal.gov/" + file;
	
	    //std::cout << "Adding: " << xrootdFile << std::endl;
	
	    chain->Add(xrootdFile.c_str());
	}
	
	std::cout << "Total files: "
	          << files.size()
	          << std::endl;
	
	makeTriggerEfficiency(chain,"HLT_PFMET120_PFMHT120_IDTight_v");


}

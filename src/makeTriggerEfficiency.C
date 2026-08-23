#include <filesystem>
#include <TChain.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1D.h>
#include <map>
#include <algorithm>
#include <iostream>
#include <string>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
namespace fs = std::filesystem;



std::vector<std::string> getEOSFiles(const std::string& eosPath)
{
    std::vector<std::string> files;

    std::string command =
        "xrdfs root://cmseos.fnal.gov ls -R " + eosPath;

    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        std::cerr << "ERROR: Could not run xrdfs" << std::endl;
        return files;
    }

    char buffer[4096];

    while (fgets(buffer, sizeof(buffer), pipe)) {
        std::string path(buffer);

        // Remove newline
        if (!path.empty() && path.back() == '\n')
            path.pop_back();

        // Only keep ROOT files
        if (path.size() >= 5 &&
            path.substr(path.size() - 5) == ".root") {
            files.push_back(path);
        }
    }

    pclose(pipe);

    return files;
}

// ------------------------------------------------------------
// Configuration
// ------------------------------------------------------------

const int NBINS = 25;
const double MET_MIN = 0.0;
const double MET_MAX = 250.;

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------

TEfficiency* makeTriggerEfficiency(TChain* chain, string triggername, string htbin, string teff_name)
{
    if (!chain) {
        std::cerr << "ERROR: Null TChain!" << std::endl;
        return nullptr;
    }

    // --------------------------------------------------------
    // Set up branches
    // --------------------------------------------------------

    float MET;
    vector<float>* jetpt = nullptr;
    bool trigger;

    chain->SetBranchAddress("Met_CPt", &MET); //corrected met
    chain->SetBranchAddress("Jet_pt", &jetpt);
    chain->SetBranchAddress(triggername.c_str(), &trigger);
   
   
    float hthi = std::stof(string(htbin.substr(htbin.find("_")+1)));
    float htlo = std::stof(string(htbin.substr(0,htbin.find("_"))));


    // --------------------------------------------------------
    // Create TEfficiency
    // --------------------------------------------------------

    TEfficiency* eff = new TEfficiency(
        teff_name.c_str(),
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

    Long64_t nEntries = chain->GetEntries();

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
        //HT calculation placeholder until can run on skims
        double ht = 0;
        for(int j = 0; j < (int)jetpt->size(); j++){
		ht += jetpt->at(j);
	}
	if(ht < htlo || ht >= hthi)
		continue;

        bool denominator = true;

        if (!denominator)
            continue;

        // ----------------------------------------------------
        // Fill efficiency
        // ----------------------------------------------------
        eff->Fill(trigger, MET);
    }

    return eff;

}


string makeTEffName(string year, string trigger, string htbin){
	std::erase(trigger,'_');
	return trigger+"_"+year+"_"+htbin;
}

int main(){

	//file dict
	map<string, vector<string>> filedirs = {};
	filedirs["18"] = {"kucmsntuple_EGamma_R18_SVHPM100_v51/EGamma/kucmsntuple_EGamma_R18_SVHPM100_v51_EGamma_MINIAOD_Run2018C-12Nov2019_UL2018-v2/"};
	filedirs["24"] = {"kucmsntuple_EGamma_R24_SVHPM100_MiniAOD_v34/EGamma0/kucmsntuple_EGamma_R24_SVHPM100_MiniAOD_v34_EGamma0_MINIAOD_Run2024C-MINIv6NANOv15-v1/"};
	

	std::string path =
	    "/store/user/lpcsusylep/jaking/KUCMSNtuple/";
	vector<string> years = {"18","24"};
	string trigger = "HLT_Ele35_WPTight_Gsf_v9";
	string htbin = "500_700";
	vector<TEfficiency*> effcurves;
	for(auto year : years){
		TChain* chain = new TChain("tree/llpgtree");
		vector<string> ntupledirs = filedirs[year];
		std::vector<string> files;
		for(auto ntupledir : ntupledirs){
			std::vector<std::string> thesefiles = getEOSFiles(path+ntupledir);
			files.insert(files.end(), thesefiles.begin(), thesefiles.end());
		}
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
		string teff_name = makeTEffName(year, trigger, htbin);			
		TEfficiency* eff = makeTriggerEfficiency(chain,trigger,htbin,teff_name);
		cout << "Passed:   " << eff->GetPassedHistogram()->GetEntries() << endl;
cout << "Total:    " << eff->GetTotalHistogram()->GetEntries() << endl;
		if(eff->GetTotalHistogram()->GetEntries() == 0) cout << "Empty eff " << teff_name << endl;
		if(eff == nullptr) continue;
		effcurves.push_back(eff);
		break; //debug
	}
    	string fout = "triggerEfficiency.root";
    	TFile* output = TFile::Open(fout.c_str(), "RECREATE");
    	if (!output || output->IsZombie()) {
    	    std::cerr << "ERROR: Could not create output file!"
    	              << std::endl;
    	    return -1;
    	}
	output->cd();
	for(auto eff : effcurves)
		eff->Write();
    	output->Close();

    	std::cout << "Saved efficiency to " << fout
    	          << std::endl;
	return 0;
}

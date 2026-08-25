#include <filesystem>
#include <TChain.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <TH1D.h>
#include <map>
#include <algorithm>
#include <iostream>
#include <string>
#include <cstdio>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::pair;
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

const int NBINS = 30;
const double MET_MIN = 0.0;
const double MET_MAX = 300.;

string makeTriggerCutstring(const vector<string>& trigs){
	string trigsel;
	for(int t = 0; t < (int)trigs.size(); t++){
	    trigsel += trigs[t] + " || ";
	}
	trigsel = trigsel.substr(0,trigsel.size()-4);
	trigsel = "("+trigsel+")";
	return trigsel;
}

string makeHTBinCutstring(const string& htbin){
    std::size_t underscore = htbin.find('_');

    if (underscore == std::string::npos)
        return "";

    std::string lower = htbin.substr(0, underscore);
    std::string upper = htbin.substr(underscore + 1);

    if (!lower.empty() && !upper.empty())
        return "((ht > " + lower + ") && (ht < " + upper + "))";

    if (!lower.empty())
        return "(ht > " + lower + ")";

    if (!upper.empty())
        return "(ht < " + upper + ")";

    return "";
}

void makeTriggerEfficiency(TChain* chain, const vector<string>& triggers, const vector<string>& preseltriggers, vector<string> htbins, string year, vector<TEfficiency*>& effs)
{
	if (!chain) {
	    std::cerr << "ERROR: Null TChain!" << std::endl;
	    return;
	}
	//do rdataframe hists for tefficiency - add ht binning here!
	string trigsel = makeTriggerCutstring(triggers);
	string trigpresel = makeTriggerCutstring(preseltriggers);

	string evtfilters = "((Flag_BadChargedCandidateFilter) && (Flag_BadPFMuonFilter) && (Flag_EcalDeadCellTriggerPrimitiveFilter) && (Flag_HBHENoiseFilter) && (Flag_HBHENoiseIsoFilter) && (Flag_ecalBadCalibFilter) && (Flag_eeBadScFilter) && (Flag_goodVertices))";//"(Flag_MetFilters == 1)"; //need pts cut when running over skims


	ROOT::RDataFrame df(*chain);
	auto df0 = df.Define("ht","ROOT::VecOps::Sum(Jet_pt)");

	map<string,pair<ROOT::RDF::RResultPtr<TH1D>, ROOT::RDF::RResultPtr<TH1D>>> dfhists;
	for(auto htbin : htbins){
		string htcutstring = makeHTBinCutstring(htbin); 
		auto df_htbin = df0.Filter(htcutstring);
		//total hist (denom)
		TH1D htotal_model("hTotal","PFMETOR_Efficiency;MET [GeV];Efficiency",NBINS,MET_MIN,MET_MAX);
		TH1D hpass_model("hPass","PFMETOR_Efficiency;MET [GeV];Efficiency",NBINS,MET_MIN,MET_MAX);
		string presel = trigpresel+" && "+evtfilters + " && "+htcutstring; //in addition to ntuple "filter" of (>= 1 SV || >= 1 photon[pt > 30]) && MET > 100
		
		auto hTotal = df0.Filter(presel).Histo1D(htotal_model,"Met_CPt");	
		//pass hist (num)
		auto hPass = df0.Filter(presel+" && "+trigsel).Histo1D(hpass_model,"Met_CPt");

		dfhists[htbin] = std::make_pair(hTotal, hPass);
	}
	df0.Report()->Print();
	for(auto it = dfhists.begin(); it != dfhists.end(); it++){
		string htbin = it->first;
		TEfficiency* eff = new TEfficiency(*(it->second.second),*(it->second.first));
		cout << "htbin " << htbin << endl;
		cout << "Passed:   " << eff->GetPassedHistogram()->GetEntries() << endl;
		cout << "Total:    " << eff->GetTotalHistogram()->GetEntries() << endl;
		if(eff->GetTotalHistogram()->GetEntries() == 0) cout << "Empty eff " << year << endl;
		if(eff == nullptr) continue;

		// Optional: choose how the uncertainty intervals are calculated
		eff->SetStatisticOption(TEfficiency::kFCP);
		eff->SetName(("PFMETOR_Efficiency_"+year+"_HTBin"+htbin).c_str());
		effs.push_back(eff);
	}

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
	vector<string> preseltriggers = {"HLT_Ele35_WPTight_Gsf_v9","HLT_Photon20_v","HLT_Mu55_v3","HLT_Mu12_v3","HLT_IsoMu27_v16","HLT_IsoMu20_v15","HLT_Ele27_WPTight_Gsf_v16"};
	vector<string> triggers = {"HLT_PFMET120_PFMHT120_IDTight_v","HLT_PFMET120_PFMHT120_IDTight_PFHT60_v","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v"};
	vector<string> htbins = {"_500","500_700","700_1000","1000_"};
	vector<TEfficiency*> effcurves;
	for(auto it = filedirs.begin(); it != filedirs.end(); it++){
		string year = it->first;
		cout << "Doing year " << year << endl;
		TChain* chain = new TChain("tree/llpgtree");
		vector<string> ntupledirs = it->second;
		std::vector<string> files;
		for(auto ntupledir : ntupledirs){
			std::vector<std::string> thesefiles = getEOSFiles(path+ntupledir);
			files.insert(files.end(), thesefiles.begin(), thesefiles.end());
		}
		cout << "Looping over " << files.size() << " files" << endl;
		for (const auto& file : files) {
			std::string xrootdFile = "root://cmseos.fnal.gov/" + file;		
		    //std::cout << "Adding: " << xrootdFile << std::endl;
		
		    chain->Add(xrootdFile.c_str());
		}
			
		std::cout << "Total files: "
		          << files.size()
		          << std::endl;
		makeTriggerEfficiency(chain,triggers,preseltriggers,htbins,year,effcurves);
		break; //2018 only while 2024 lepton triggers are sorted out
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

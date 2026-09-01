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

const int NBINS = 30;
const double MET_MIN = 0.0;
const double MET_MAX = 600.;

void makeTriggerEfficiency(string skimfile, const vector<string>& triggers, const vector<string>& preseltriggers, vector<string> htbins, string year, vector<TEfficiency*>& effs)
{
	//do rdataframe hists for tefficiency - add ht binning here!
	string trigsel = makeTriggerCutstring(triggers);
	string trigpresel = makeTriggerCutstring(preseltriggers);

	string evtfilters = "((Flag_BadChargedCandidateFilter) && (Flag_BadPFMuonFilter) && (Flag_EcalDeadCellTriggerPrimitiveFilter) && (Flag_HBHENoiseFilter) && (Flag_HBHENoiseIsoFilter) && (Flag_ecalBadCalibFilter) && (Flag_eeBadScFilter) && (Flag_goodVertices))";//"(Flag_MetFilters == 1)"; //need pts cut when running over skims


	ROOT::RDataFrame df("kuSkimTree",skimfile);
	auto df0 = df.Define("ht","ROOT::VecOps::Sum(selJetPt)");

	map<string,pair<ROOT::RDF::RResultPtr<TH1D>, ROOT::RDF::RResultPtr<TH1D>>> dfhists;
	for(auto htbin : htbins){
		string htcutstring = makeHTBinCutstring(htbin);
		if(htcutstring == "")
			htcutstring = "(true)";
		auto df_htbin = df0.Filter(htcutstring);
		//total hist (denom)
		TH1D htotal_model("hTotal","PFMETOR_Efficiency;MET [GeV];Efficiency",NBINS,MET_MIN,MET_MAX);
		TH1D hpass_model("hPass","PFMETOR_Efficiency;MET [GeV];Efficiency",NBINS,MET_MIN,MET_MAX);
		string presel = trigpresel+" && "+evtfilters + " && "+htcutstring; //in addition to ntuple "filter" of (>= 1 SV || >= 1 photon[pt > 30]) && MET > 100
	
		presel += " && ((nBaseLinePhotons > 0) || ((SV_nLeptonic > 0) || (SV_nHadronic > 0)))";
	
		auto hTotal = df0.Filter(presel).Histo1D(htotal_model,"selCMet");	
		//pass hist (num)
		auto hPass = df0.Filter(presel+" && "+trigsel).Histo1D(hpass_model,"selCMet");

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
	map<string, string> filedirs = {};
	filedirs["18"] = "EGamma_Run2018C-12Nov2019_UL2018-v2__SVHPM100_v34__rjrskim_v51.root";
	filedirs["24"] = "EGamma0_Run2024C-MINIv6NANOv15-v1__SVHPM100_v34__rjrskim_v51.root";
	

	std::string path =
	    "root://cmseos.fnal.gov//store/group/lpcsusylep/malazaro/KUCMSSkims/skims_v51/";
	vector<string> preseltriggers = {"Trigger_Ele35_WPTight_Gsf","Trigger_Mu55","Trigger_Mu12","Trigger_IsoMu27","Trigger_IsoMu20","Trigger_Ele27_WPTight_Gsf"};
	vector<string> triggers = {"Trigger_PFMET120_PFMHT120_IDTight","Trigger_PFMET120_PFMHT120_IDTight_PFHT60","Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight","Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"};
	vector<string> htbins = {"_500","500_700","700_1000","1000_","_"};
	vector<TEfficiency*> effcurves;
	for(auto it = filedirs.begin(); it != filedirs.end(); it++){
		string year = it->first;
		cout << "Doing year " << year << endl;
		string skimfile = path+it->second;	
		makeTriggerEfficiency(skimfile,triggers,preseltriggers,htbins,year,effcurves);
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

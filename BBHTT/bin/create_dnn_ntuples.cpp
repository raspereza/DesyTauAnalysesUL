
#include <iostream>
#include <map>
#include "TString.h"
#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1.h"
#include "TList.h"
#include "TSystem.h"
#include "DesyTauAnalyses/BBHTT/interface/settings_for_eras.h"

bool applyPreselection = true;

//define the main function to make this code executable. argc gives the number of arguments and argv contains those arguments. 

int main(int argc, char * argv[]) {

  bool TEST = false;

  TString process(argv[1]);
  TString era(argv[2]);
  TString channel(argv[3]);

  TString Sample = process + "_" + era ;

if (argc!=4) {
    cout << "Usage of the program : CreateDNN sample era channel" << endl;
    exit(EXIT_FAILURE);
  }

  if (era!="2018"&&era!="2017"&&era!="2016") {
    cout << "ERROR : specified era " << era << " is unknown " << endl;
    exit(EXIT_FAILURE);
  }
  bool sampleFound = false;
  for (auto name_sample : map_sample) {
    if (name_sample.first==Sample) {
      sampleFound = true;
      break;
    }
  }
  if (!sampleFound) {
    cout << "ERROR : specified sample " << Sample << " is unknown " << endl;
    cout << "Available samples ===> " << endl;
    for (auto name_sample : map_sample) {
      cout << "     " << name_sample.first << endl;
    }
    exit(EXIT_FAILURE);
  }

//void create_dnn_ntuples( TString era = "2016" ){

  // Some definitions
  //TString channel = "em";
  double luminosity = 0;
  float qcd_ss_os_iso_relaxed_ratio = 0;
  float trigger_filter_efficiency = 1;
  float embedded_trigger_weight = 1.0;
  float embedded_tracking_weight = 1.0;
  TString input_dir;
  bool IsEra2016 = 0;
  bool IsEra2017 = 0;
  bool IsEra2018 = 0;

// Mapping of subsamples to output root-file
  map< TString , vector<TString> > samples_map;
  const map<TString, double>  *xsec_map    = 0;
  const map<TString, TString> *process_map = 0;
  
  //TString output_dir = "";
  samples_map[channel + "-NOMINAL_ntuple_"+Sample     ] = map_sample.at(Sample);
  if(channel=="em")
      input_dir = "/nfs/dust/cms/user/makou/Higgs_Tau_Tau/CMSSW_10_2_22/src/DesyTauAnalyses/NTupleMaker/test/emu_MSSM/2016";
    //input_dir= "/nfs/dust/cms/user/rasp/Run/emu_MSSM/Feb10/" + era;
//input_dir ="/nfs/dust/cms/user/rasp/storage/cardinia/SynchNTuples/mutau_June2-PAS/" + era ;
  else 
    input_dir ="/nfs/dust/cms/user/rasp/storage/cardinia/SynchNTuples/etau_Aug27/" + era ;
 
 

  // Mapping of subsamples to output root-file
  //map< TString , vector<TString> > samples_map;
  //const map<TString, double>  *xsec_map    = 0;
  //const map<TString, TString> *process_map = 0;
  //map<TString , int> n_events_per_sample;
  // TH1F *hNonClosureCorrection_njet1 = NULL;
  // TH1F *hNonClosureCorrection_njetgt1 = NULL;
  // TH1F *hNonClosureCorrection_njet0 = NULL;
  // TH2F *hIsolationCorrection= NULL;
  // TH2F *OS_SS_transfer_factors_deltar2to4= NULL;
  // TH2F *OS_SS_transfer_factors_deltarlt2= NULL;
  // TH2F *OS_SS_transfer_factors_deltarltp8= NULL;
  // TH2F *OS_SS_transfer_factors_deltargt4= NULL;

  //TH2F *hNonClosureCorrection = NULL;
  //TH2F *hIsolationCorrection= NULL;
  //TF1 *hOS_SS_transfer_factors_njetgt1= NULL;
  //TF1 *hOS_SS_transfer_factors_njet1= NULL;
  //TF1 *hOS_SS_transfer_factors_njet0= NULL;
  
  TString output_dir= "/nfs/dust/cms/user/makou/bbh_analysis_NTuples/HTauTau_emu/Inputs/NTuples_" + era; 
  //TString output_dir= "./NTuples_" + era; 
  gSystem -> Exec("mkdir " + output_dir);

if (era == "2018"){
     xsec_map    = &xsec_map_2018;
     process_map = &process_map_2018;
     luminosity  = 59740;
     trigger_filter_efficiency = 1.0;
     qcd_ss_os_iso_relaxed_ratio = 1.89; //number from Janek's talk in TauPOG meeting (10.04.19)
     embedded_trigger_weight  = 1.00;
     embedded_tracking_weight = 1.00; 
     /*IsEra2018 = 1;
     TFile *f= new TFile("/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/HTauTau_emu/QCDEstimation/output_2018_withuncertainties_v3/QCDweights.root","READ");
     hNonClosureCorrection=(TH2F*)f->Get("NonClosureCorrection");
     hIsolationCorrection=(TH2F*)f->Get("IsolationCorrection");
     hOS_SS_transfer_factors_njetgt1=(TF1*)f->Get("OS_SS_transfer_factors_njetgt1");
     hOS_SS_transfer_factors_njet1=(TF1*)f->Get("OS_SS_transfer_factors_njet1");
     hOS_SS_transfer_factors_njet0=(TF1*)f->Get("OS_SS_transfer_factors_njet0");
     // OS_SS_transfer_factors_deltarltp8=(TH2F*)f->Get("OS_SS_transfer_factors_deltarltp8");
     // OS_SS_transfer_factors_deltar2to4=(TH2F*)f->Get("OS_SS_transfer_factors_deltar2to4");
     // OS_SS_transfer_factors_deltargt4=(TH2F*)f->Get("OS_SS_transfer_factors_deltargt4");
     // OS_SS_transfer_factors_deltarlt2=(TH2F*)f->Get("OS_SS_transfer_factors_deltarlt2");
     // hNonClosureCorrection_njet0=(TH1F*)f->Get("NonClosureCorrection_njet0");
     // hNonClosureCorrection_njet1=(TH1F*)f->Get("NonClosureCorrection_njet1");
     // hNonClosureCorrection_njetgt1=(TH1F*)f->Get("NonClosureCorrection_njetgt1");

     samples_map[channel + "-NOMINAL_ntuple_MuonEG"   ] = MuonEG_Run2018;
     samples_map[channel + "-NOMINAL_ntuple_Embedded" ] = Embedded_2018;
     samples_map[channel + "-NOMINAL_ntuple_DYJets"   ] = DYJets_2018;
     samples_map[channel + "-NOMINAL_ntuple_WJets"    ] = WJets_2018;
     samples_map[channel + "-NOMINAL_ntuple_TTbar"    ] = TTbar_2018;
     samples_map[channel + "-NOMINAL_ntuple_SingleTop"] = SingleTop_2018;
     samples_map[channel + "-NOMINAL_ntuple_Diboson"  ] = Diboson_2018;
     samples_map[channel + "-NOMINAL_ntuple_ggH"      ] = GluGluHToTauTau_2018;
     samples_map[channel + "-NOMINAL_ntuple_VBFH"     ] = VBFHToTauTau_2018;
     samples_map[channel + "-NOMINAL_ntuple_ZH"       ] = ZHToTauTau_2018;
     samples_map[channel + "-NOMINAL_ntuple_WH"       ] = WHToTauTau_2018;
     samples_map[channel + "-NOMINAL_ntuple_ggHWW"    ] = ggHToWW_2018;
     samples_map[channel + "-NOMINAL_ntuple_VBFHWW"   ] = VBFHToWW_2018;
     samples_map[channel + "-NOMINAL_ntuple_ttH"      ] = ttH_2018;
     samples_map[channel + "-NOMINAL_ntuple_ZHWW"       ] = ZHToWW_2018;
     samples_map[channel + "-NOMINAL_ntuple_WHWW"       ] = WHToWW_2018;
     // samples_map[channel + "-NOMINAL_ntuple_ggh_1J_PTH0to120"      ] = GluGluHToTauTau_STXS1p1_Bin104to105_2018;
     // samples_map[channel + "-NOMINAL_ntuple_ggh_2J"                ] = GluGluHToTauTau_STXS1p1_Bin107to109_2018;
     // samples_map[channel + "-NOMINAL_ntuple_ggh_1J_PTH120to200"    ] = GluGluHToTauTau_STXS1p1_Bin106_2018;
     // samples_map[channel + "-NOMINAL_ntuple_vbftopo"               ] = GluGluHToTauTau_STXS1p1_Bin110to113_2018; //noch weiter trennen?
     // samples_map[channel + "-NOMINAL_ntuple_ggh_PTHGT200"          ] = GluGluHToTauTau_STXS1p1_Bin101_2018;

     // samples_map[channel + "-NOMINAL_ntuple_qqh_2J"                ] = VBFHToTauTau_STXS1p1_Bin203to205_2018;
     // samples_map[channel + "-NOMINAL_ntuple_qqh_PTHGT200"          ] = VBFHToTauTau_STXS1p1_Bin206_2018;
     // samples_map[channel + "-NOMINAL_ntuple_vbftopo"               ] = VBFHToTauTau_STXS1p1_Bin207to210_2018;  //noch weiter trennen?
     input_dir="/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/master/CMSSW_10_2_22/src/DesyTauAnalyses/NTupleMaker/test/HTauTau_EMu_2018_all_eras/";  */
  }
  else if(era == "2017"){
    xsec_map    = &xsec_map_2017;
    process_map = &process_map_2017;
    //n_events_per_sample = n_events_per_sample_2017;
    luminosity  = 41530;                 // Take number from LUMI twiki : https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#SummaryTable
    trigger_filter_efficiency = 1.0;
    qcd_ss_os_iso_relaxed_ratio = 2.38;
    embedded_trigger_weight  = 1.00;
    embedded_tracking_weight = 0.99;
    /*IsEra2017 = 1;
    TFile *f= new TFile("/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/HTauTau_emu/QCDEstimation/output_2017_withuncertainties_v3/QCDweights.root","READ");
    hNonClosureCorrection=(TH2F*)f->Get("NonClosureCorrection");
    hIsolationCorrection=(TH2F*)f->Get("IsolationCorrection");
    hOS_SS_transfer_factors_njetgt1=(TF1*)f->Get("OS_SS_transfer_factors_njetgt1");
    hOS_SS_transfer_factors_njet1=(TF1*)f->Get("OS_SS_transfer_factors_njet1");
    hOS_SS_transfer_factors_njet0=(TF1*)f->Get("OS_SS_transfer_factors_njet0");
    samples_map[channel + "-NOMINAL_ntuple_MuonEG"   ] = MuonEG_Run2017;
    samples_map[channel + "-NOMINAL_ntuple_Embedded" ] = Embedded_2017;
    samples_map[channel + "-NOMINAL_ntuple_DYJets"   ] = DYJets_2017;
    samples_map[channel + "-NOMINAL_ntuple_WJets"    ] = WJets_2017;
    samples_map[channel + "-NOMINAL_ntuple_TTbar"    ] = TTbar_2017;
    samples_map[channel + "-NOMINAL_ntuple_SingleTop"] = SingleTop_2017;
    samples_map[channel + "-NOMINAL_ntuple_Diboson"  ] = Diboson_2017;
    samples_map[channel + "-NOMINAL_ntuple_ggH"      ] = GluGluHToTauTau_2017;
    samples_map[channel + "-NOMINAL_ntuple_VBFH"     ] = VBFHToTauTau_2017;
    samples_map[channel + "-NOMINAL_ntuple_ZH"       ] = ZHToTauTau_2017;
    samples_map[channel + "-NOMINAL_ntuple_WH"       ] = WHToTauTau_2017;
    samples_map[channel + "-NOMINAL_ntuple_ggHWW"    ] = ggHToWW_2017;
    samples_map[channel + "-NOMINAL_ntuple_VBFHWW"   ] = VBFHToWW_2017;
    samples_map[channel + "-NOMINAL_ntuple_ttH"      ] = ttH_2017;
    samples_map[channel + "-NOMINAL_ntuple_ZHWW"       ] = ZHToWW_2017;
    samples_map[channel + "-NOMINAL_ntuple_WHWW"       ] = WHToWW_2017;
    // samples_map[channel + "-NOMINAL_ntuple_ggh_1J_PTH0to120"      ] = GluGluHToTauTau_STXS1p1_Bin104to105_2017;
    // samples_map[channel + "-NOMINAL_ntuple_ggh_2J"                ] = GluGluHToTauTau_STXS1p1_Bin107to109_2017;
    // samples_map[channel + "-NOMINAL_ntuple_ggh_1J_PTH120to200"    ] = GluGluHToTauTau_STXS1p1_Bin106_2017;
    // samples_map[channel + "-NOMINAL_ntuple_vbftopo"               ] = GluGluHToTauTau_STXS1p1_Bin110to113_2017; //noch weiter trennen?
    // samples_map[channel + "-NOMINAL_ntuple_ggh_PTHGT200"          ] = GluGluHToTauTau_STXS1p1_Bin101_2017;
    // samples_map[channel + "-NOMINAL_ntuple_qqh_2J"                ] = VBFHToTauTau_STXS1p1_Bin203to205_2017;
    // samples_map[channel + "-NOMINAL_ntuple_qqh_PTHGT200"          ] = VBFHToTauTau_STXS1p1_Bin206_2017;
    // samples_map[channel + "-NOMINAL_ntuple_vbftopo"               ] = VBFHToTauTau_STXS1p1_Bin207to210_2017;  //noch weiter trennen?
    //input_dir="/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/newMETv2/CMSSW_9_4_9/src/DesyTauAnalyses/NTupleMaker/test/HTauTau_EMu_2017_all_eras/";
    input_dir="/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/master/CMSSW_10_2_22/src/DesyTauAnalyses/NTupleMaker/test/HTauTau_EMu_2017_all_eras/";*/
  }
  else if(era == "2016"){
    xsec_map    = &xsec_map_2016;
    process_map = &process_map_2016;
    //n_events_per_sample = n_events_per_sample_2016;
    luminosity  = 35920;               // Take number from LUMI twiki : https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#SummaryTable
    trigger_filter_efficiency = 0.979;
    qcd_ss_os_iso_relaxed_ratio = 2.3;
    embedded_trigger_weight  = 1.0;
    embedded_tracking_weight = 1.0;
   /* IsEra2016 = 1;
    TFile *f= new TFile("/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/HTauTau_emu/QCDEstimation/output_2016_withuncertainties_v3/QCDweights.root","READ");
    hNonClosureCorrection=(TH2F*)f->Get("NonClosureCorrection");
    hIsolationCorrection=(TH2F*)f->Get("IsolationCorrection");
    hOS_SS_transfer_factors_njetgt1=(TF1*)f->Get("OS_SS_transfer_factors_njetgt1");
    hOS_SS_transfer_factors_njet1=(TF1*)f->Get("OS_SS_transfer_factors_njet1");
    hOS_SS_transfer_factors_njet0=(TF1*)f->Get("OS_SS_transfer_factors_njet0");

    samples_map[channel + "-NOMINAL_ntuple_MuonEG"   ] = MuonEG_Run2016;
    samples_map[channel + "-NOMINAL_ntuple_Embedded" ] = Embedded_2016;
    samples_map[channel + "-NOMINAL_ntuple_DYJets"   ] = DYJets_2016;
    samples_map[channel + "-NOMINAL_ntuple_WJets"    ] = WJets_2016;
    samples_map[channel + "-NOMINAL_ntuple_TTbar"    ] = TTbar_2016;
    samples_map[channel + "-NOMINAL_ntuple_SingleTop"] = SingleTop_2016;
    samples_map[channel + "-NOMINAL_ntuple_Diboson"  ] = Diboson_2016;

    samples_map[channel + "-NOMINAL_ntuple_ggH"      ] = GluGluHToTauTau_2016;
    //samples_map[channel + "-NOMINAL_ntuple_VBFH"     ] = VBFHToTauTau_2016;
    //samples_map[channel + "-NOMINAL_ntuple_ZH"       ] = ZHToTauTau_2016;
    //samples_map[channel + "-NOMINAL_ntuple_WH"       ] = WHToTauTau_2016;
    //samples_map[channel + "-NOMINAL_ntuple_ggHWW"    ] = ggHToWW_2016;
    //samples_map[channel + "-NOMINAL_ntuple_VBFHWW"   ] = VBFHToWW_2016;
    //samples_map[channel + "-NOMINAL_ntuple_ttH"      ] = ttH_2016;
    //samples_map[channel + "-NOMINAL_ntuple_ZHWW"       ] = ZHToWW_2016;
    //samples_map[channel + "-NOMINAL_ntuple_WHWW"       ] = WHToWW_2016;
    samples_map[channel + "-NOMINAL_ntuple_BBH"       ] = BBH_2016;
    input_dir="/nfs/dust/cms/user/rasp/Run/emu_MSSM/Feb10/2016/";*/
  }

  // Needed for stitching
  double xsecWIncl = 0;
  double xsecW1Jets = 0;
  double xsecW2Jets = 0;
  double xsecW3Jets = 0;
  double xsecW4Jets = 0;
  double xsecDYIncl     = 0;
  double xsecDY1Jets    = 0;
  double xsecDY2Jets    = 0;
  double xsecDY3Jets    = 0;
  double xsecDY4Jets    = 0;

  double neventsWIncl = 0;
  double neventsW1Jets  = 0;
  double neventsW2Jets  = 0;
  double neventsW3Jets  = 0;
  double neventsW4Jets  = 0;
  double neventsDYIncl  = 0;
  double neventsDY1Jets = 0;
  double neventsDY2Jets = 0;
  double neventsDY3Jets = 0;
  double neventsDY4Jets = 0;
  double neventsGluGluHIncl = 0;
  double neventsGluGluHBin104to105   = 0;
  double neventsGluGluHBin107to109   = 0; 
  double neventsGluGluHBin106        = 0; 
  double neventsGluGluHBin110to113   = 0; 
  double neventsGluGluHBin101        = 0; 
  double neventsVBFIncl              = 0; 
  double neventsVBFBin203to205    = 0; 
  double neventsVBFBin206         = 0; 
  double neventsVBFBin207to210    = 0; 

  if (Sample.Contains("WJets")){
     neventsWIncl   = getNEventsProcessed(input_dir+"/"+process_map->at("WJets")+".root");
     neventsW1Jets  = getNEventsProcessed(input_dir+"/"+process_map->at("W1Jets")+".root");
     neventsW2Jets  = getNEventsProcessed(input_dir+"/"+process_map->at("W2Jets")+".root");
     neventsW3Jets  = getNEventsProcessed(input_dir+"/"+process_map->at("W3Jets")+".root");
     neventsW4Jets  = getNEventsProcessed(input_dir+"/"+process_map->at("W4Jets")+".root");
     xsecWIncl      = xsec_map->at(process_map->at("WJets"));
     xsecW1Jets     = xsec_map->at(process_map->at("W1Jets"));
     xsecW2Jets     = xsec_map->at(process_map->at("W2Jets"));
     xsecW3Jets     = xsec_map->at(process_map->at("W3Jets"));
     xsecW4Jets     = xsec_map->at(process_map->at("W4Jets"));
  }
  if (Sample.Contains("DYJets")){
     neventsDYIncl  = getNEventsProcessed(input_dir+"/"+process_map->at("DYJets")+".root");
     neventsDY1Jets = getNEventsProcessed(input_dir+"/"+process_map->at("DY1Jets")+".root");
     neventsDY2Jets = getNEventsProcessed(input_dir+"/"+process_map->at("DY2Jets")+".root");
     neventsDY3Jets = getNEventsProcessed(input_dir+"/"+process_map->at("DY3Jets")+".root");
     neventsDY4Jets = getNEventsProcessed(input_dir+"/"+process_map->at("DY4Jets")+".root");
     xsecDYIncl     = xsec_map->at(process_map->at("DYJets"));
     xsecDY1Jets    = xsec_map->at(process_map->at("DY1Jets"));
     xsecDY2Jets    = xsec_map->at(process_map->at("DY2Jets"));
     xsecDY3Jets    = xsec_map->at(process_map->at("DY3Jets"));
     xsecDY4Jets    = xsec_map->at(process_map->at("DY4Jets"));
  }/*
  else{
     neventsWIncl   = n_events_per_sample.at(process_map->at("WJets"));
     neventsW1Jets  = n_events_per_sample.at(process_map->at("W1Jets"));
     neventsW2Jets  = n_events_per_sample.at(process_map->at("W2Jets"));
     neventsW3Jets  = n_events_per_sample.at(process_map->at("W3Jets"));
     neventsW4Jets  = n_events_per_sample.at(process_map->at("W4Jets"));
     neventsDYIncl  = n_events_per_sample.at(process_map->at("DYJets"));
     neventsDY1Jets = n_events_per_sample.at(process_map->at("DY1Jets"));
     neventsDY2Jets = n_events_per_sample.at(process_map->at("DY2Jets"));
     neventsDY3Jets = n_events_per_sample.at(process_map->at("DY3Jets"));
     neventsDY4Jets = n_events_per_sample.at(process_map->at("DY4Jets"));
  }
  if (era!="2016"){
     neventsGluGluHIncl        = getNEventsProcessed(input_dir+"/GluGluHToTauTau_M125_13TeV_powheg_pythia8.root");
     neventsGluGluHBin104to105 = getNEventsProcessed(input_dir+"/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin104to105_M125.root");
     neventsGluGluHBin107to109 = getNEventsProcessed(input_dir+"/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin107to109_M125.root");
     neventsGluGluHBin106      = getNEventsProcessed(input_dir+"/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin106_M125.root");     
     neventsGluGluHBin110to113 = getNEventsProcessed(input_dir+"/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin110to113_M125.root");
     neventsGluGluHBin101      = getNEventsProcessed(input_dir+"/GluGluHToTauTau_HTXSFilter_STXS1p1_Bin101_M125.root");     
     neventsVBFIncl            = getNEventsProcessed(input_dir+"/VBFHToTauTau_M125_13TeV_powheg_pythia8.root");           
     neventsVBFBin203to205  = getNEventsProcessed(input_dir+"/VBFHToTauTau_HTXSFilter_STXS1p1_Bin203to205_M125.root");
     neventsVBFBin206       = getNEventsProcessed(input_dir+"/VBFHToTauTau_HTXSFilter_STXS1p1_Bin206_M125.root");      
     neventsVBFBin207to210  = getNEventsProcessed(input_dir+"/VBFHToTauTau_HTXSFilter_STXS1p1_Bin207to210_M125.root");
  }*/
  //TString output_dir = "NTuples_" + era;
  //gSystem -> Exec("mkdir " + output_dir);

  // Loop over all samples

  cout <<"loop start " << endl;
  for (auto const& sample : samples_map){
    cout <<"loop start " << endl;
    cout << endl << sample.first << "  :  " << endl ;

    TFile *outFile = new TFile(output_dir + "/" + sample.first + ".root","RECREATE");
    TTree *outTree = new TTree("TauCheck", "tree created as DNN input");
    bool firstTree = true;
    TList* treeList = new TList();

    
    for(TString const& subsample: sample.second) {
       bool FoundSampleInMap=false;
      cout << "  - " << subsample << " : ";

      TFile *inFile  = new TFile( input_dir + "/" + subsample + ".root" ,"READ");
      TTree *inTree  = (TTree*) inFile -> Get("TauCheck");
      double nevents=0.;
      /*
	for (auto const& sample_evt : n_events_per_sample){
         if (subsample==sample_evt.first){
            nevents=n_events_per_sample.at(subsample);
            FoundSampleInMap=true;
         }
      }*/
      //if (!FoundSampleInMap){
      nevents = getNEventsProcessed( input_dir + "/" + subsample + ".root" );
      //}
      
      // SetBranchAddress for variables that need are needed for preselection or stitching
      int gen_noutgoing;
      float iso_1;
      float iso_2;
      int extraelec_veto;
      int extramuon_veto;
      float pt_1;
      float pt_2;
      float pzeta;
      //float mTdileptonMET;
      bool passedAllMetFilters;
      bool trg_muonelectron;
      bool trg_muhigh_elow;
      bool trg_ehigh_mulow;
      unsigned int run;
      int njets;
      float mjj;
      float jdeta;
      float pt_tt;
      /* float pt_ttjj;
	 float dijetpt;*/
      float jpt_1;
      float jpt_2;
      float jeta_1;
      float jeta_2;
      float weightEMu;
      float qcdweight;
      //int htxs_stage1p1cat;
      float dr_tt;
      int os;
      inTree->SetBranchAddress("gen_noutgoing",&gen_noutgoing);
      inTree->SetBranchAddress("iso_1",&iso_1);
      inTree->SetBranchAddress("iso_2",&iso_2);
      inTree->SetBranchAddress("extraelec_veto",&extraelec_veto);
      inTree->SetBranchAddress("extramuon_veto",&extramuon_veto);
      inTree->SetBranchAddress("pt_1",&pt_1);
      inTree->SetBranchAddress("pt_2",&pt_2);
      inTree->SetBranchAddress("pzeta",&pzeta);
      //inTree->SetBranchAddress("mTdileptonMET",&mTdileptonMET);
      inTree->SetBranchAddress("passedAllMetFilters",&passedAllMetFilters);
      inTree->SetBranchAddress("trg_muhigh_elow",&trg_muhigh_elow);
      inTree->SetBranchAddress("trg_ehigh_mulow", &trg_ehigh_mulow);
      inTree->SetBranchAddress("run",&run);
      inTree->SetBranchAddress("njets",&njets);
      inTree->SetBranchAddress("mjj",&mjj);
      inTree->SetBranchAddress("jdeta",&jdeta);
      inTree->SetBranchAddress("pt_tt",&pt_tt);
      //inTree->SetBranchAddress("pt_ttjj",&pt_ttjj);
      // inTree->SetBranchAddress("dijetpt",&dijetpt);
      inTree->SetBranchAddress("jpt_1",&jpt_1);
      inTree->SetBranchAddress("jpt_2",&jpt_2);
      inTree->SetBranchAddress("jeta_1",&jeta_1);
      inTree->SetBranchAddress("jeta_2",&jeta_2);
      // inTree->SetBranchAddress("htxs_stage1p1cat",&htxs_stage1p1cat);
      inTree->SetBranchAddress("dr_tt",&dr_tt);
      inTree->SetBranchAddress("os",&os);
      inTree->SetBranchAddress("weightEMu",&weightEMu);
      inTree->SetBranchAddress("qcdweight",&qcdweight);

      outFile->cd();
      
      TTree *currentTree = new TTree(subsample,"temporary tree");

      // Create a branch for xsec_lumi_weight
      float xsec_lumi_weight;
      float qcd_correction;
      float trigger_filter_weight;
      float embedded_stitching_weight;
      float embedded_rate_weight;
      int htxs_reco_flag_ggh;
      int htxs_reco_flag_qqh;
      bool era2016;
      bool era2017;
      bool era2018;
      float qcdweightml;
      float qcdweightml_1;
      float qcdweightml_2;
      float qcdweightml_3;
      
      if(firstTree){
	outTree    = inTree->CloneTree(0);
	outTree->Branch("xsec_lumi_weight", &xsec_lumi_weight, "xsec_lumi_weight/F");
	outTree->Branch("qcd_correction", &qcd_correction, "qcd_correction/F");
	outTree->Branch("trigger_filter_weight", &trigger_filter_weight, "trigger_filter_weight/F");
	outTree->Branch("embedded_stitching_weight", &embedded_stitching_weight, "embedded_stitching_weight/F");
	outTree->Branch("embedded_rate_weight", &embedded_rate_weight, "embedded_rate_weight/F");
	outTree->Branch("htxs_reco_flag_ggh", &htxs_reco_flag_ggh, "htxs_reco_flag_ggh/I");
	outTree->Branch("htxs_reco_flag_qqh", &htxs_reco_flag_qqh, "htxs_reco_flag_qqh/I");
  outTree->Branch("era2016", &era2016, "era2016/O");
  outTree->Branch("era2017", &era2017, "era2017/O");
  outTree->Branch("era2018", &era2018, "era2018/O");
  outTree->Branch("qcdweightml", &qcdweightml, "qcdweightml/F");
  outTree->Branch("qcdweightml_1", &qcdweightml_1, "qcdweightml_1/F");
  outTree->Branch("qcdweightml_2", &qcdweightml_2, "qcdweightml_2/F");
  outTree->Branch("qcdweightml_3", &qcdweightml_3, "qcdweightml_3/F");
	firstTree  = false;
      }
      currentTree = inTree->CloneTree(0);
      currentTree->Branch("xsec_lumi_weight", &xsec_lumi_weight, "xsec_lumi_weight/F");
      currentTree->Branch("qcd_correction", &qcd_correction, "qcd_correction/F");
      currentTree->Branch("trigger_filter_weight", &trigger_filter_weight, "trigger_filter_weight/F");
      currentTree->Branch("embedded_stitching_weight", &embedded_stitching_weight, "embedded_stitching_weight/F");
      currentTree->Branch("embedded_rate_weight", &embedded_rate_weight, "embedded_rate_weight/F");
      currentTree->Branch("htxs_reco_flag_ggh", &htxs_reco_flag_ggh, "htxs_reco_flag_ggh/I");
      currentTree->Branch("htxs_reco_flag_qqh", &htxs_reco_flag_qqh, "htxs_reco_flag_qqh/I");
      currentTree->Branch("era2016", &era2016, "era2016/O");
      currentTree->Branch("era2017", &era2017, "era2017/O");
      currentTree->Branch("era2018", &era2018, "era2018/O");
      currentTree->Branch("qcdweightml", &qcdweightml, "qcdweightml/F");
      currentTree->Branch("qcdweightml_1", &qcdweightml_1, "qcdweightml_1/F");
      currentTree->Branch("qcdweightml_2", &qcdweightml_2, "qcdweightml_2/F");
      currentTree->Branch("qcdweightml_3", &qcdweightml_3, "qcdweightml_3/F");
      // lumi-xsec-weight added
      if( xsec_map->find(subsample) == xsec_map->end() && !sample.first.Contains("MuonEG")  && !sample.first.Contains("Embedded")){
	cout << endl << endl << "Sample " << subsample << " is missing in xsec_map. Exit code." << endl << endl ;
	exit(-1);
      }
      float xsec = 1;
      if(!sample.first.Contains("MuonEG") && !sample.first.Contains("Embedded")) xsec = xsec_map->at(subsample);

      for (int i=0; i<inTree->GetEntries(); i++) {
         inTree->GetEntry(i);
         trg_muonelectron = (trg_muhigh_elow && pt_2 > 24) || (trg_ehigh_mulow && pt_1 > 24);
	// Add here preselection if necessary
	if(applyPreselection){
	  if( iso_1 > 0.5 )                continue;
	  if( iso_2 > 0.5 )                continue;
	  if( extraelec_veto > 0.5 )       continue;
	  if( extramuon_veto > 0.5 )       continue;
         
	  if( pt_1 < 10 )                  continue;
	  if( pt_2 < 10 )                  continue;
	  if( TMath::Max(pt_1,pt_2) < 20 ) continue;
	  //if( pzeta< -50 )                 continue;
	  //if( mTdileptonMET > 90 )         continue;
	  //if( passedAllMetFilters < 0.5 )           continue;
	  if( trg_muonelectron < 0.5 )     continue;
	}

	xsec_lumi_weight = xsec*luminosity/nevents;
	qcd_correction = qcd_ss_os_iso_relaxed_ratio;
	trigger_filter_weight = trigger_filter_efficiency;
        
  //cout <<"here"<<endl;
  era2016=IsEra2016;
  era2017=IsEra2017;
  era2018=IsEra2018;

  
  qcdweightml=1.0;
  qcdweightml_1=1.0;
  qcdweightml_2=1.0;
  qcdweightml_3=1.0;

  // if (dr_tt < 0.8){
  //     if (pt_1>pt_2){
  //        qcdweightml =OS_SS_transfer_factors_deltarltp8->GetBinContent(OS_SS_transfer_factors_deltarltp8->GetXaxis()->FindBin(pt_2),OS_SS_transfer_factors_deltarltp8->GetYaxis()->FindBin(pt_1));
  //        qcdweightml_1 =OS_SS_transfer_factors_deltarltp8->GetBinContent(OS_SS_transfer_factors_deltarltp8->GetXaxis()->FindBin(pt_2),OS_SS_transfer_factors_deltarltp8->GetYaxis()->FindBin(pt_1));
  //     }
  //     else {
  //         qcdweightml =OS_SS_transfer_factors_deltarltp8->GetBinContent(OS_SS_transfer_factors_deltarltp8->GetXaxis()->FindBin(pt_1),OS_SS_transfer_factors_deltarltp8->GetYaxis()->FindBin(pt_2));
  //        qcdweightml_1 =O
   //          S_SS_transfer_factors_deltarltp8->GetBinContent(OS_SS_transfer_factors_deltarltp8->GetXaxis()->FindBin(pt_1),OS_SS_transfer_factors_deltarltp8->GetYaxis()->FindBin(pt_2));
   //    }
   // }

   // if (dr_tt < 2 && dr_tt >= 0.8){
   //    if (pt_1>pt_2){
   //       qcdweightml =OS_SS_transfer_factors_deltarlt2->GetBinContent(OS_SS_transfer_factors_deltarlt2->GetXaxis()->FindBin(pt_2),OS_SS_transfer_factors_deltarlt2->GetYaxis()->FindBin(pt_1));
   //       qcdweightml_1 =OS_SS_transfer_factors_deltarlt2->GetBinContent(OS_SS_transfer_factors_deltarlt2->GetXaxis()->FindBin(pt_2),OS_SS_transfer_factors_deltarlt2->GetYaxis()->FindBin(pt_1));
   //    }
   //    else {
   //        qcdweightml =OS_SS_transfer_factors_deltarlt2->GetBinContent(OS_SS_transfer_factors_deltarlt2->GetXaxis()->FindBin(pt_1),OS_SS_transfer_factors_deltarlt2->GetYaxis()->FindBin(pt_2));
   //       qcdweightml_1 =OS_SS_transfer_factors_deltarlt2->GetBinContent(OS_SS_transfer_factors_deltarlt2->GetXaxis()->FindBin(pt_1),OS_SS_transfer_factors_deltarlt2->GetYaxis()->FindBin(pt_2));
   //    }
   // }
   // if (dr_tt >= 2 && dr_tt<4){
   //    if (pt_1>pt_2){
   //       qcdweightml =OS_SS_transfer_factors_deltar2to4->GetBinContent(OS_SS_transfer_factors_deltar2to4->GetXaxis()->FindBin(pt_2),OS_SS_transfer_factors_deltar2to4->GetYaxis()->FindBin(pt_1));
   //       qcdweightml_1 =OS_SS_transfer_factors_deltar2to4->GetBinContent(OS_SS_transfer_factors_deltar2to4->GetXaxis()->FindBin(pt_2),OS_SS_transfer_factors_deltar2to4->GetYaxis()->FindBin(pt_1));
   //    }
   //    else {
   //       qcdweightml =OS_SS_transfer_factors_deltar2to4->GetBinContent(OS_SS_transfer_factors_deltar2to4->GetXaxis()->FindBin(pt_1),OS_SS_transfer_factors_deltar2to4->GetYaxis()->FindBin(pt_2));
   //       qcdweightml_1 =OS_SS_transfer_factors_deltar2to4->GetBinContent(OS_SS_transfer_factors_deltar2to4->GetXaxis()->FindBin(pt_1),OS_SS_transfer_factors_deltar2to4->GetYaxis()->FindBin(pt_2));
   //    }
   // }
   // if (dr_tt >= 4 ){
   //    if (pt_1>pt_2){
   //       qcdweightml =OS_SS_transfer_factors_deltargt4->GetBinContent(OS_SS_transfer_factors_deltargt4->GetXaxis()->FindBin(pt_2),OS_SS_transfer_factors_deltargt4->GetYaxis()->FindBin(pt_1));
   //       qcdweightml_1 =OS_SS_transfer_factors_deltargt4->GetBinContent(OS_SS_transfer_factors_deltargt4->GetXaxis()->FindBin(pt_2),OS_SS_transfer_factors_deltargt4->GetYaxis()->FindBin(pt_1));
   //    }
   //    else {
   //      qcdweightml =OS_SS_transfer_factors_deltargt4->GetBinContent(OS_SS_transfer_factors_deltargt4->GetXaxis()->FindBin(pt_1),OS_SS_transfer_factors_deltargt4->GetYaxis()->FindBin(pt_2));
   //       qcdweightml_1 =OS_SS_transfer_factors_deltargt4->GetBinContent(OS_SS_transfer_factors_deltargt4->GetXaxis()->FindBin(pt_1),OS_SS_transfer_factors_deltargt4->GetYaxis()->FindBin(pt_2));
   //    }
   // }
  
   // if(njets==0){
   //    qcdweightml =hNonClosureCorrection_njet0->GetBinContent(hNonClosureCorrection_njet0->GetXaxis()->FindBin(1.0));
   //    qcdweightml_2 =hNonClosureCorrection_njet0->GetBinContent(hNonClosureCorrection_njet0->GetXaxis()->FindBin(1.0));
   // }
   // else if(njets ==1) {
   //    qcdweightml = hNonClosureCorrection_njet1->GetBinContent(hNonClosureCorrection_njet1->GetXaxis()->FindBin(1.0));
   //    qcdweightml_2 = hNonClosureCorrection_njet1->GetBinContent(hNonClosureCorrection_njet1->GetXaxis()->FindBin(1.0));
   // }
   // else {
   //   qcdweightml = hNonClosureCorrection_njetgt1->GetBinContent(hNonClosureCorrection_njetgt1->GetXaxis()->FindBin(1.0));
   //    qcdweightml_2 = hNonClosureCorrection_njetgt1->GetBinContent(hNonClosureCorrection_njetgt1->GetXaxis()->FindBin(1.0));
   // }
   // qcdweightml = qcdweightml*hIsolationCorrection->GetBinContent(hIsolationCorrection->GetXaxis()->FindBin(pt_2),hIsolationCorrection->GetYaxis()->FindBin(pt_1));
   // qcdweightml_3 = hIsolationCorrection->GetBinContent(hIsolationCorrection->GetXaxis()->FindBin(pt_2),hIsolationCorrection->GetYaxis()->FindBin(pt_1));
 

  /* if (njets==0){
     qcdweightml =hOS_SS_transfer_factors_njet0->Eval(dr_tt);
     qcdweightml_1 =hOS_SS_transfer_factors_njet0->Eval(dr_tt);
  }
  if (njets==1){
        qcdweightml =hOS_SS_transfer_factors_njet1->Eval(dr_tt);
        qcdweightml_1 =hOS_SS_transfer_factors_njet1->Eval(dr_tt);
  }
  if (njets>1){
        qcdweightml =hOS_SS_transfer_factors_njetgt1->Eval(dr_tt);
        qcdweightml_1 =hOS_SS_transfer_factors_njetgt1->Eval(dr_tt);
  }
  qcdweightml = qcdweightml*hNonClosureCorrection->GetBinContent(hNonClosureCorrection->GetXaxis()->FindBin(pt_2),hNonClosureCorrection->GetYaxis()->FindBin(pt_1));
  qcdweightml_2 = hNonClosureCorrection->GetBinContent(hNonClosureCorrection->GetXaxis()->FindBin(pt_2),hNonClosureCorrection->GetYaxis()->FindBin(pt_1));
  qcdweightml = qcdweightml*hIsolationCorrection->GetBinContent(hIsolationCorrection->GetXaxis()->FindBin(pt_2),hIsolationCorrection->GetYaxis()->FindBin(pt_1));
  qcdweightml_3 = hIsolationCorrection->GetBinContent(hIsolationCorrection->GetXaxis()->FindBin(pt_2),hIsolationCorrection->GetYaxis()->FindBin(pt_1));*/
 




  // Replace jet variables to have an effectie cut of jetpt > 30 GeV
	if(njets < 2){
	  // jdeta   = -10;
	  // m_jj     = -10;
	  // dijetpt = -10;
	  // pt_ttjj = -10;
	  jpt_2   = -10;
	  jeta_2  = -10;
	  if(njets < 1){
	    jpt_1 = -10;
	    jeta_1= -10;
	  }
	}
	
	// Stitching only for wjets MC in n-jet binned samples in npartons
	if( subsample.Contains("W") && subsample.Contains("JetsToLNu") ){
	  if(gen_noutgoing== 1)      xsec_lumi_weight = luminosity / ( neventsW1Jets/xsecW1Jets + neventsWIncl/xsecWIncl );
	  else if(gen_noutgoing== 2) xsec_lumi_weight = luminosity / ( neventsW2Jets/xsecW2Jets + neventsWIncl/xsecWIncl );
	  else if(gen_noutgoing== 3) xsec_lumi_weight = luminosity / ( neventsW3Jets/xsecW3Jets + neventsWIncl/xsecWIncl );
	  else if(gen_noutgoing== 4) xsec_lumi_weight = luminosity / ( neventsW4Jets/xsecW4Jets + neventsWIncl/xsecWIncl );
	  else                   xsec_lumi_weight = luminosity / ( neventsWIncl/xsecWIncl );
	}
	else if( subsample.Contains("DY") && subsample.Contains("JetsToLL_M-50") ){
	  if(gen_noutgoing== 1)      xsec_lumi_weight = luminosity / ( neventsDY1Jets/xsecDY1Jets + neventsDYIncl/xsecDYIncl );
	  else if(gen_noutgoing== 2) xsec_lumi_weight = luminosity / ( neventsDY2Jets/xsecDY2Jets + neventsDYIncl/xsecDYIncl );
	  else if(gen_noutgoing== 3) xsec_lumi_weight = luminosity / ( neventsDY3Jets/xsecDY3Jets + neventsDYIncl/xsecDYIncl );
	  else if(gen_noutgoing== 4) xsec_lumi_weight = luminosity / ( neventsDY4Jets/xsecDY4Jets + neventsDYIncl/xsecDYIncl );
	  else                   xsec_lumi_weight = luminosity / ( neventsDYIncl/xsecDYIncl );
	}
   // Stitching for STXS binned signal samples
	/* if( subsample.Contains("GluGluHToTauTau") && era!="2016"){
      if (htxs_stage1p1cat == 104 || htxs_stage1p1cat == 105) {
         xsec_lumi_weight = luminosity * xsec_map->at("GluGluHToTauTau_M125_13TeV_powheg_pythia8") / (neventsGluGluHIncl + neventsGluGluHBin104to105 / xsec_map->at("GluGluHToTauTau_HTXSFilter_STXS1p1_Bin104to105_M125"));
      }
      else if (htxs_stage1p1cat == 107 || htxs_stage1p1cat == 108 || htxs_stage1p1cat == 109) {
         xsec_lumi_weight = luminosity * xsec_map->at("GluGluHToTauTau_M125_13TeV_powheg_pythia8") / (neventsGluGluHIncl + neventsGluGluHBin107to109 / xsec_map->at("GluGluHToTauTau_HTXSFilter_STXS1p1_Bin107to109_M125"));
      }
      else if (htxs_stage1p1cat == 106) {
         xsec_lumi_weight = luminosity * xsec_map->at("GluGluHToTauTau_M125_13TeV_powheg_pythia8") / (neventsGluGluHIncl + neventsGluGluHBin106 / xsec_map->at("GluGluHToTauTau_HTXSFilter_STXS1p1_Bin106_M125"));
      }
      else if (htxs_stage1p1cat >=110 && htxs_stage1p1cat <=113) {
         xsec_lumi_weight = luminosity * xsec_map->at("GluGluHToTauTau_M125_13TeV_powheg_pythia8") / (neventsGluGluHIncl + neventsGluGluHBin110to113 / xsec_map->at("GluGluHToTauTau_HTXSFilter_STXS1p1_Bin110to113_M125"));
      }
      else if (htxs_stage1p1cat == 101) {
         xsec_lumi_weight = luminosity * xsec_map->at("GluGluHToTauTau_M125_13TeV_powheg_pythia8") / (neventsGluGluHIncl + neventsGluGluHBin101 / xsec_map->at("GluGluHToTauTau_HTXSFilter_STXS1p1_Bin101_M125"));
      }
      else {
         xsec_lumi_weight = luminosity / ( neventsGluGluHIncl / xsec_map->at("GluGluHToTauTau_M125_13TeV_powheg_pythia8"));
      }
   }

   if( subsample.Contains("VBFHToTauTau") && era!="2016"){
      if (htxs_stage1p1cat == 206) {
            xsec_lumi_weight = luminosity * xsec_map->at("VBFHToTauTau_M125_13TeV_powheg_pythia8") / (neventsVBFIncl + neventsVBFBin206 / xsec_map->at("VBFHToTauTau_HTXSFilter_STXS1p1_Bin206_M125"));
      }
      else if (htxs_stage1p1cat == 203 || htxs_stage1p1cat == 204 || htxs_stage1p1cat == 205) {
         xsec_lumi_weight = luminosity * xsec_map->at("VBFHToTauTau_M125_13TeV_powheg_pythia8") / (neventsVBFIncl + neventsVBFBin203to205 / xsec_map->at("VBFHToTauTau_HTXSFilter_STXS1p1_Bin203to205_M125"));
      }
      else if (htxs_stage1p1cat >=207 && htxs_stage1p1cat <=210) {
         xsec_lumi_weight = luminosity * xsec_map->at("VBFHToTauTau_M125_13TeV_powheg_pythia8") / (neventsVBFIncl + neventsVBFBin207to210 / xsec_map->at("VBFHToTauTau_HTXSFilter_STXS1p1_Bin207to210_M125"));
      }
      else {
         xsec_lumi_weight = luminosity / ( neventsVBFIncl / xsec_map->at("VBFHToTauTau_M125_13TeV_powheg_pythia8"));
      }
      }*/

	if( sample.first.Contains("MuonEG") || sample.first.Contains("Embedded")){

	  trigger_filter_weight = 1.;
	}
	xsec_lumi_weight = xsec_lumi_weight * weightEMu;
	embedded_stitching_weight = 1.;

	embedded_rate_weight = embedded_trigger_weight * embedded_tracking_weight;

	// add flags for cut categories which correspond to htxs_stage1p1cats
	htxs_reco_flag_ggh = 0;
	htxs_reco_flag_qqh = 0;
	// reco bins corresponding to stage 1.1 categories (see: https://indico.cern.ch/event/788748/contributions/3492755/attachments/1882178/3101573/SMHTT_ML_2017_status.pdf)
	/*	if(njets==0)                    htxs_reco_flag_ggh = 100;
	else if(njets==1 && pt_tt<120)  htxs_reco_flag_ggh = 101;
	else if(njets==1 && pt_tt>=120) htxs_reco_flag_ggh = 102;

	else if(njets>=2)               htxs_reco_flag_ggh = 103;

	if(njets<=1)                                htxs_reco_flag_qqh = 200;
	else if(njets>=2 && m_jj<350)                htxs_reco_flag_qqh = 201;
	else if(njets>=2 && m_jj>=350 && pt_tt<200)  htxs_reco_flag_qqh = 202;
	else if(njets>=2 && m_jj>=350 && pt_tt>=200) htxs_reco_flag_qqh = 203;

	// Select hadronic and leptonic part of VH sample
	if( subsample.Contains("VH") || subsample.Contains("WplusH") || subsample.Contains("WminusH") || subsample.Contains("HZJ") || subsample.Contains("HWminusJ") || subsample.Contains("HWplusJ")){
	  if( sample.first.Contains("VBFH") && (htxs_stage1p1cat>210 || htxs_stage1p1cat<200) ) continue;
	  if( (sample.first.Contains("WH") || sample.first.Contains("ZH")) && htxs_stage1p1cat<=210 && htxs_stage1p1cat>=200 ) continue;
	}
  // select hadronic and leptonic part of ggZH sample
  if(subsample.Contains("GluGluZH") || subsample.Contains("ggZH")){
	  if( (sample.first.Contains("GluGluH") || sample.first.Contains("ggH"))&& (htxs_stage1p1cat>116 || htxs_stage1p1cat< 100) ) continue;
	  if( (sample.first.Contains("ZH")) && htxs_stage1p1cat<=116 && htxs_stage1p1cat>=100 ) continue;
	  }*/

	currentTree->Fill();
      }
      cout<<currentTree->GetEntries()<<endl;
      treeList->Add(currentTree);
      inFile->Delete();
    }
    outTree = TTree::MergeTrees(treeList);
    cout<<"Together : "<<outTree->GetEntries()<<endl;
    outTree  -> Write( "" , TObject::kOverwrite );
    treeList -> Delete();
    outFile  -> Close();
  }
  cout << endl; 
}

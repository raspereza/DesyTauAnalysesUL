#ifndef SETTINGS_FOR_ERAS_H
#define SETTINGS_FOR_ERAS_H

using namespace std;
// **************************************************************************************************
// Define the subsamples that belong to a certain process
// 2018
const vector<TString> MuonEG_2018       = { "MuonEG_Run2018A" ,
                                               "MuonEG_Run2018B", 
                                               "MuonEG_Run2018C" ,
                                               "MuonEG_Run2018D" };
const vector<TString> Embedded_2018        = { "EmbeddedElMu_Run2018A" ,
                                               "EmbeddedElMu_Run2018B" ,
                                               "EmbeddedElMu_Run2018C" , 
                                               "EmbeddedElMu_Run2018D" };
const vector<TString> DYJets_2018          = { "DY1JetsToLL_M-50" ,
                                               "DY2JetsToLL_M-50" ,
                                               "DY3JetsToLL_M-50" ,
                                               "DY4JetsToLL_M-50" ,
                                               "DYJetsToLL_M-50",
                                               //"EWKZ2Jets_ZToNuNu_TuneCP5_PSweights_13TeV-madgraph-pythia8" ,
                                               //"EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8" , 
                                               //"DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8"
 };
const vector<TString> WJets_2018           = { "W1JetsToLNu" ,
                                               "W2JetsToLNu" ,
                                               "W3JetsToLNu" ,
                                               "W4JetsToLNu" ,
                                               "WJetsToLNu" ,
                                              };
const vector<TString> TTbar_2018           = { "TTTo2L2Nu" ,
                                               "TTToHadronic" ,
                                               "TTToSemiLeptonic" };
const vector<TString> SingleTop_2018       = { "ST_t-channel_antitop_4f" ,
                                               "ST_t-channel_top_4f" ,
                                               "ST_tW_antitop_5f" ,
                                               "ST_tW_top_5f" };
const vector<TString> Diboson_2018         = { "WZTo2L2Q" ,
                                               "WZTo3LNu" ,
                                               "ZZTo2L2Q" ,
                                               "ZZTo4L" ,
                                               "VVTo2L2Nu" };
const vector<TString> GluGluHToTauTau_2018 = {"SUSYGluGluToHToTauTau_powheg_M125" };
const vector<TString> BBH_2018 =  {"SUSYGluGluToBBHToTauTau_powheg_M125"};
const vector<TString> VBFHToTauTau_2018    = { "VBFHToTauTau_M125"};
const vector<TString> ZHToTauTau_2018      = { "ZHToTauTau_M125_13TeV_powheg_pythia8", 
                                              "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8",
                                              "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8",
                                              "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8" };
const vector<TString> WHToTauTau_2018      = { "WplusHToTauTau_M125_13TeV_powheg_pythia8" , "WminusHToTauTau_M125_13TeV_powheg_pythia8" };
const vector<TString> ggHToWW_2018         = { "GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8",
                                               "GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5_PSweights"};
const vector<TString> VBFHToWW_2018        = { "VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8",
                                               "HZJ_HToWW_M125_13TeV_powheg_jhugen714_pythia8_TuneCP5",
                                               "HWminusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5",
                                               "HWplusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5"};
const vector<TString> ZHToWW_2018      = { "HZJ_HToWW_M125_13TeV_powheg_jhugen714_pythia8_TuneCP5", "GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5_PSweights" };
const vector<TString> WHToWW_2018      = { "HWminusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5" , "HWplusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5" };
const vector<TString> ttH_2018             = { "ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin104to105_2018 = { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin104to105_M125" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin107to109_2018 = { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin107to109_M125" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin106_2018      = { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin106_M125" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin110to113_2018 = { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin110to113_M125" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin101_2018      = { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin101_M125" };
// const vector<TString> VBFHToTauTau_STXS1p1_Bin203to205_2018    = { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin203to205_M125" };
// const vector<TString> VBFHToTauTau_STXS1p1_Bin206_2018         = { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin206_M125" };
// const vector<TString> VBFHToTauTau_STXS1p1_Bin207to210_2018    = { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin207to210_M125" };

// 2017
const vector<TString> MuonEG_2017       = {"MuonEG_Run2017B",
  "MuonEG_Run2017C",
  "MuonEG_Run2017D",
  "MuonEG_Run2017E",
  "MuonEG_Run2017F" };
const vector<TString> Embedded_2017        = { "EmbeddedElMu_Run2017B",
  "EmbeddedElMu_Run2017C",
  "EmbeddedElMu_Run2017D",
  "EmbeddedElMu_Run2017E",
  "EmbeddedElMu_Run2017F", };
const vector<TString> DYJets_2017          = { "DY1JetsToLL_M-50",
					       "DY2JetsToLL_M-50" ,
					       "DY3JetsToLL_M-50" ,
                      "DY4JetsToLL_M-50",
                      "DYJetsToLL_M-50",};
const vector<TString> WJets_2017           = { "W1JetsToLNu" ,
					       "W2JetsToLNu" ,
					       "W3JetsToLNu" ,
					       "W4JetsToLNu" ,
					       "WJetsToLNu" ,
					       };
const vector<TString> TTbar_2017           = { "TTTo2L2Nu" ,
					       "TTToHadronic" ,
                                               "TTToSemiLeptonic"};
const vector<TString> SingleTop_2017       = { "ST_t-channel_antitop_4f" ,
					       "ST_t-channel_top_4f" ,
					       "ST_tW_antitop_5f" ,
					       "ST_tW_top_5f" };
const vector<TString> Diboson_2017         = {  "VVTo2L2Nu", "WZTo2L2Q", "WZTo3LNu", "ZZTo2L2Q","ZZTo4L"}; 

const vector<TString> GluGluHToTauTau_2017 = { "SUSYGluGluToHToTauTau_powheg_M125" };

const vector<TString> BBH_2017 =  {"SUSYGluGluToBBHToTauTau_powheg_M125"};

const vector<TString> VBFHToTauTau_2017    = { "VBFHToTauTau_M125"};


const vector<TString> ZHToTauTau_2017      = { "ZHToTauTau_M125_13TeV_powheg_pythia8",
                                               "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8",
                                               "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8",
                                               "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8"};
const vector<TString> WHToTauTau_2017      = { "WplusHToTauTau_M125_13TeV_powheg_pythia8" , "WminusHToTauTau_M125_13TeV_powheg_pythia8" };
const vector<TString> ggHToWW_2017         = { "GluGluHToWWTo2L2Nu", "GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5"};
const vector<TString> VBFHToWW_2017        = { "VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8", "HZJ_HToWW_M125_13TeV_powheg_jhugen714_pythia8_TuneCP5", "HWminusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5","HWplusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5"};
const vector<TString> ZHToWW_2017      = { "HZJ_HToWW_M125_13TeV_powheg_jhugen714_pythia8_TuneCP5", "GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5" };
const vector<TString> WHToWW_2017      = { "HWminusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5", "HWplusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5" };
const vector<TString> ttH_2017             = { "ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin104to105_2017 = { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin104to105_M125" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin107to109_2017 = { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin107to109_M125" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin106_2017 =      { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin106_M125" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin110to113_2017 = { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin110to113_M125" };
// const vector<TString> GluGluHToTauTau_STXS1p1_Bin101_2017 =      { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin101_M125" };
// const vector<TString> VBFHToTauTau_STXS1p1_Bin203to205_2017 =    { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin203to205_M125" };
// const vector<TString> VBFHToTauTau_STXS1p1_Bin206_2017 =         { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin206_M125" };
// const vector<TString> VBFHToTauTau_STXS1p1_Bin207to210_2017 =    { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin207to210_M125" };

// 2016
//mb implemented the changes from Alexeis settings.h
const vector<TString> MuonEG_2016       = { "MuonEG_Run2016B" , "MuonEG_Run2016C" , "MuonEG_Run2016D" , "MuonEG_Run2016E" , "MuonEG_Run2016F" , "MuonEG_Run2016G" , "MuonEG_Run2016H" };
const vector<TString> Embedded_2016        = { //"Embedding_Run2016"
"EmbeddedElMu_Run2016B",
  "EmbeddedElMu_Run2016C",
  "EmbeddedElMu_Run2016D",
  "EmbeddedElMu_Run2016E",
  "EmbeddedElMu_Run2016F",
  "EmbeddedElMu_Run2016G",
  "EmbeddedElMu_Run2016H" };
const vector<TString> DYJets_2016          = { "DY1JetsToLL_M-50" , "DY2JetsToLL_M-50" , "DY3JetsToLL_M-50" , "DY4JetsToLL_M-50" , "DYJetsToLL_M-50" , //"DYJetsToLL_M-10to50" ,
 //"EWKZ2Jets_ZToLL", 
//"EWKZ2Jets_ZToNuNu" 
};
const vector<TString> WJets_2016           = { "W1JetsToLNu" , "W2JetsToLNu" , "W3JetsToLNu" , "W4JetsToLNu" , "WJetsToLNu" //, "WGToLNuG" , "WGstarToLNuEE" , "WGstarToLNuMuMu" , "EWKWPlus2Jets" , "EWKWMinus2Jet"
 };
const vector<TString> TTbar_2016           = { "TTTo2L2Nu", "TTToHadronic", "TTToSemiLeptonic" };
const vector<TString> SingleTop_2016       = { "ST_t-channel_antitop_4f" , "ST_t-channel_top_4f" , "ST_tW_antitop_5f" , "ST_tW_top_5f" };
//const vector<TString> Diboson_2016         = { "VVTo2L2Nu" , "WZTo1L1Nu2Q" , "WZTo1L3Nu" , "WZTo2L2Q" , "ZZTo2L2Q" , "ZZTo4L" , "WWToLNuQQ" };
const vector<TString> Diboson_2016         = {  "VVTo2L2Nu", "WZTo2L2Q", "WZTo3LNu", "ZZTo2L2Q","ZZTo4L"};


/*const vector<TString> GluGluHToTauTau_2016 = { "GluGluHToTauTau_M125", "GluGluHToTauTau_M125_ext1", "GluGluHToTauTau_M125_ext2","ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8", "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8", "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8" };*/

const vector<TString> BBH_2016 =  {"SUSYGluGluToBBHToTauTau_powheg_M125"};
const vector<TString> GluGluHToTauTau_2016 = {"SUSYGluGluToHToTauTau_powheg_M125"}; //GluGluHToTauTau_M125"};
const vector<TString> VBFHToTauTau_2016    = { "VBFHToTauTau_M125" //"VBFHToTauTau_M125" , //"VBFHToTauTau_M125_ext2", "ZHToTauTau_M125" , "WplusHToTauTau_M125" , "WminusHToTauTau_M125"
};
const vector<TString> ZHToTauTau_2016      = { "ZHToTauTau_M125" , "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8", "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8", "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8"};
const vector<TString> WHToTauTau_2016      = { "WplusHToTauTau_M125" , "WminusHToTauTau_M125" };
const vector<TString> ggHToWW_2016         = { "GluGluHToWWTo2L2Nu_M125" , "GluGluZH_HToWW_M125_13TeV_powheg_pythia8"};
const vector<TString> VBFHToWW_2016        = { "VBFHToWWTo2L2Nu_M125", "HZJ_HToWW_M125_13TeV_powheg_pythia8", "HWminusJ_HToWW_M125_13TeV_powheg_pythia8","HWplusJ_HToWW_M125_13TeV_powheg_pythia8" };
const vector<TString> ZHToWW_2016      = { "HZJ_HToWW_M125_13TeV_powheg_pythia8" , "GluGluZH_HToWW_M125_13TeV_powheg_pythia8"};
const vector<TString> WHToWW_2016      = { "HWminusJ_HToWW_M125_13TeV_powheg_pythia8", "HWplusJ_HToWW_M125_13TeV_powheg_pythia8" };
const vector<TString> ttH_2016             = { "ttHJetToTT_M125" };

// **************************************************************************************************
// **************************************************************************************************
// Cross-section map
// 2018 (checked with AN-19-177)
const map<TString, double> xsec_map_2018 = {

 {"WJetsToLNu"  , 61526.7}, // NNLO (1)
  {"W1JetsToLNu" , 9370.5}, // NNLO (2)
  {"W2JetsToLNu" , 3170.9}, // NNLO (3)
  {"W3JetsToLNu" , 1132.5}, // NNLO (4)
  {"W4JetsToLNu" , 631.5 }, // NNLO (5)
  {"DYJetsToLL_M-50"       , 6077.22},  // NNLO (20)
  {"DY1JetsToLL_M-50"      , 977.1}, // NNLO (20a)
  {"DY2JetsToLL_M-50"      , 347.3}, // NNLO (20b)
  {"DY3JetsToLL_M-50"      , 126.1}, // NNLO (20c)
  {"DY4JetsToLL_M-50"      , 71.67}, // NNLO (20d)
  {"TTTo2L2Nu"        , 88.29},  // NNLO (21)
  {"TTToHadronic"     , 377.96}, // NNLO (22)
  {"TTToSemiLeptonic" , 365.35}, // NNLO (23)
  {"ST_t-channel_top_4f"      , 136.02}, // ? (24) -> could be improved
  {"ST_t-channel_antitop_4f"  , 80.95}, // ? (25) -> could be improved
  {"ST_tW_top_5f"             , 35.85}, // ? (26) -> could be improved
  {"ST_tW_antitop_5f"         , 35.85}, // ? (27) -> could be improved
  {"VVTo2L2Nu"                , 13.84},
  {"WWToLNuQQ"                , 49.997},
  {"WZTo2L2Q"                 , 5.52},
  {"WZTo1L1Nu2Q"              , 10.71},
  {"WZTo1L3Nu"                , 3.05},
  {"WZJToLLLNu"               , 4.708},
  {"WZTo3LNu"                 , 4.43},
  {"ZZTo4L"                   , 1.26},
  {"ZZTo2L2Q"                 , 3.38},
  {"SUSYGluGluToHToTauTau_powheg_M125"     , 3.00},
  {"VBFHToTauTau_M125"        , 0.237},
  {"WplusHToTauTau_M125"      , 0.0527},
  {"WminusHToTauTau_M125"     , 0.0334},
  {"ZHToTauTau_M125_13TeV"    , 0.0477},
  {"GluGluHToWWTo2L2Nu_M125"  , 1.09},
  {"VBFHToWWTo2L2Nu_M125"     , 0.0850},
  {"HWminusJ_HToWW_M125"      , 0.114},
  {"HWplusJ_HToWW_M125"       , 0.18},
  {"ZHJ_HToWW_M125"           , 0.163},

  {"SUSYGluGluToBBHToTauTau_powheg_M125",  0.0307},

};


 /*  { "DYJetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8"  , 6077.22 }, // NNLO xsec -> k-factor=1.1137 (see mail from Izaak on tauid hn 14.10.19)
   { "DY1JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" , 877.8*1.165 },
   { "DY2JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" , 304.4*1.165 },
   { "DY3JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" , 111.5*1.165 },
   { "DY4JetsToLL_M-50_TuneCP5_13TeV_madgraphMLM_pythia8" , 44.03*1.165 },
   { "EWKZ2Jets_ZToNuNu_TuneCP5_PSweights_13TeV-madgraph-pythia8" , 10.66 },
   { "EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8" , 4.321 },
   { "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8" , 21658.0 },
   { "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"  , 61526.7 },   // NNLO -> k-factor=1.162
   { "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 8104.0*1.079 },
   { "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 2793.0*1.079 },
   { "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 992.5*1.079 },
   { "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 544.3*1.079 },
   { "EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8" , 23.24},
   { "EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8" , 29.59},
   { "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8" , 464.4},
   { "TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8"        , 88.29 },
   { "TTToHadronic_TuneCP5_13TeV_powheg_pythia8"     , 377.96 },
   { "TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8" , 365.35 },
   { "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8" , 80.95 },
   { "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"     , 136.02 },
   { "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"                  , 35.85 },
   { "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"                      , 35.85 },
   { "WW_TuneCP5_13TeV-pythia8" , 118.7 },  // NNLO
   { "WZ_TuneCP5_13TeV-pythia8" , 27.57 },
   { "ZZ_TuneCP5_13TeV-pythia8" , 12.14 },
   { "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8" , 5.52},
   { "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8" , 4.43},
   { "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8" , 3.38},
   { "ZZTo4L_TuneCP5_13TeV-amcatnloFXFX-pythia8" , 1.26},
   { "VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8" , 13.84},
   { "GluGluHToTauTau_M125_13TeV_powheg_pythia8" , 3.08 }, //bottom quark associated production with 1% of the gluon fusion cross section is included in ggh
   { "VBFHToTauTau_M125_13TeV_powheg_pythia8" , 3.782*0.0627 },
   { "ZHToTauTau_M125_13TeV_powheg_pythia8" ,  0.04774},
   { "WplusHToTauTau_M125_13TeV_powheg_pythia8" , 0.0527 },
   { "WminusHToTauTau_M125_13TeV_powheg_pythia8" , 0.0334 },
   { "GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8" , 1.091343 }, // 48.58 (Production) * 0.214 (H->WW) * 0.104976 (WW->2l2Nu), from Dennis Roy (Heavy H->WW analysis), SM H->WW analysis still uses old ggh cross section (according to their AN), other values are the same
   { "VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8" , 0.084962 },// 3.782 * 0.214 * 0.104976, from Dennis Roy (Heavy H->WW analysis)
   { "ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8" , 0.5071*0.0627 },
   { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin104to105_M125",  0.2211}, //for STXS binned samples this is the fraction = xsec_subsample/xsec_incl (LO)
   { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin107to109_M125",  0.0738},
   { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin106_M125",  0.0140},
   { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin110to113_M125", 0.0156 },
   { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin101_M125",  0.0139},
   { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin203to205_M125", 0.1608},
   { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin206_M125", 0.0398 },
   { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin207to210_M125", 0.3338},
   { "HZJ_HToWW_M125_13TeV_powheg_jhugen714_pythia8_TuneCP5" , 0.163},
   { "HWminusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5", 0.114}, 
   { "HWplusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5", 0.18},
   { "GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5_PSweights", 0.0262},
   { "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8",  0.00538017},
   { "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8", 0.00153915},
   { "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8", 0.00077719}
};*/
// 2017 (checked again - AN-19-177 -> DY xsec updated )
const map<TString, double> xsec_map_2017 = {
  {"WJetsToLNu"  , 61526.7}, // NNLO (1)
  {"W1JetsToLNu" , 9370.5}, // NNLO (2)
  {"W2JetsToLNu" , 3170.9}, // NNLO (3)
  {"W3JetsToLNu" , 1132.5}, // NNLO (4)
  {"W4JetsToLNu" , 631.5 }, // NNLO (5)
  {"DYJetsToLL_M-50"       , 6077.22},  // NNLO (20)
  {"DY1JetsToLL_M-50"      , 977.1}, // NNLO (20a)
  {"DY2JetsToLL_M-50"      , 347.3}, // NNLO (20b)
  {"DY3JetsToLL_M-50"      , 126.1}, // NNLO (20c)
  {"DY4JetsToLL_M-50"      , 71.67}, // NNLO (20d)
  {"TTTo2L2Nu"        , 88.29},  // NNLO (21)
  {"TTToHadronic"     , 377.96}, // NNLO (22)
  {"TTToSemiLeptonic" , 365.35}, // NNLO (23)
  {"ST_t-channel_top_4f"      , 136.02}, // ? (24) -> could be improved
  {"ST_t-channel_antitop_4f"  , 80.95}, // ? (25) -> could be improved
  {"ST_tW_top_5f"             , 35.85}, // ? (26) -> could be improved
  {"ST_tW_antitop_5f"         , 35.85}, // ? (27) -> could be improved
  {"VVTo2L2Nu"                , 13.84},
  {"WWToLNuQQ"                , 49.997},
  {"WZTo2L2Q"                 , 5.52},
  {"WZTo1L1Nu2Q"              , 10.71},
  {"WZTo1L3Nu"                , 3.05},
  {"WZJToLLLNu"               , 4.708},
  {"WZTo3LNu"                 , 4.43},
  {"ZZTo4L"                   , 1.26},
  {"ZZTo2L2Q"                 , 3.38},
  {"SUSYGluGluToHToTauTau_powheg_M125"     , 3.00},
  {"VBFHToTauTau_M125"        , 0.237},
  {"WplusHToTauTau_M125"      , 0.0527},
  {"WminusHToTauTau_M125"     , 0.0334},
  {"ZHToTauTau_M125_13TeV"    , 0.0477},
  {"GluGluHToWWTo2L2Nu_M125"  , 1.09},
  {"VBFHToWWTo2L2Nu_M125"     , 0.0850},
  {"HWminusJ_HToWW_M125"      , 0.114},
  {"HWplusJ_HToWW_M125"       , 0.18},
  {"ZHJ_HToWW_M125"           , 0.163},

  {"SUSYGluGluToBBHToTauTau_powheg_M125",  0.0307},

};

























  /*{ "DYJetsToLL_M-10to50_13TeV-12Apr2018"      , 21658.0 },
  { "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"  , 6077.22 }, // NNLO xsec -> k-factor=1.137 (see mail from Izaak on tauid hn 14.10.19)
  { "DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8" , 877.8*1.165 },
  { "DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8" , 304.4*1.165 },
  { "DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8" , 111.5*1.165 },
  { "DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8" , 44.04*1.165 },
  { "EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8" , 4.321 },
  { "EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8" , 10.66 },
  { "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"  , 61526.7 },   // NNLO -> k-factor=1.162
  { "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 8104.0*1.079 },
  { "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 2793.0*1.079 },
  { "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 992.5*1.079 },
  { "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 544.3*1.079 },
  { "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8"    , 464.4 },
  { "EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8" , 23.24 },
  { "EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8" , 29.59 },
  { "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8"        , 88.29 },
  { "TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8"     , 377.96 },
  { "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8" , 365.35 },
  { "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8" , 80.95 },
  { "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8"     , 136.02 },
  { "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"                  , 35.85 },
  { "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"                      , 35.85 },
  { "WW_TuneCP5_13TeV-pythia8" , 118.7 },  // NNLO
  { "WZ_TuneCP5_13TeV-pythia8" , 27.57 },
  { "ZZ_TuneCP5_13TeV-pythia8" , 12.14 },
  { "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8" , 5.52},
  { "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8" , 4.43},
  { "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8" , 3.38} ,
  { "ZZTo4L_TuneCP5_13TeV-amcatnloFXFX-pythia8" , 1.26},
  { "VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8", 13.84},
  { "GluGluHToTauTau_M125_13TeV_powheg_pythia8" , 3.08 }, //bottom quark associated production with 1% of the gluon fusion cross section is included in ggh
  { "VBFHToTauTau_M125_13TeV_powheg_pythia8"    , 3.782*0.0627 },
  { "ZHToTauTau_M125_13TeV_powheg_pythia8"      ,  0.04774 },
  { "WplusHToTauTau_M125_13TeV_powheg_pythia8"  , 0.0527 },
  { "WminusHToTauTau_M125_13TeV_powheg_pythia8" , 0.0334 },
  { "GluGluHToWWTo2L2Nu"                        , 1.091343 },  // see 2018
  { "VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8" , 0.084962 }, // see 2018
  { "ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8"         , 0.5071*0.0627 },  // from: see links above
  { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin104to105_M125", 0.2216 }, //for STXS binned samples this is the fraction = xsec_subsample/xsec_incl (LO)
  { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin107to109_M125", 0.0738 },
  { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin106_M125", 0.0141 },
  { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin110to113_M125", 0.0157},
  { "GluGluHToTauTau_HTXSFilter_STXS1p1_Bin101_M125", 0.0139 },
  { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin203to205_M125", 0.1606},
  { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin206_M125", 0.0398 },
  { "VBFHToTauTau_HTXSFilter_STXS1p1_Bin207to210_M125", 0.3341},
  { "HZJ_HToWW_M125_13TeV_powheg_jhugen714_pythia8_TuneCP5" , 0.163},
  { "HWminusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5", 0.114}, 
  { "HWplusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5", 0.18},
  { "GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5", 0.0262},
  { "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8",  0.00538017},
  { "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8", 0.00153915},
  { "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8", 0.00077719}*/


// 2016 (taken from AN-19-177)
const map<TString, double> xsec_map_2016 = {
 
  {"WJetsToLNu"  , 61526.7}, // NNLO (1)
  {"W1JetsToLNu" , 9370.5}, // NNLO (2)
  {"W2JetsToLNu" , 3170.9}, // NNLO (3)
  {"W3JetsToLNu" , 1132.5}, // NNLO (4)
  {"W4JetsToLNu" , 631.5 }, // NNLO (5)
  {"DYJetsToLL_M-50"       , 6077.22},  // NNLO (20)
  {"DY1JetsToLL_M-50"      , 977.1}, // NNLO (20a)
  {"DY2JetsToLL_M-50"      , 347.3}, // NNLO (20b)
  {"DY3JetsToLL_M-50"      , 126.1}, // NNLO (20c)
  {"DY4JetsToLL_M-50"      , 71.67}, // NNLO (20d)
  {"TTTo2L2Nu"        , 88.29},  // NNLO (21)
  {"TTToHadronic"     , 377.96}, // NNLO (22)
  {"TTToSemiLeptonic" , 365.35}, // NNLO (23)
  {"ST_t-channel_top_4f"      , 136.02}, // ? (24) -> could be improved
  {"ST_t-channel_antitop_4f"  , 80.95}, // ? (25) -> could be improved
  {"ST_tW_top_5f"             , 35.85}, // ? (26) -> could be improved
  {"ST_tW_antitop_5f"         , 35.85}, // ? (27) -> could be improved
  {"VVTo2L2Nu"                , 13.84},
  {"WWToLNuQQ"                , 49.997},
  {"WZTo2L2Q"                 , 5.52},
  {"WZTo1L1Nu2Q"              , 10.71},
  {"WZTo1L3Nu"                , 3.05},
  {"WZJToLLLNu"               , 4.708},
  {"WZTo3LNu"                 , 4.43},
  {"ZZTo4L"                   , 1.26},
  {"ZZTo2L2Q"                 , 3.38},
  {"SUSYGluGluToHToTauTau_powheg_M125"     , 3.00},
  {"VBFHToTauTau_M125"        , 0.237},
  {"WplusHToTauTau_M125"      , 0.0527},
  {"WminusHToTauTau_M125"     , 0.0334},
  {"ZHToTauTau_M125_13TeV"    , 0.0477},
  {"GluGluHToWWTo2L2Nu_M125"  , 1.09},
  {"VBFHToWWTo2L2Nu_M125"     , 0.0850},
  {"HWminusJ_HToWW_M125"      , 0.114},
  {"HWplusJ_HToWW_M125"       , 0.18},
  {"ZHJ_HToWW_M125"           , 0.163},
  //{"SUSYGluGluToBBHToTauTau_powheg_M125",  0.0307}

  {"SUSYGluGluToBBHToTauTau_powheg_M125",  0.0648}








  /*  { "DYJetsToLL_M-10to50"      , 21658.0 },
    { "DYJetsToLL_M-50"          , 6077.22 }, // NNLO xsec -> k-factor=1.227 (see mail from Izaak on tauid hn 14.10.19)
    { "DY1JetsToLL_M-50"         , 1.227*1012.5 },// to do: check!
    { "DY2JetsToLL_M-50"         , 1.227*332.8 },// to do: check!
    { "DY3JetsToLL_M-50"         , 1.227*101.8 },// to do: check!
    { "DY4JetsToLL_M-50"         , 1.227*54.8 },// to do: check!
    { "WJetsToLNu"               , 61526.7 },
    { "W1JetsToLNu"              , 1.221*9644.5 },// to do: check!
    { "W2JetsToLNu"              , 1.221*3144.5 },// to do: check!
    { "W3JetsToLNu"              , 1.221*954.8 },// to do: check!
    { "W4JetsToLNu"              , 1.221*485.6 },// to do: check!
    { "TTbar"                    , 831.6 },
    { "TTTo2L2Nu"                , 88.29}, 
    { "TTToHadronic"             , 377.96},
    { "TTToSemiLeptonic"         , 365.35},
    { "ST_t-channel_antitop_4f"     , 80.95 },
    { "ST_t-channel_top_4f"         , 136.02 },
    { "ST_tW_antitop_5f"            , 35.85 },
    { "ST_tW_top_5f"                , 35.85 },
    // { "VVTo2L2Nu"                , 11.95 },
    // { "WWToLNuQQ"                , 49.997 },
    // { "WZTo2L2Q"                 , 5.595 },
    // { "WZTo1L1Nu2Q"              , 10.71 },
    // { "WZTo1L3Nu"                , 3.05 },
    // { "WZJToLLLNu"               , 5.26 },
    // { "ZZTo4L"                   , 1.212 },
    // { "ZZTo2L2Q"                 , 3.22 },
    {"WW", 118.7},
    {"WZ", 27.57},
    {"ZZ", 12.14},
    {"VVTo2L2Nu", 13.84}, 
    {"WZTo2L2Q", 5.52},
    {"WZTo3LNu", 4.43},
    {"ZZTo2L2Q", 3.38},
    {"ZZTo4L", 1.26},
    { "WGToLNuG"                 , 178.4 }, // xsdb
    { "WGstarToLNuMuMu"          , 2.793 },
    { "WGstarToLNuEE"            , 3.526 },
    { "EWKWPlus2Jets"            , 29.59 },
    { "EWKWMinus2Jet"            , 23.24 },
    { "EWKZ2Jets_ZToLL"          , 4.321 },
    { "EWKZ2Jets_ZToNuNu"        , 10.66 },
    { "GluGluHToTauTau_M125"     , 3.08 },//bottom quark associated production with 1% of the gluon fusion cross section is included in ggh
    { "VBFHToTauTau_M125"        , 3.782*0.0627 },
    { "GluGluHToTauTau_M125_ext1" , 3.08 },//bottom quark associated production with 1% of the gluon fusion cross section is included in ggh
    //{ "VBFHToTauTau_M125_ext1"   , 3.782*0.0627 },
    { "GluGluHToTauTau_M125_ext2" , 3.08 },//bottom quark associated production with 1% of the gluon fusion cross section is included in ggh
    //{ "VBFHToTauTau_M125_ext2"   , 3.782*0.0627 }, 
    {"VBFHToTauTau_M125", 0.237     }, //MB added from Alexeis code
    { "ZHToTauTau_M125"          ,  0.04774 },
    { "WplusHToTauTau_M125"      , 0.0527 },
    { "WminusHToTauTau_M125"     , 0.0334 },
    { "GluGluHToWWTo2L2Nu_M125"  , 1.091343 }, // see 2018
    { "VBFHToWWTo2L2Nu_M125"     , 0.084962 }, // see 2018
    { "HZJ_HToWW_M125_13TeV_powheg_pythia8" , 0.163},
    { "HWminusJ_HToWW_M125_13TeV_powheg_pythia8" , 0.114},em-NOMINAL_ntuple_GluGluHToTauTau_2016.root

    { "HWplusJ_HToWW_M125_13TeV_powheg_pythia8" , 0.18},
    { "GluGluZH_HToWW_M125_13TeV_powheg_pythia8" , 0.0262},
    { "ttHJetToTT_M125"          , 0.5071*0.0627 },// see 2017
    { "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8",  0.00538017},
    { "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8", 0.00153915},
    { "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8", 0.00077719},*/
     //MB: here i got the error that is missing couldnt find it in Alexeis code
  };
// **************************************************************************************************
// **************************************************************************************************
// maps short name to filenames - needed for stitching


const map<TString, vector<TString> > map_sample = {
  {"MuonEG_2018", MuonEG_2018},
  {"MuonEG_2017", MuonEG_2017},
  {"MuonEG_2016", MuonEG_2016},
  

  {"Embedded_2018", Embedded_2018},
  {"Embedded_2017", Embedded_2017},
  {"Embedded_2016", Embedded_2016},

  {"DYJets_2018",DYJets_2018},
  {"DYJets_2017",DYJets_2017},
  {"DYJets_2016",DYJets_2016},

  {"WJets_2018",WJets_2018},
  {"WJets_2017",WJets_2017},
  {"WJets_2016",WJets_2016},

  {"TTbar_2018",TTbar_2018},
  {"TTbar_2017",TTbar_2017},
  {"TTbar_2016",TTbar_2016},

  {"SingleTop_2018",SingleTop_2018},
  {"SingleTop_2017",SingleTop_2017},
  {"SingleTop_2016",SingleTop_2016},

  {"Diboson_2018",Diboson_2018},
  {"Diboson_2017",Diboson_2017},
  {"Diboson_2016",Diboson_2016},

  {"GluGluHToTauTau_2018",GluGluHToTauTau_2018},
  {"GluGluHToTauTau_2017",GluGluHToTauTau_2017},
  {"GluGluHToTauTau_2016",GluGluHToTauTau_2016},

  {"BBH_2018",BBH_2018},
  {"BBH_2017",BBH_2017},
  {"BBH_2016",BBH_2016}
};
const map<TString , TString> process_map_2018 = {
  { "WJets"   , "WJetsToLNu"},
  { "W1Jets"  , "W1JetsToLNu"},
  { "W2Jets"  , "W2JetsToLNu"},
  { "W3Jets"  , "W3JetsToLNu"},
  { "W4Jets"  , "W4JetsToLNu"},
  { "DYJets"  , "DYJetsToLL_M-50"},
  { "DY1Jets" , "DY1JetsToLL_M-50"},
  { "DY2Jets" , "DY2JetsToLL_M-50"},
  { "DY3Jets" , "DY3JetsToLL_M-50"},
  { "DY4Jets" , "DY4JetsToLL_M-50"},
};

const map<TString , TString> process_map_2017 = {
  { "WJets"   , "WJetsToLNu"},
  { "W1Jets"  , "W1JetsToLNu"},
  { "W2Jets"  , "W2JetsToLNu"},
  { "W3Jets"  , "W3JetsToLNu"},
  { "W4Jets"  , "W4JetsToLNu"},
  { "DYJets"  , "DYJetsToLL_M-50"},
  { "DY1Jets" , "DY1JetsToLL_M-50"},
  { "DY2Jets" , "DY2JetsToLL_M-50"},
  { "DY3Jets" , "DY3JetsToLL_M-50"},
  { "DY4Jets" , "DY4JetsToLL_M-50"},
};

const map<TString , TString> process_map_2016 = {
  { "WJets"   , "WJetsToLNu"},
  { "W1Jets"  , "W1JetsToLNu"},
  { "W2Jets"  , "W2JetsToLNu"},
  { "W3Jets"  , "W3JetsToLNu"},
  { "W4Jets"  , "W4JetsToLNu"},
  { "DYJets"  , "DYJetsToLL_M-50"},
  { "DY1Jets" , "DY1JetsToLL_M-50"},
  { "DY2Jets" , "DY2JetsToLL_M-50"},
  { "DY3Jets" , "DY3JetsToLL_M-50"},
  { "DY4Jets" , "DY4JetsToLL_M-50"},
};

// **************************************************************************************************
// **************************************************************************************************
// maps for number of events in sample, needed for proper normalization in 2016 & 2017
map<TString, int> n_events_per_sample_2017 = {
{"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"  , 74635450},
{"W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 54988117},
{"W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 32368249},
{"W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 19700377},
{"W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" , 11333705},
{"ZZ_TuneCP5_13TeV-pythia8" , 1949768},
{"WW_TuneCP5_13TeV-pythia8" , 7791498},
{"WZ_TuneCP5_13TeV-pythia8" , 3928630},
{"DYJetsToLL_M-10to50_13TeV-12Apr2018" , 39521230},
{"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"     , 97800939},
{"DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"    , 34859434},
{"DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"    , 9790490},
{"DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"    , 6897933},
{"DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"    , 4346952},
{"EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8" , 3000000},
{"EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8" , 954336},
{"EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8" , 4449714},
{"EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8" , 4453013},
{"WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8" , 6283083}
};
//TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8, TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8, TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8 , ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8, ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8, ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8, ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8, GluGluHToTauTau_M125_13TeV_powheg_pythia8, VBFHToTauTau_M125_13TeV_powheg_pythia8, exclusive diboson samples can be normalized using GetSumOfWeights()

map<TString, int> n_events_per_sample_2016 = {
{"WJetsToLNu"  , 86916455},
{"W1JetsToLNu" , 45283121},
{"W2JetsToLNu" , 30374504},
{"W3JetsToLNu" , 39501912},
{"W4JetsToLNu" , 20824737},
{"ZZ" , 1988098},
{"WW" , 7982180},
{"WZ" , 3997571},
{"DYJetsToLL_M-10to50" , 35114961},
{"DYJetsToLL_M-50"     , 146280395},
{"DY1JetsToLL_M-50"    , 63730337},
{"DY2JetsToLL_M-50"    , 19879279},
{"DY3JetsToLL_M-50"    , 5857441},
{"DY4JetsToLL_M-50"    , 4197868},
{"TTbar" , 76915549},
{"ST_t-channel_top"     , 67105876},
{"ST_t-channel_antitop" , 38811017},
{"ST_tW_top"            , 6952830},
{"ST_tW_antitop"        , 6933094},
{"EWKWPlus2Jets"            ,  4943595},
{"EWKWMinus2Jet"            ,  4782075},
{"EWKZ2Jets_ZToLL"          ,  999451 },
{"EWKZ2Jets_ZToNuNu"        , 2985280 },
{"WGToLNuG"                 ,  27511845},
{"WGstarToLNuMuMu"          ,  2218596},
{"WGstarToLNuEE"            ,  2149729},
};
//signal samples, exclusive ttbar, exclusive diboson samples  can be normalized using GetSumOfWeights()


double getNEventsProcessed(TString filename)
{
  TFile * file = new TFile(filename);
  TH1D * nWeightedEvents = (TH1D*)file->Get("nWeightedEvents");
  if(!nWeightedEvents){
    cout << endl << endl << "Histogram nWeightedEvents doesn't exist in the file "<< filename <<". Quit program." << endl << endl;
    exit(-1);
  }
  double nevents = nWeightedEvents->GetSumOfWeights();
  std::cout<<"nevents= "<<nevents<<std::endl;
  file -> Close();
  delete file;
  return nevents;
}


#endif

void CheckMC(TString era = "2017",
	     TString version = "KIT",
	     double pt_e = 40.0,
	     double eta_e = -0.43,
	     double iso_e = 0.13,
	     double pt_m = 25.,
	     double eta_m = 1.,
	     double iso_m = 0.) {

  //  TString dir1 = "./";
  TString dir1 = "/nfs/dust/cms/user/rasp/CMSSW/Update/CMSSW_10_2_22/src/DesyTauAnalyses/NTupleMaker/data/";
  TString dir2 = "/nfs/dust/cms/user/rasp/CMSSW/Update/CMSSW_10_2_22/src/DesyTauAnalyses/NTupleMaker/data/CorrectionWS_KIT/";
  TString name1 = dir1+"/htt_scalefactors_legacy_"+era+".root";
  TString name2 = dir2+"/htt_scalefactors_legacy_"+era+".root";
  TFile * file1 = new TFile(name1);
  TFile * file2 = new TFile(name2);
  RooWorkspace * w1 = (RooWorkspace*)file1->Get("w");
  RooWorkspace * w2 = (RooWorkspace*)file2->Get("w");

  w1->var("m_pt")->setVal(pt_m);
  w1->var("m_eta")->setVal(eta_m);
  w1->var("m_iso")->setVal(iso_m);
  w1->var("e_pt")->setVal(pt_e);
  w1->var("e_eta")->setVal(eta_e);
  w1->var("e_iso")->setVal(iso_e);

  w2->var("m_pt")->setVal(pt_m);
  w2->var("m_eta")->setVal(eta_m);
  w2->var("m_iso")->setVal(iso_m);
  w2->var("e_pt")->setVal(pt_e);
  w2->var("e_eta")->setVal(eta_e);
  w2->var("e_iso")->setVal(iso_e);

  double isoweight_e_1 = w1->function("e_id90_kit_ratio")->getVal()*w1->function("e_iso_kit_ratio")->getVal();
  double isoweight_m_1 = w1->function("m_looseiso_ic_ratio")->getVal()*w1->function("m_id_kit_ratio")->getVal();
  double isoweight_e_2 = w2->function("e_id90_kit_ratio")->getVal()*w2->function("e_iso_kit_ratio")->getVal();
  double isoweight_m_2 = w2->function("m_looseiso_ic_ratio")->getVal()*w2->function("m_id_kit_ratio")->getVal();

  double eff_data_trig_ehigh_1 = w1->function("e_trg_23_ic_data")->getVal();
  double eff_data_trig_elow_1 = w1->function("e_trg_12_ic_data")->getVal();
  double eff_mc_trig_ehigh_1 = w1->function("e_trg_23_ic_mc")->getVal();
  double eff_mc_trig_elow_1 = w1->function("e_trg_12_ic_mc")->getVal();

  double eff_data_trig_mhigh_1 = w1->function("m_trg_23_ic_data")->getVal();
  double eff_data_trig_mlow_1 = w1->function("m_trg_8_ic_data")->getVal();
  double eff_mc_trig_mhigh_1 = w1->function("m_trg_23_ic_mc")->getVal();
  double eff_mc_trig_mlow_1 = w1->function("m_trg_8_ic_mc")->getVal();

  double eff_data_trig_ehigh_2 = w2->function("e_trg_23_ic_data")->getVal();
  double eff_data_trig_elow_2 = w2->function("e_trg_12_ic_data")->getVal();
  double eff_mc_trig_ehigh_2 = w2->function("e_trg_23_ic_mc")->getVal();
  double eff_mc_trig_elow_2 = w2->function("e_trg_12_ic_mc")->getVal();

  double eff_data_trig_mhigh_2 = w2->function("m_trg_23_ic_data")->getVal();
  double eff_data_trig_mlow_2 = w2->function("m_trg_8_ic_data")->getVal();
  double eff_mc_trig_mhigh_2 = w2->function("m_trg_23_ic_mc")->getVal();
  double eff_mc_trig_mlow_2 = w2->function("m_trg_8_ic_mc")->getVal();

  std::cout << std::endl;
  std::cout << "e_idiso = " << isoweight_e_1 << " : " << isoweight_e_2 << std::endl;
  std::cout << "m_idiso = " << isoweight_m_1 << " : " << isoweight_m_2 << std::endl;
  std::cout << std::endl;
  std::cout << "e_trg_23_data = " << eff_data_trig_ehigh_1 << " : " << eff_data_trig_ehigh_2 << std::endl;
  std::cout << "m_trg_23_data = " << eff_data_trig_mhigh_1 << " : " << eff_data_trig_mhigh_2 << std::endl;
  std::cout << std::endl;
  std::cout << "e_trg_12_data = " << eff_data_trig_elow_1 << " : " << eff_data_trig_elow_2 << std::endl;
  std::cout << "m_trg_8_data  = " << eff_data_trig_mlow_1 << " : " << eff_data_trig_mlow_2 << std::endl;
  std::cout << std::endl;
  std::cout << "e_trg_23_mc = " << eff_mc_trig_ehigh_1 << " : " << eff_mc_trig_ehigh_2 << std::endl;
  std::cout << "m_trg_23_mc = " << eff_mc_trig_mhigh_1 << " : " << eff_mc_trig_mhigh_2 << std::endl;
  std::cout << std::endl;
  std::cout << "e_trg_12_mc = " << eff_mc_trig_elow_1 << " : " << eff_mc_trig_elow_2 << std::endl;
  std::cout << "m_trg_8_mc  = " << eff_mc_trig_mlow_1 << " : " << eff_mc_trig_mlow_2 << std::endl;

}

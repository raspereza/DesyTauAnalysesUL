void CheckEmbed(TString era = "2016",
		double pt_e = 40.0,
		double eta_e = -0.43,
		double iso_e = 0.16,
		double pt_m = 17.,
		double eta_m = 1.1,
		double iso_m = 0.,
		double gt1_pt = 60,
		double gt1_eta = 1.2,
		double gt2_pt = 60,
		double gt2_eta = 1.2) {

  //  TString dir1 = "./";
  TString dir1 = "/nfs/dust/cms/user/rasp/CMSSW/Update/CMSSW_10_2_22/src/DesyTauAnalyses/NTupleMaker/data/";
  TString dir2 = "/nfs/dust/cms/user/rasp/CMSSW/Update/CMSSW_10_2_22/src/DesyTauAnalyses/NTupleMaker/data/CorrectionWS_IC/";
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

  double isoweight_e_1 = w1->function("e_iso_ratio_emb")->getVal()*w1->function("e_id_ratio_emb")->getVal();
  double isoweight_m_1 = w1->function("m_looseiso_ic_embed_ratio")->getVal()*w1->function("m_id_ratio_emb")->getVal();
  double isoweight_e_2 = w1->function("e_idiso_ratio_emb")->getVal();
  double isoweight_m_2 = w1->function("m_idlooseiso_binned_ic_embed_ratio")->getVal();


  std::cout << "e_idiso = " << isoweight_e_1 << " : " << isoweight_e_2 << std::endl;
  std::cout << "m_idiso = " << isoweight_m_1 << " : " << isoweight_m_2 << std::endl;
  std::cout << std::endl;

  w2->var("gt_pt")->setVal(gt1_pt);
  w2->var("gt_eta")->setVal(gt1_eta);
  double id1_embed = w2->function("m_sel_id_ic_ratio")->getVal();
  w2->var("gt_pt")->setVal(gt2_pt);
  w2->var("gt_eta")->setVal(gt2_eta);
  double id2_embed = w2->function("m_sel_id_ic_ratio")->getVal();
  w2->var("gt1_pt")->setVal(gt1_pt);
  w2->var("gt2_pt")->setVal(gt2_pt);
  w2->var("gt1_eta")->setVal(gt1_eta);
  w2->var("gt2_eta")->setVal(gt2_eta);
  double trg_emb = w2->function("m_sel_trg_ic_ratio")->getVal();
  double emWeight_IC = id1_embed * id2_embed * trg_emb;

  w1->var("gt_pt")->setVal(gt1_pt);
  w1->var("gt_eta")->setVal(gt1_eta);
  double id1_embed_x = w1->function("m_sel_idEmb_ratio")->getVal();
  w1->var("gt_pt")->setVal(gt2_pt);
  w1->var("gt_eta")->setVal(gt2_eta);
  double id2_embed_x = w1->function("m_sel_idEmb_ratio")->getVal();
  w1->var("gt1_pt")->setVal(gt1_pt);
  w1->var("gt2_pt")->setVal(gt2_pt);
  w1->var("gt1_eta")->setVal(gt1_eta);
  w1->var("gt2_eta")->setVal(gt2_eta);
  double trg_emb_x = 1.0;
  trg_emb_x = w1->function("m_sel_trg_kit_ratio")->getVal();
  double emWeight_KIT = id1_embed_x * id2_embed_x * trg_emb_x;
  std::cout << std::endl;
  std::cout << "EmbWeight =  KIT : " << emWeight_KIT 
	    << "  IC : " << emWeight_IC << std::endl;


}

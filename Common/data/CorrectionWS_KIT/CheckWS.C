void CheckWS(TString era = "2016",
	     double pt = 20.,
	     double eta = 0.1,
	     double iso = 0.15) {

  TString name = "htt_scalefactors_legacy_"+era+".root";
  TFile * file = new TFile(name);
  RooWorkspace * w = (RooWorkspace*)file->Get("w");
  w->var("m_pt")->setVal(pt);
  w->var("m_eta")->setVal(eta);
  w->var("m_iso")->setVal(iso);
  double m_trg_23_data = w->function("m_trg_23_ic_data")->getVal();
  double m_trg_23_mc = w->function("m_trg_23_ic_mc")->getVal();
  double m_trg_23_emb = w->function("m_trg_23_ic_embed")->getVal();
  double m_trg_8_data = w->function("m_trg_8_ic_data")->getVal();
  double m_trg_8_mc = w->function("m_trg_8_ic_mc")->getVal();
  double m_trg_8_emb = w->function("m_trg_8_ic_embed")->getVal();
  std::cout << "m_trg_23_data = " << m_trg_23_data << std::endl;
  std::cout << "m_trg_23_mc   = " << m_trg_23_mc << std::endl;
  std::cout << "m_trg_23_emb  = " << m_trg_23_emb << std::endl;
  std::cout << "m_trg_8_data = " << m_trg_8_data << std::endl;
  std::cout << "m_trg_8_mc   = " << m_trg_8_mc << std::endl;
  std::cout << "m_trg_8_emb  = " << m_trg_8_emb << std::endl;
  w->var("e_pt")->setVal(pt);
  w->var("e_eta")->setVal(eta);
  w->var("e_iso")->setVal(iso);
  double e_trg_23_data = w->function("e_trg_23_ic_data")->getVal();
  double e_trg_23_mc = w->function("e_trg_23_ic_mc")->getVal();
  double e_trg_23_emb = w->function("e_trg_23_ic_embed")->getVal();
  double e_trg_12_data = w->function("e_trg_12_ic_data")->getVal();
  double e_trg_12_mc = w->function("e_trg_12_ic_mc")->getVal();
  double e_trg_12_emb = w->function("e_trg_12_ic_embed")->getVal();
  std::cout << "e_trg_23_data = " << e_trg_23_data << std::endl;
  std::cout << "e_trg_23_mc   = " << e_trg_23_mc << std::endl;
  std::cout << "e_trg_23_emb  = " << e_trg_23_emb << std::endl;
  std::cout << "e_trg_12_data = " << e_trg_12_data << std::endl;
  std::cout << "e_trg_12_mc   = " << e_trg_12_mc << std::endl;
  std::cout << "e_trg_12_emb  = " << e_trg_12_emb << std::endl;
  double m_idiso = 1.0;
  double m_idiso_emb = 1.0;
  if (era=="2016") {
    m_idiso_emb = w->function("m_looseiso_binned_ic_embed_ratio")->getVal()*
      w->function("m_id_ratio_emb")->getVal();
    m_idiso = w->function("m_idlooseiso_binned_ic_ratio")->getVal();
  }
  if (era=="2017") {
    m_idiso_emb = w->function("m_looseiso_ic_embed_ratio")->getVal()*
      w->function("m_id_embed_kit_ratio")->getVal();
    m_idiso = w->function("m_looseiso_ic_ratio")->getVal()*w->function("m_id_kit_ratio")->getVal();
  }
  if (era=="2018") {
    m_idiso_emb = w->function("m_looseiso_embed_ratio")->getVal()*
      w->function("m_id_embed_kit_ratio")->getVal();
    m_idiso = w->function("m_looseiso_ic_ratio")->getVal()*w->function("m_id_kit_ratio")->getVal();
  }
  std::cout << std::endl;
  std::cout << "m_idiso = " << m_idiso << std::endl;
  std::cout << "m_idiso_emb = " << m_idiso_emb << std::endl;

}

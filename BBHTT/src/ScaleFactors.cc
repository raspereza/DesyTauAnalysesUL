#include "DesyTauAnalyses/BBHTT/interface/ScaleFactors.h"

ScaleFactors::ScaleFactors(RooWorkspace * w_, int era_ , bool isEmbedded_) {

  correctionWS = w_;
  era = era_;
  isEmbedded = isEmbedded_;

  isoweight_1 = 1.0;
  isoweight_2 = 1.0;
  trkeffweight_1 = 1.0;
  trkeffweight_2 = 1.0;
  trig_emu = 1.0;

}

void ScaleFactors::setLep1( double pt1, double eta1, double iso1 ) {
  pt_1 = pt1;
  eta_1 = eta1;
  iso_1 = TMath::Min(iso1,0.149999);
  correctionWS->var("e_pt")->setVal(pt_1);
  correctionWS->var("e_eta")->setVal(eta_1);
  correctionWS->var("e_iso")->setVal(iso_1);
}

void ScaleFactors::setLep2( double pt2, double eta2, double iso2 ) {
  pt_2 = pt2;
  eta_2 = eta2;
  iso_2 = TMath::Min(iso2,0.19999);  
  correctionWS->var("m_pt")->setVal(pt_2);
  correctionWS->var("m_eta")->setVal(eta_2);
  correctionWS->var("m_iso")->setVal(iso_2);

}

void ScaleFactors::setLeptons(double pt1, double eta1, double iso1,
			      double pt2, double eta2, double iso2) {

  setLep1(pt1,eta1,iso1);
  setLep2(pt2,eta2,iso2);
  computeSFs();

}

void ScaleFactors::computeSFs() {

  TString suffix = "mc";
  TString suffixRatio = "ratio";
  if (isEmbedded) {suffix = "embed"; suffixRatio = "embed_ratio";}
 
  isoweight_1 = 1.0;
  isoweight_2 = 1.0;
  trkeffweight_1 = 1.0;
  trkeffweight_2 = 1.0;
  trig_emu = 1.0;

  // id/iso/trk scale factors (from KIT)
  if (era==2016){
    if (isEmbedded) {
      isoweight_1 = correctionWS->function("e_idiso_ratio_emb")->getVal();
      isoweight_2 = correctionWS->function("m_idlooseiso_binned_ic_embed_ratio")->getVal();
    }
    else {
      isoweight_1 = correctionWS->function("e_idiso_ratio")->getVal();
      isoweight_2 = correctionWS->function("m_idlooseiso_binned_ic_ratio")->getVal();
    }
  }
  else{
    if (isEmbedded) {
      if (era==2017) {
	isoweight_1 = correctionWS->function("e_iso_binned_embed_kit_ratio")->getVal()*correctionWS->function("e_id90_embed_kit_ratio")->getVal();
	isoweight_2 = correctionWS->function("m_looseiso_binned_ic_embed_ratio")->getVal()*correctionWS->function("m_id_embed_kit_ratio")->getVal();
      }
      else {
	isoweight_1 = correctionWS->function("e_iso_binned_embed_kit_ratio")->getVal() * correctionWS->function("e_id90_embed_kit_ratio")->getVal();
	isoweight_2 = correctionWS->function("m_looseiso_binned_embed_ratio")->getVal()*correctionWS->function("m_id_embed_kit_ratio")->getVal();
      }
    }
    else {
      isoweight_1 = correctionWS->function("e_id90_kit_ratio")->getVal() * correctionWS->function("e_iso_kit_ratio")->getVal();
      isoweight_2 = correctionWS->function("m_looseiso_ic_ratio")->getVal()*correctionWS->function("m_id_kit_ratio")->getVal();
    }
  }
  if (!isEmbedded){
    if (era == 2018) trkeffweight_1 = correctionWS->function("e_trk_ratio")->getVal();
    if (era==2016 || era==2018) 
      trkeffweight_2 = correctionWS->function("m_trk_ratio")->getVal();
  }
  if (era == 2017) trkeffweight_1 = correctionWS->function("e_trk_ratio")->getVal();

  double eff_data_trig_mhigh = correctionWS->function("m_trg_23_ic_data")->getVal();
  double eff_data_trig_mlow = correctionWS->function("m_trg_8_ic_data")->getVal();
  double eff_mc_trig_mhigh = correctionWS->function("m_trg_23_ic_"+suffix)->getVal();
  double eff_mc_trig_mlow = correctionWS->function("m_trg_8_ic_"+suffix)->getVal();

  double eff_data_trig_ehigh = correctionWS->function("e_trg_23_ic_data")->getVal();
  double eff_data_trig_elow = correctionWS->function("e_trg_12_ic_data")->getVal();
  double eff_mc_trig_ehigh = correctionWS->function("e_trg_23_ic_"+suffix)->getVal();
  double eff_mc_trig_elow = correctionWS->function("e_trg_12_ic_"+suffix)->getVal();

  double eff_emu_data = 
    eff_data_trig_mhigh*eff_data_trig_elow + 
    eff_data_trig_mlow*eff_data_trig_ehigh -
    eff_data_trig_mhigh*eff_data_trig_ehigh;
  double eff_emu_mc = 
    eff_mc_trig_mhigh*eff_mc_trig_elow + 
    eff_mc_trig_mlow*eff_mc_trig_ehigh -
    eff_mc_trig_mhigh*eff_mc_trig_ehigh;	

  if (eff_emu_mc<1e-4||eff_emu_data<1e-4)
    trig_emu = 0.0;
  else
    trig_emu = eff_emu_data/eff_emu_mc;

  if (era==2016) 
    trig_emu *= 0.98;

}

double ScaleFactors::getIdIso1_SF() {
  return isoweight_1;
}

double ScaleFactors::getIdIso2_SF() {
  return isoweight_2;
}

double ScaleFactors::getTrk1_SF() {
  return trkeffweight_1;
}

double ScaleFactors::getTrk2_SF() {
  return trkeffweight_2;
}

double ScaleFactors::getTrigger_SF() {
  return trig_emu;
}

ScaleFactors::~ScaleFactors() {
}

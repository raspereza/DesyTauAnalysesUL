#ifndef ScaleFactors_h
#define ScaleFactors_h

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

class ScaleFactors {
 public:

  ScaleFactors(RooWorkspace * w_, int era_, bool isEmbedded_);
  void setLeptons(double pt1, double eta1, double iso1,
		  double pt2, double eta2, double iso2);
  double getIdIso1_SF();
  double getIdIso2_SF();
  double getTrk1_SF();
  double getTrk2_SF();
  double getTrigger_SF();
  ~ScaleFactors();

 private:

  void setLep1(double pt1, double eta1, double iso1);
  void setLep2(double pt2, double eta2, double iso2);
  void computeSFs();

  RooWorkspace * correctionWS;
  int era;
  bool isEmbedded;
  double pt_1;
  double eta_1;
  double iso_1;
  double pt_2;
  double eta_2;
  double iso_2;

  double trig_emu;
  double isoweight_1;
  double isoweight_2;
  double trkeffweight_1;
  double trkeffweight_2;

};

#endif

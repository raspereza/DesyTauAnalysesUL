#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"

#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TH1D.h>
#include "TCanvas.h"

void Compare()
{
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);

  TFile* f1 = new TFile("PileupMC.root");
  TH1D * histo1 = (TH1D*)f1->Get("MCPU");
  TFile* f2 = new TFile("DataPileupHistogram.root");
  TH1D * histo2 = (TH1D*)f2->Get("pileup");

  TH1D * _NhistoM = new TH1D("NhistoM","NhistoM",100,0,100);
  TH1D * _NhistoD = new TH1D("NhistoD","NhistoD",100,0,100);

  for(int i=1;i<100;i++)
   {
     _NhistoM->SetBinContent(i,histo1->GetBinContent(i));
     _NhistoD->SetBinContent(i,histo2->GetBinContent(i));
   }

  _NhistoM->Scale(1/_NhistoM->GetSumOfWeights());
  _NhistoD->Scale(1/_NhistoD->GetSumOfWeights());

  _NhistoD->SetLineColor(2);
  _NhistoD->Draw("HIST same");
  _NhistoM->Draw("HIST same");

  for(int i=1;i<100;i++)
    cout<<_NhistoD->GetBinContent(i)/_NhistoM->GetBinContent(i)<<",";
  cout<<endl;

}
void Convert() {

  TFile * file_mc = new TFile("PileupMC.root");
  TH1D * histPU_mc = (TH1D*)file_mc->Get("MCPU");
  TFile * file_mc_out = new TFile("pileUp_MC_UL18.root","recreate");
  file_mc_out->cd("");
  histPU_mc->Write("pileup");
  file_mc_out->Close();

  TFile * file_data = new TFile("DataPileupHistogram.root");
  TH1D * histPU_data = (TH1D*)file_data->Get("pileup");
  TFile * file_data_out = new TFile("pileUp_data_UL18.root","recreate");
  file_data_out->cd("");
  histPU_data->Write("pileup");
  file_data_out->Close();



}

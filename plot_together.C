template <class T>
void shift_waveform(T *h, Int_t new_max){
  Int_t npts = h->GetNbinsX();
  Int_t old_max = h->GetMaximumBin();
  Int_t old_ref = old_max - new_max;
  TH1D *htemp = (TH1D*)h->Clone("htemp");
  Double_t temp;
  if(old_ref<0){
    // cout << " case lower" << endl;
    old_ref = npts-(new_max-old_max);
  }
  for(Int_t i = 1; i<npts-(old_ref); i++){
    temp = htemp->GetBinContent(old_ref+i);
    h->SetBinContent(i,temp);
  }
  Int_t aux = 1;
  for(Int_t i = npts-(old_ref); i<=npts; i++){
    temp = htemp->GetBinContent(aux);
    h->SetBinContent(i,temp);
    aux++;
  }
  delete htemp;
}


void plot_together(){
  TFile *fmu = new TFile("ITALIAN_WVF/averaged_waveforms.root","READ");
  TFile *fspe1 = new TFile("ITALIAN_WVF/sphe_waveforms_Ch1.root","READ");
  // TFile *fspe2 = new TFile("ITALIAN_WVF/SPE_HPK_AW.root","READ");
  // TFile *flas = new TFile("ITALIAN_WVF/LASER_HPK_AW.root","READ");
  // TFile *falp = new TFile("ITALIAN_WVF/ALPHA_HPK_AW.root","READ");

  TGraph *gmu = (TGraph*)fmu->Get("average_normalized_ch1");
  TGraph *gspe1 = (TGraph*)fspe1->Get("mean_ch1");
  // TH1D* halp = (TH1D*)falp->Get("Run27_ch0_ADC0_0V_ScintProfFirstSignalBin_5_"); // get Histogram
  // TH1F* hspe2 = (TH1F*)fspe2->Get("Run3_ch0_ADC0_0V_ScintProfFirstSignalBin_17_"); // get Histogram
  
  // Double_t FactorHalp = halp->GetMaximum();
  // Double_t FactorHalp = 1.;
  // Double_t FactorHspe2 = hspe2->GetMaximum();
  // halp->Scale(FactorHalp/halp->GetEntries());
  // hspe2->Scale(FactorHspe2/hspe2->GetEntries());
  // TH1F *hspe2 = (TH1F*)fspe2->Get("Run27_ch0_ADC0_0V_ScintProfFirstSignalBin_5_");
  // TH1D *hlas = (TH1D*)flas->Get("mean_ch1");
  // TH1D *halp = (TH1D*)falp->Get("mean_ch1");
  
  Int_t n = gmu->GetN();
  Double_t dtime = 4e-9;
  
  TH1D *hmu = new TH1D("hmu","hmu",n,0,n*dtime);
  TH1D *hspe1 = new TH1D("hspe1","hspe1",n,0,n*dtime);
  // TH1D *hspe2 = new TH1D("hspe2","hspe2",n,0,n*dtime);
  // TH1D *hlas = new TH1D("hlas","hlas",n,0,n*dtime);
  // TH1D *halp = new TH1D("halp","halp",n,0,n*dtime);

  Double_t normsphe = -1e12;
  for(Int_t i = 0; i<n; i++){
    hmu->SetBinContent(i+1,*(gmu->GetY()+i));
    Double_t tmp = *(gspe1->GetY()+i);
    if(tmp>=normsphe){
      normsphe = tmp;
    }
  }

  for(Int_t i = 0; i<n; i++){
    hspe1->SetBinContent(i+1,*(gspe1->GetY()+i)/normsphe);
  }

  shift_waveform(hspe1,5986/4);

  TFile *f = new TFile("ITALIAN_WVFs.root", "RECREATE");
  hmu->Write();
  hspe1->Write(); // Here write the three histograms on the file and close the file
  f->Close();

  hspe1->SetLineWidth(2);
  hmu->SetLineWidth(2);
  // halp->SetLineWidth(2);
  hspe1->SetLineColor(kBlue);
  // hspe1->SetLineColor(kBlack);
  // halp->SetLineColor(kRed);
  hspe1->SetTitle("spe1");
  hmu->SetTitle("Muon");
  hmu->GetXaxis()->SetTitle("Time in s");
  hmu->GetYaxis()->SetTitle("Amplitude (A.U.)");
  hmu->Draw("same");
  hspe1->Draw("same");
  // halp->Draw(); 
  // hspe2->Draw("same");
  
}


/*
  TFile *fmu = new TFile("averaged_waveforms.root","READ");
  TGraph *gmu = (TGraph*)fmu->Get("average_normalized_ch1");
  Int_t n = gmu->GetN();
  Double_t dtime = 4;
  TH1D *hmu = new TH1D("hmu","hmu",n,0,n*dtime);

  for(Int_t i = 0; i<n; i++){
    hmu->SetBinContent(i+1,*(gmu->GetY()+i));
  }

  hmu->Draw();
*/

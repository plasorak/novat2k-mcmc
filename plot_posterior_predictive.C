
void plot_posterior_predictive() {
  auto file_in = TFile("posterior_predictive.root");
  auto file_out = TFile("for_artur_pp.root", "RECREATE");

  TTree* bins = (TTree*) file_in.Get("bin_values");
  std::vector<double>* nue_fhc    = nullptr;
  std::vector<double>* nue_rhc    = nullptr;
  std::vector<double>* numu_fhc_1 = nullptr;
  std::vector<double>* numu_fhc_2 = nullptr;
  std::vector<double>* numu_fhc_3 = nullptr;
  std::vector<double>* numu_fhc_4 = nullptr;
  std::vector<double>* numu_rhc_1 = nullptr;
  std::vector<double>* numu_rhc_2 = nullptr;
  std::vector<double>* numu_rhc_3 = nullptr;
  std::vector<double>* numu_rhc_4 = nullptr;
  
  bins->SetBranchAddress("nue_fhc_bins"   , &nue_fhc   );
  bins->SetBranchAddress("nue_rhc_bins"   , &nue_rhc   );
  bins->SetBranchAddress("numu_fhc_1_bins", &numu_fhc_1);
  bins->SetBranchAddress("numu_fhc_2_bins", &numu_fhc_2);
  bins->SetBranchAddress("numu_fhc_3_bins", &numu_fhc_3);
  bins->SetBranchAddress("numu_fhc_4_bins", &numu_fhc_4);
  bins->SetBranchAddress("numu_rhc_1_bins", &numu_rhc_1);
  bins->SetBranchAddress("numu_rhc_2_bins", &numu_rhc_2);
  bins->SetBranchAddress("numu_rhc_3_bins", &numu_rhc_3);
  bins->SetBranchAddress("numu_rhc_4_bins", &numu_rhc_4);

  TH1D* data_nue_fhc    = (TH1D*) file_in.Get("asimov_nue_fhc");
  TH1D* data_nue_rhc    = (TH1D*) file_in.Get("asimov_nue_rhc");
  TH1D* data_numu_fhc_1 = (TH1D*) file_in.Get("asimov_numu_fhc_1");
  TH1D* data_numu_fhc_2 = (TH1D*) file_in.Get("asimov_numu_fhc_2");
  TH1D* data_numu_fhc_3 = (TH1D*) file_in.Get("asimov_numu_fhc_3");
  TH1D* data_numu_fhc_4 = (TH1D*) file_in.Get("asimov_numu_fhc_4");
  TH1D* data_numu_rhc_1 = (TH1D*) file_in.Get("asimov_numu_rhc_1");
  TH1D* data_numu_rhc_2 = (TH1D*) file_in.Get("asimov_numu_rhc_2");
  TH1D* data_numu_rhc_3 = (TH1D*) file_in.Get("asimov_numu_rhc_3");
  TH1D* data_numu_rhc_4 = (TH1D*) file_in.Get("asimov_numu_rhc_4");

  TH2D* prediction_nue_fhc    = new TH2D("pp_nue_fhc"   , "#nu_{e} FHC;E_{rec} / Analysis bin;n events", data_nue_fhc   ->GetXaxis()->GetNbins(), data_nue_fhc->GetXaxis()->GetXmin(),data_nue_fhc->GetXaxis()->GetXmax(), 1000, 0., data_nue_fhc   ->GetMaximum()*1.2);
  TH2D* prediction_nue_rhc    = new TH2D("pp_nue_rhc"   , "#nu_{e} RHC;E_{rec} / Analysis bin;n events", data_nue_rhc   ->GetXaxis()->GetNbins(), data_nue_rhc->GetXaxis()->GetXmin(),data_nue_rhc->GetXaxis()->GetXmax(), 1000, 0., data_nue_rhc   ->GetMaximum()*1.2);
  
  TH2D* prediction_numu_fhc_1 = new TH2D("pp_numu_fhc_1", "#nu_{#mu} FHC Q1;E_{rec} GeV;n events", data_numu_fhc_1->GetXaxis()->GetNbins(), data_numu_fhc_1->GetXaxis()->GetXbins()->GetArray(), 1000, 0., data_numu_fhc_1->GetMaximum()*1.2);
  TH2D* prediction_numu_fhc_2 = new TH2D("pp_numu_fhc_2", "#nu_{#mu} FHC Q2;E_{rec} GeV;n events", data_numu_fhc_2->GetXaxis()->GetNbins(), data_numu_fhc_2->GetXaxis()->GetXbins()->GetArray(), 1000, 0., data_numu_fhc_2->GetMaximum()*1.2);
  TH2D* prediction_numu_fhc_3 = new TH2D("pp_numu_fhc_3", "#nu_{#mu} FHC Q3;E_{rec} GeV;n events", data_numu_fhc_3->GetXaxis()->GetNbins(), data_numu_fhc_3->GetXaxis()->GetXbins()->GetArray(), 1000, 0., data_numu_fhc_3->GetMaximum()*1.2);
  TH2D* prediction_numu_fhc_4 = new TH2D("pp_numu_fhc_4", "#nu_{#mu} FHC Q4;E_{rec} GeV;n events", data_numu_fhc_4->GetXaxis()->GetNbins(), data_numu_fhc_4->GetXaxis()->GetXbins()->GetArray(), 1000, 0., data_numu_fhc_4->GetMaximum()*1.2);
  TH2D* prediction_numu_rhc_1 = new TH2D("pp_numu_rhc_1", "#nu_{#mu} RHC Q1;E_{rec} GeV;n events", data_numu_rhc_1->GetXaxis()->GetNbins(), data_numu_rhc_1->GetXaxis()->GetXbins()->GetArray(), 1000, 0., data_numu_rhc_1->GetMaximum()*1.2);
  TH2D* prediction_numu_rhc_2 = new TH2D("pp_numu_rhc_2", "#nu_{#mu} RHC Q2;E_{rec} GeV;n events", data_numu_rhc_2->GetXaxis()->GetNbins(), data_numu_rhc_2->GetXaxis()->GetXbins()->GetArray(), 1000, 0., data_numu_rhc_2->GetMaximum()*1.2);
  TH2D* prediction_numu_rhc_3 = new TH2D("pp_numu_rhc_3", "#nu_{#mu} RHC Q3;E_{rec} GeV;n events", data_numu_rhc_3->GetXaxis()->GetNbins(), data_numu_rhc_3->GetXaxis()->GetXbins()->GetArray(), 1000, 0., data_numu_rhc_3->GetMaximum()*1.2);
  TH2D* prediction_numu_rhc_4 = new TH2D("pp_numu_rhc_4", "#nu_{#mu} RHC Q4;E_{rec} GeV;n events", data_numu_rhc_4->GetXaxis()->GetNbins(), data_numu_rhc_4->GetXaxis()->GetXbins()->GetArray(), 1000, 0., data_numu_rhc_4->GetMaximum()*1.2);

  prediction_nue_fhc   ->SetStats(0);
  prediction_nue_rhc   ->SetStats(0);
  prediction_numu_fhc_1->SetStats(0);
  prediction_numu_fhc_2->SetStats(0);
  prediction_numu_fhc_3->SetStats(0);
  prediction_numu_fhc_4->SetStats(0);
  prediction_numu_rhc_1->SetStats(0);
  prediction_numu_rhc_2->SetStats(0);
  prediction_numu_rhc_3->SetStats(0);
  prediction_numu_rhc_4->SetStats(0);


  for (int i=0; i<bins->GetEntries(); i++) {
    bins->GetEntry(i);
    if (i<2000) continue;
    if (i%10000==0) std::cout << "\033[32mi : " << i << "\033[0m\n";
    for (size_t bx=0; bx<nue_fhc->size(); ++bx) {
      double bc = prediction_nue_fhc->GetXaxis()->GetBinCenter(bx);
      prediction_nue_fhc->Fill(bc,nue_fhc->at(bx));
    }
    for (size_t bx=0; bx<nue_rhc->size(); ++bx) {
      double bc = prediction_nue_rhc->GetXaxis()->GetBinCenter(bx);
      prediction_nue_rhc->Fill(bc,nue_rhc->at(bx));
    }
    for (size_t bx=0; bx<numu_rhc_4->size(); ++bx) {
      double bc = prediction_numu_rhc_4->GetXaxis()->GetBinCenter(bx);
      prediction_numu_rhc_4->Fill(bc,numu_rhc_4->at(bx));
    }
    for (size_t bx=0; bx<numu_rhc_1->size(); ++bx) {
      double bc = prediction_numu_rhc_1->GetXaxis()->GetBinCenter(bx);
      prediction_numu_rhc_1->Fill(bc,numu_rhc_1->at(bx));
    }
    for (size_t bx=0; bx<numu_rhc_2->size(); ++bx) {
      double bc = prediction_numu_rhc_2->GetXaxis()->GetBinCenter(bx);
      prediction_numu_rhc_2->Fill(bc,numu_rhc_2->at(bx));
    }
    for (size_t bx=0; bx<numu_rhc_3->size(); ++bx) {
      double bc = prediction_numu_rhc_3->GetXaxis()->GetBinCenter(bx);
      prediction_numu_rhc_3->Fill(bc,numu_rhc_3->at(bx));
    }
    for (size_t bx=0; bx<numu_fhc_4->size(); ++bx) {
      double bc = prediction_numu_fhc_4->GetXaxis()->GetBinCenter(bx);
      prediction_numu_fhc_4->Fill(bc,numu_fhc_4->at(bx));
    }
    for (size_t bx=0; bx<numu_fhc_1->size(); ++bx) {
      double bc = prediction_numu_fhc_1->GetXaxis()->GetBinCenter(bx);
      prediction_numu_fhc_1->Fill(bc,numu_fhc_1->at(bx));
    }
    for (size_t bx=0; bx<numu_fhc_2->size(); ++bx) {
      double bc = prediction_numu_fhc_2->GetXaxis()->GetBinCenter(bx);
      prediction_numu_fhc_2->Fill(bc,numu_fhc_2->at(bx));
    }
    for (size_t bx=0; bx<numu_fhc_3->size(); ++bx) {
      double bc = prediction_numu_fhc_3->GetXaxis()->GetBinCenter(bx);
      prediction_numu_fhc_3->Fill(bc,numu_fhc_3->at(bx));
    }
  }

  TCanvas c;
  c.Print("posterior_predictive.pdf[");
  prediction_nue_fhc->Draw("COLZ");
  data_nue_fhc->SetLineColor(kBlack);
  data_nue_fhc->SetLineWidth(2);
  data_nue_fhc->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  prediction_nue_rhc->Draw("COLZ");
  data_nue_rhc->SetLineColor(kBlack);
  data_nue_rhc->SetLineWidth(2);
  data_nue_rhc->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  
  prediction_numu_fhc_1->Draw("COLZ");
  data_numu_fhc_1->SetLineColor(kBlack);
  data_numu_fhc_1->SetLineWidth(2);
  data_numu_fhc_1->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  prediction_numu_fhc_2->Draw("COLZ");
  data_numu_fhc_2->SetLineColor(kBlack);
  data_numu_fhc_2->SetLineWidth(2);
  data_numu_fhc_2->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  prediction_numu_fhc_3->Draw("COLZ");
  data_numu_fhc_3->SetLineColor(kBlack);
  data_numu_fhc_3->SetLineWidth(2);
  data_numu_fhc_3->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  prediction_numu_fhc_4->Draw("COLZ");
  data_numu_fhc_4->SetLineColor(kBlack);
  data_numu_fhc_4->SetLineWidth(2);
  data_numu_fhc_4->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  
  prediction_numu_rhc_1->Draw("COLZ");
  data_numu_rhc_1->SetLineColor(kBlack);
  data_numu_rhc_1->SetLineWidth(2);
  data_numu_rhc_1->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  prediction_numu_rhc_2->Draw("COLZ");
  data_numu_rhc_2->SetLineColor(kBlack);
  data_numu_rhc_2->SetLineWidth(2);
  data_numu_rhc_2->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  prediction_numu_rhc_3->Draw("COLZ");
  data_numu_rhc_3->SetLineColor(kBlack);
  data_numu_rhc_3->SetLineWidth(2);
  data_numu_rhc_3->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  prediction_numu_rhc_4->Draw("COLZ");
  data_numu_rhc_4->SetLineColor(kBlack);
  data_numu_rhc_4->SetLineWidth(2);
  data_numu_rhc_4->Draw("SAME");
  c.Print("posterior_predictive.pdf");
  
  c.Print("posterior_predictive.pdf]");


  file_out.cd();  
  prediction_nue_fhc    ->Write();
  data_nue_fhc          ->Write();
  prediction_nue_rhc    ->Write();
  data_nue_rhc          ->Write();
  prediction_numu_fhc_1 ->Write();
  data_numu_fhc_1       ->Write();
  prediction_numu_fhc_2 ->Write();
  data_numu_fhc_2       ->Write();
  prediction_numu_fhc_3 ->Write();
  data_numu_fhc_3       ->Write();
  prediction_numu_fhc_4 ->Write();
  data_numu_fhc_4       ->Write();
  prediction_numu_rhc_1 ->Write();
  data_numu_rhc_1       ->Write();
  prediction_numu_rhc_2 ->Write();
  data_numu_rhc_2       ->Write();
  prediction_numu_rhc_3 ->Write();
  data_numu_rhc_3       ->Write();
  prediction_numu_rhc_4 ->Write();
  data_numu_rhc_4       ->Write();

}

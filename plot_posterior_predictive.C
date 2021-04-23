#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <map>

enum data_sample {
  nue_fhc=0,  
  nue_rhc,
  numu_fhc_1,
  numu_fhc_2,
  numu_fhc_3,
  numu_fhc_4,
  numu_rhc_1,
  numu_rhc_2,
  numu_rhc_3,
  numu_rhc_4,
  last
};

static const std::map<data_sample, std::string> sample_name {{nue_fhc   , "nue_fhc"   },
                                                             {nue_rhc   , "nue_rhc"   },
                                                             {numu_fhc_1, "numu_fhc_1"},
                                                             {numu_fhc_2, "numu_fhc_2"},
                                                             {numu_fhc_3, "numu_fhc_3"},
                                                             {numu_fhc_4, "numu_fhc_4"},
                                                             {numu_rhc_1, "numu_rhc_1"},
                                                             {numu_rhc_2, "numu_rhc_2"},
                                                             {numu_rhc_3, "numu_rhc_3"},
                                                             {numu_rhc_4, "numu_rhc_4"}};

void plot_posterior_predictive(std::string input, std::string output, std::string datafile="", bool first_step_best_fit=false) {
  TFile* file_in = new TFile(input.c_str(), "READ");
  TFile* data_in = nullptr;
  if (datafile!= "")
    data_in = new TFile(datafile.c_str(), "READ");
  TFile* file_out = new TFile((output+".root").c_str(), "RECREATE");

  TTree* bins = (TTree*) file_in->Get("bin_values");

  std::map<data_sample, std::vector<double>*> samples;
  std::map<data_sample, TH1D*> data;
  std::map<data_sample, TH1D*> best_fit;
  std::map<data_sample, TH2D*> predictions;
  
  for (int ii=nue_fhc; ii<last; ii++) {
    data_sample i = static_cast<data_sample>(ii);
    std::string sname = sample_name.at(i);
    
    // deal with ttree branches
    samples[i] = nullptr;
    std::string branch_name = sname+"_bins";
    bins->SetBranchAddress(branch_name.c_str(), &samples[i]);

    std::string asimov_name = "asimov_"+sname;
    // get the data
    if (data_in == nullptr) {
      data[i] = (TH1D*) file_in->Get(asimov_name.c_str());
    } else {
      data[i] = (TH1D*) data_in->Get(sname.c_str());
    }
    
    best_fit[i] = (TH1D*) file_in->Get(asimov_name.c_str());
    best_fit[i]->SetDirectory(0);

    // Add the poissonian error
    for (int ibin=1; ibin<=data[i]->GetXaxis()->GetNbins(); ++ibin) {
      double bc = data[i]->GetBinContent(ibin);
      if (i>2) {
        data[i]->SetBinContent(ibin, bc             /data[i]->GetXaxis()->GetBinWidth(ibin)*0.1);
        data[i]->SetBinError  (ibin, TMath::Sqrt(bc)/data[i]->GetXaxis()->GetBinWidth(ibin)*0.1);
      } else {
        data[i]->SetBinContent(ibin, bc             );
        data[i]->SetBinError  (ibin, TMath::Sqrt(bc));
      }
      best_fit[i]->SetBinContent(ibin, 0.);
      best_fit[i]->SetBinError  (ibin, 0.);
    }

    // create the th2d for the posterior predictive
    std::string histname = "pp_"+sname;
    int b = data[i]->GetMaximumBin();
    double max = (data[i]->GetBinContent(b)+data[i]->GetBinError(b))*1.2;

    if (sname.find("nue")!=std::string::npos) {
      std::string hc = sname.substr(4,3);
      std::transform(hc.begin(), hc.end(),hc.begin(), ::toupper);
      
      std::string title = "#nu_{e} "+hc+";E_{rec} / Analysis bin;Events";
      
      predictions[i] = new TH2D(histname.c_str(), title.c_str(),
                                data[i]->GetXaxis()->GetNbins(), data[i]->GetXaxis()->GetXmin(), data[i]->GetXaxis()->GetXmax(), 100, 0., max);
    } else {
      std::string q  = sname.substr(9,1);
      std::string hc = sname.substr(5,3);
      std::transform(hc.begin(), hc.end(),hc.begin(), ::toupper);

      std::string title = "#nu_{#mu} " +hc+" Q"+q+";E_{rec} GeV;Events / 0.1 GeV";
      
      predictions[i] = new TH2D(histname.c_str(), title.c_str(),
                                data[i]->GetXaxis()->GetNbins(), data[i]->GetXaxis()->GetXbins()->GetArray(), 100, 0., max);
    }
    predictions[i]->SetStats(0);
  }
  
  std::cout << "\033[32mdone initiating\033[0m\n";
  for (int i=0; i<bins->GetEntries(); i++) {
    bins->GetEntry(i);
    if (i==0 and first_step_best_fit) {
      for (int ii=nue_fhc; ii<last; ii++) {
        data_sample i = static_cast<data_sample>(ii);
        for (size_t bx=0; bx<samples[i]->size(); ++bx) {
          double bc = best_fit[i]->GetXaxis()->GetBinCenter(bx);
          if (i>2) {
            best_fit[i]->SetBinContent(bx,samples[i]->at(bx)/data[i]->GetXaxis()->GetBinWidth(bx)*0.1);
          } else {
            best_fit[i]->SetBinContent(bx,samples[i]->at(bx));
          }
        }
      }
      continue;
    }
    
    if (i%10000==0) std::cout << "\033[32mi : " << i << "\033[0m\n";

    for (int ii=nue_fhc; ii<last; ii++) {
      data_sample i = static_cast<data_sample>(ii);
      for (size_t bx=0; bx<samples[i]->size(); ++bx) {
        double bc = predictions[i]->GetXaxis()->GetBinCenter(bx);
        if (i>2) {
          predictions[i]->Fill(bc,samples[i]->at(bx)/data[i]->GetXaxis()->GetBinWidth(bx)*0.1);
        } else {
          predictions[i]->Fill(bc,samples[i]->at(bx));
        }
      }
    }
  }

  TCanvas c;
  gPad->SetRightMargin(1.2*gPad->GetRightMargin());
  c.Print((output+".pdf[").c_str());
  for (int ii=nue_fhc; ii<last; ++ii) {
    data_sample i = static_cast<data_sample>(ii);
    CenterTitles(predictions[i]);
    predictions[i]->Draw("COLZ");
    data       [i]->SetLineColor  (kOrange-6);
    data       [i]->SetMarkerColor(kOrange-6);
    data       [i]->SetLineStyle(2);
    data       [i]->SetLineWidth(2);
    data       [i]->Draw("SAME");
    if (first_step_best_fit) {
      best_fit   [i]->SetLineColor  (kRed);
      best_fit   [i]->SetMarkerColor(kRed);
      best_fit   [i]->SetLineWidth(2);
      best_fit   [i]->Draw("SAME");
    }
    gPad->RedrawAxis();
    c.Print((output+".pdf").c_str());
  }
  c.Print((output+".pdf]").c_str());

  file_out->cd();
  for (int ii=nue_fhc; ii<last; ++ii) {
    data_sample i = static_cast<data_sample>(ii);
    predictions[i]->Write();
    data       [i]->Write();
  }
  file_out->Close();

}

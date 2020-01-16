
#include "feu.h"
#include "Event.h"


int Run_feu_File_Series()
{
  gSystem->Load("libfeu.so");
  
  TH1F *h_cent;
  TH1F *hcent0 = new TH1F("hcent0", "Cumulative Centroid Distr. (Centered on Zero)", 10000, -50, 50);
  hcent0->GetXaxis()->SetTitle("Centroid (mm)");
  hcent0->GetYaxis()->SetTitle("Counts");
  
  float hcent_mean_bin = 0.0;

  char rootline[256];

  for(int i=0; i<102; i++)
    {
      if (i==1 || i==74 || i==86) continue;
      int run = 2100 + i;
      sprintf(rootline, "pfileopen(\"/data/eic/bobs_lab/pos_scan/dream-0000%d-0000.evt\")", run);
      cout<<rootline<<endl;
      gROOT->ProcessLine(rootline);
      gROOT->ProcessLine("prun()");
      gROOT->ProcessLine("pclose()");
      
      //h_cent->Reset();
      h_cent = (TH1F*)gDirectory->Get("hcent");
      h_cent->Fit("gaus");
      float mean = h_cent->GetFunction("gaus")->GetParameter(1);
      //hcent_mean_bin = h_cent->FindBin(h_cent->GetMean()); cout<<"cent_mean: "<<hcent_mean_bin<<endl;
      hcent_mean_bin = h_cent->FindBin(mean); cout<<"cent_mean: "<<hcent_mean_bin<<endl;

      for (int i=1; i<h_cent->GetNbinsX()-2; i++)
	{
	  //if (hcent_mean_bin>=0) hcent0->SetBinContent(i,h_cent->GetBinContent(i+(hcent_mean_bin-5000)));
	  int mbin = i+(hcent_mean_bin-5000); //this bin value should not exceed h_cent histo bin limits
	  if (mbin >= 0 && mbin <= h_cent->GetNbinsX()) hcent0->Fill(h_cent->GetBinCenter(i),h_cent->GetBinContent(mbin));
	}
  
      // TCanvas *can1 = new TCanvas("can1", "can1");
      // can1->Divide(2,1);
      // can1->cd(1);
      // TH1F *h_sed_100 = (TH1F*)gDirectory->Get("sed_100");
      // h_sed_100->Draw("colz");
      // can1->cd(2);
      // TH1F *h_sed_proj_100 = (TH1F*)gDirectory->Get("sed_proj_100");
      // h_sed_proj_100->Draw();
      
      
      TCanvas *can2 = new TCanvas ("can2", "can2");
      can2->Divide(2,2);
      can2->cd(1);
      h_cent->Draw();
      can2->cd(2);
      hcent0->Draw();
      can2->cd(3);
      TH1F *h_corr = (TH1F*)gDirectory->Get("hcorr");
      h_corr->Draw("colz");
      can2->cd(4);
      TH1F *h_res_corr = (TH1F*)gDirectory->Get("hres_corr");
      h_res_corr->Draw("colz");
    }
  return 0;
}



{
  gROOT->ProcessLine("#include \"feu.h\" ");
  gSystem->Load("libfeu.so");
  pfileopen( "/data/eic/bobs_lab/pos_scan/dream-00002190-0000.evt" );
  prun(2000);

  TCanvas *can = new TCanvas("main","main",800,800);
  if(0) { // cnt_debug
    TH2F *feu100 = (TH2F*) gDirectory->FindObject("sed_100");
    can->cd(1);
    feu100->Draw("lego");
  } else {
    can->Divide(3,3);
    TH1F *slf = (TH1F*) gDirectory->FindObject("slf");
    TH2F *sed = (TH2F*) gDirectory->FindObject("sed");
    TH2F *ped = (TH2F*) gDirectory->FindObject("ped");
    TH2F *sgn = (TH2F*) gDirectory->FindObject("sgn");
    TH1F *max = (TH1F*) gDirectory->FindObject("max");
    TH1F *ampl = (TH1F*) gDirectory->FindObject("ampl");
    TH1F *gx = (TH1F*) gDirectory->FindObject("gx");
    can->cd(1); slf->Draw();
    can->cd(2); sed->Draw("colz");
    can->cd(3); ped->Draw("colz");
    can->cd(4); sgn->Draw("colz");
    can->cd(5); max->Draw();
    can->cd(6); ampl->Draw();
    can->cd(7); gx->Draw();
    can->cd(8); ampl->Draw();
    TF1 *fit = new TF1("fit","[0]*TMath::Gaus(x,[1],[2],1)");
    gx->Fit(fit);
    cout << " "<< fit->GetParameter(1) << " " << fit->GetParameter(2) << endl;
  }
}



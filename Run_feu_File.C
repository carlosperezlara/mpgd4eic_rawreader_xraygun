
//#include "feu.h"  

{
  gROOT->ProcessLine("#include \"feu.h\" ");
  gSystem->Load("libfeu.so");
  pfileopen("/data/eic/bobs_lab/pos_scan/dream-00002190-0000.evt");
  prun(2);

  TCanvas can;
  can.Divide(2,1);
  can.cd(1);
  sed_100->Draw("colz");
  can.cd(2);
  sed_proj_100->Draw();
  
}



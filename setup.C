{
  int xoffset=550;
  int yoffset=0;
  TCanvas *c0 = new TCanvas("c0", "c0", 50,50,400,500);

  c0->Divide(1,2);
  c0->cd(1)->SetLogy();
  h2->Draw();
  c0->cd(2)->SetLogy();
  h4->Draw();
  sleep(3);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 50+xoffset, 50+yoffset, 600, 400);
  c1->Divide(2,1);
  c1->cd(1);
  sed_100->Draw("Lego2");
  c1->cd(2);
  sed_119->Draw("Lego2");


  pupdate(c0,10);
  pupdate(c1,5);


}

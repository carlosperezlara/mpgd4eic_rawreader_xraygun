

int Run_feu_File(int filename_num, int Nevt)
{
  gROOT->ProcessLine("#include \"feu.h\" ");
  
  gSystem->Load("libfeu.so");

  
  char filename[256];
  sprintf(filename, "/data/eic/bobs_lab/pos_scan/dream-0000%d-0000.evt", filename_num);             

  cout<<"Filename: "<<filename<<endl;

  char rootline[256];

  sprintf(rootline, "pfileopen(\"%s\")", filename);
  gROOT->ProcessLine(rootline);

  sprintf(rootline, "prun(%d)",Nevt);
  gROOT->ProcessLine(rootline);

  gROOT->ProcessLine("sed_100->Draw(\"colz\");");
  
  return 0;
}



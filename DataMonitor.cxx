#include <TApplication.h>
#include <TGFrame.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TStyle.h>
#include <TGTextEdit.h>
#include <TGNumberEntry.h>
#include <TGIcon.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>

#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TGraph.h>
#include <TTimer.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH2D.h>
#include <TTimeStamp.h>
#include <TLatex.h>
#include <TImage.h>

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>

#include <MappingTableCollection.h>
#include "DataMonitor.h"

const int __window_wd__ = 1200;
const int __window_ht__ = 600;

ClassImp(DataMonitor);

//====================
void DataMonitor::Merge() {
  // Merges event-by-event information stored
  // in fChannel into persistent histograms
  if(!fReady) return; // executes only if all
  if(fLearning) {     // objects were created
    // loading up pedestal-subtracted waveforms
    for(int i=0; i!=kNumberOfBoards; ++i) {
      int binsch = fScan[i]->GetXaxis()->GetNbins();
      int binssa = fScan[i]->GetYaxis()->GetNbins();
      int kmin = fIntMean[i] - fIntHWindow[i]-1;
      int kmax = fIntMean[i] + fIntHWindow[i]+1;
      int kwid = kmax-kmin-1;
      for(int j=0; j!=binsch; ++j) {
	double ped = fPedestals[i][j];
	if(1) {
	  ped=0;
	  for(int k=kNumberOfSamples-fIntHWindow[i]; k!=kNumberOfSamples; ++k) {
	    ped += fScan[i]->GetBinContent(j+1,k+1);
	  }
	  ped /= fIntHWindow[i];
	}
	double sum = 0;
	//std::cout << "DBG " << j << ": " << ped << " " << kmin << " " << kmax << " " << kwid << std::endl;
	for(int k=0; k!=binssa; ++k) {
	  double val = fScan[i]->GetBinContent(j+1,k+1);
	  fScanAdc[i]->Fill(double(j),double(k),val);
	  if(j==75) fScanWav[i]->Fill(double(k),val);
	  val -= ped;
	  fScanPed[i]->Fill(double(j),double(k),val);
	  if((k>kmin)&&(k<kmax)) {
	    sum += val;
	  }
	}
	fScanSgn[i]->Fill(double(j),sum);
      }
    }
    if(fNoEventsSampled>kNumberOfEventsLearning) {
      fLearning = kFALSE;
    }
    if(fNoEventsSampled%kMergeRefresh==0) {
      RefreshAll();
      std::cout << fNoEventsSampled << std::endl;
    }
  }

  return;

  std::vector<std::pair<unsigned, unsigned> > myhits;
  for(int i=0; i!=kNumberOfBoards; ++i) {
    for(int j=0; j!=fDREAMChannels[i]; ++j) {
      Int_t bpea = fChannel[i][j]->GetMaximumBin();
      Int_t bmin = bpea-fIntHWindow[i];
      Int_t bmax = bpea+fIntHWindow[i];
      if(bmin<0) bmin=0;
      if(bmax>kNumberOfSamples) bmax=kNumberOfSamples;
      Double_t sum = fChannel[i][j]->Integral( bmin, bmax  );
      sum -= (bmax-bmin)*fPedestals[i][j];
      //Double_t hei = fChannel[i][j]->GetMaximum() - fPedestals[i][j];
      Double_t hei = fChannel[i][j]->GetBinContent( bpea ) - fPedestals[i][j];
      fHeight[i][j]->Fill( hei );
      fWidth[i][j]->Fill( sum );
      for(int b=0; b!=kNumberOfSamples; ++b) {
	Double_t xc = fChannel[i][j]->GetXaxis()->GetBinCenter(b+1);
	Double_t yc = fChannel[i][j]->GetBinContent(b+1);
	fSignal[i][j]->Fill(xc,yc,hei);
	fTimeSummary[i]->Fill(double(b), yc );
      }
      fHitSummary[i]->Fill(double(j-fDREAMChannels[i]/2),sum);
      myhits.push_back( std::make_pair( fDREAMChannel[i][j], sum ) );
    }
    //std::cout << " myhits " << myhits.size() << " || ";
    /*
    for(unsigned i=0; i!=myhits.size(); ++i) {
      std::cout << myhits[i].first << " ";
      std::cout << myhits[i].second << " ";
      std::cout << " || ";
    }
    */
    //std::cout << std::endl;
    
    const MappingPlane *mplane = GetMappingPlane(i, GetLocalCell(i) );
    if(!mplane) continue;
    std::vector<Cluster> myclusters = mplane->GetClusters( myhits );
      
    //std::cout << " BOARD " << i << " || myclusters " << myclusters.size() << " || WD: ";
    fClusters_Num[i]->Fill( double(myclusters.size()) );
    for(unsigned kk=0; kk!=myclusters.size(); ++kk) {
      fClusters_Wid[i]->Fill( myclusters[kk].mWidth );
      fClusters_Amp[i]->Fill( myclusters[kk].mAmplitude );
      fClusters_xCE[i]->Fill( myclusters[kk].mCentroid );
      //std::cout << myclusters[kk].mWidth << " AMP: ";
      //std::cout << myclusters[kk].mAmplitude << " xCENT: ";
      //std::cout << myclusters[kk].mCentroid << " ";
      //std::cout << " || ";
    }
    //std::cout << std::endl;
  }
  if(fNoEventsSampled%kMergeRefresh==0) {
    RefreshAll();
  }
}
//====================
TString DataMonitor::GetLocalCell(Int_t bd) {
  // returns the cell name as
  // string from board index
  Double_t x, y, xp, yp;
  GetXYFromCell(fClosestCell,x,y);
  TString cell = GetCellName(x,y);
  return cell;
}
//====================
void DataMonitor::StampRun(TGCompositeFrame *mf) {
  TGLabel *labs;
  fThisRun = new TGLabel(mf, "Run : -1" );
  //TGFont *font = fThisRun->GetFont();
  mf->AddFrame(fThisRun, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX, 5, 5, 2, 2));
  labs = new TGLabel(mf, Form("Position : ( %.1f, %.1f )", fPosX, fPosY) );
  mf->AddFrame(labs, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX, 5, 5, 2, 2));
  labs = new TGLabel(mf, Form("Closest Cell : %s", fClosestCell.Data()) );
  mf->AddFrame(labs, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX, 5, 5, 2, 2));
  fEventsSampled = new TGLabel(mf, "Events sampled: ---" );
  mf->AddFrame(fEventsSampled, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX, 5, 5, 2, 2));
  fSamplingFraction = new TGLabel(mf, "Sampling fraction: ----" );
  mf->AddFrame(fSamplingFraction, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX, 5, 5, 2, 2));
}
//====================
void DataMonitor::NewRun(Int_t run) {
  // FIRST SAVE
  static bool awaking = true;
  if(!awaking) {
    // saving info
    std::cout << " * NEW RUN saving statistics of previus run" << std::endl;
    TTimeStamp timestamp;
    TString ts = timestamp.AsString("s");;
    std::ofstream fout("Run_Data/history.log",std::ofstream::out|std::ofstream::app);
    fout << ts.Data() << "  " << kNumberOfBoards << "  ";
    for(int bd=0;bd!=kNumberOfBoards; ++bd) {
      Double_t SUM = 0;
      for(int i=0; i!=fDREAMChannels[bd]; ++i) {
	SUM += fWidth[bd][i]->Integral();
      }
      fout << Form("  %e",SUM);
    }
    fout << std::endl;
    fout.close();
  } else {
    awaking = false;
  }

  // SECOND CLEAR
  std::cout << " * NEW RUN reset histograms" << std::endl;
  for(int i=0; i!=kNumberOfBoards; ++i) {
    fHitSummary[i]->Reset();
    for(int j=0; j!=kNumberOfChannels; ++j) {
      fHeight[i][j]->Reset();
      fWidth[i][j]->Reset();
    }
  }
  fNoEventsSampled = 0;
  ReadPosition();
  //ReadConfig();
  //ConfigureChannels();

  return;
  //THIRD REFRESH
  std::cout << " * NEW RUN refresh" << std::endl;
  if(fThisRun) fThisRun->SetText( Form("Run : %d",run) );
  std::ifstream fin("Run_Data/history.log");
  TString timeS1, timeS2;
  Int_t numb;
  Double_t ddd[kNumberOfBoards];
  for(;;) {
    fin >> timeS1 >> timeS2 >> numb;
    if(!fin.good()) break;
    for(int i=0; i<numb; ++i) {
      fin >> ddd[i];
    }
    if(!fin.good()) break;
    for(int i=0; i<numb; ++i) {
      fHistory[i]->SetPoint( fHistory[i]->GetN(), fHistory[i]->GetN(), ddd[i] );
    }
  }
  fin.close();
}
//====================
void DataMonitor::NewEvent(Int_t evr) {
  //fReady = true;
  for(int i=0; i!=kNumberOfBoards; ++i) {
    for(int j=0; j!=kNumberOfChannels; ++j) {
      fScan[i]->Reset();
      if(!fChannel[i][j]) {
	fReady = false;
	break;
      }
      fChannel[i][j]->Reset();
    }
    if(!fReady) break;
  }
  //std::cout << "NewEvent" << std::endl;
  if(!fReady) {
    std::cout << "DM: not ready to take data now." << std::endl;
    return;
  }
  fThisEvent = evr;
  fEventsReaded += 1;
  fNoEventsSampled += 1;
  return;
  static bool __up__ = false;
  //std::cout << "AA " << evr << std::endl;
  if(fNoEventsSampled%kNewEventRefresh==0) {
    TString tmp = Form("Events sampled : %d",fNoEventsSampled);
    fEventsSampled->SetText( tmp.Data() );
    if(fThisEvent>100) {
      TString tmp = Form("Sampling fraction : %.3f",fNoEventsSampled/double(fThisEvent));
      fSamplingFraction->SetText( tmp.Data() );
    }
    if(!__up__) {
      __up__ = true;
      std::cout << "DM: I am ready."<< std::endl;
    }
  }
  //std::cout << "BB " << evr << std::endl;
}
//====================
void DataMonitor::CreateClusters(TGCompositeFrame *mf) {
  TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("clustersplot",mf,
								500,500,
								kSunkenFrame);
  Int_t cId = embeddedCanvas->GetCanvasWindowId();
  fCanvasMapCL = new TCanvas("CanvasMapCL", 10, 10, cId);
  embeddedCanvas->AdoptCanvas(fCanvasMapCL);
  mf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasMapCL->SetTopMargin(0.1);
  fCanvasMapCL->SetBottomMargin(0.1);
  fCanvasMapCL->SetLeftMargin(0.2);
  fCanvasMapCL->SetRightMargin(0.1);
  fCanvasMapCL->Divide(4,kNumberOfBoards);
  for(int bd=0; bd!=kNumberOfBoards; ++bd) {
    TVirtualPad *tmp1 = fCanvasMapCL->cd(bd*4+1);
    fClusters_Num[bd]->Draw("HIST");
    TVirtualPad *tmp2 = fCanvasMapCL->cd(bd*4+2);
    fClusters_Amp[bd]->Draw("HIST");
    TVirtualPad *tmp3 = fCanvasMapCL->cd(bd*4+3);
    fClusters_Wid[bd]->Draw("HIST");
    TVirtualPad *tmp4 = fCanvasMapCL->cd(bd*4+4);
    fClusters_xCE[bd]->Draw("HIST");
  }
  fCanvasMapCL->SetEditable(kFALSE);
}
//====================
void DataMonitor::CreateChannels(TGCompositeFrame *mf) {
  TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("channelsplot",mf,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId = embeddedCanvas->GetCanvasWindowId();
  fCanvasMap = new TCanvas("CanvasMap", 10, 10, cId);
  embeddedCanvas->AdoptCanvas(fCanvasMap);
  mf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasMap->SetTopMargin(0.1);
  fCanvasMap->SetBottomMargin(0.1);
  fCanvasMap->SetLeftMargin(0.2);
  fCanvasMap->SetRightMargin(0.1);
  fCanvasMap->Divide(kNumberOfChannels,kNumberOfBoards,0,0);
  //TH2D *axis = new TH2D("axis","",100,-15,270,100,-50,1099);
  TH2D *axis = new TH2D("axis","",100,-kNumberOfSamples*0.1,kNumberOfSamples*1.1,100,-50,1099);
  axis->SetStats(0);
  axis->GetYaxis()->SetLabelSize(0.1);
  axis->GetXaxis()->SetLabelSize(0.1);
  //axis->SetTitleSize(50);
  axis->GetXaxis()->SetNdivisions(508);
  axis->GetYaxis()->SetNdivisions(508);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    for(int j=0; j!=kNumberOfChannels; ++j) {
      TVirtualPad *tmp = fCanvasMap->cd(i*kNumberOfChannels+j+1);
      axis->SetTitle( Form("S%dP%d",i,j) );
      axis->DrawCopy();
      fChannel[i][j]->Draw("BSAME");
    }
  }
  fCanvasMap->SetEditable(kFALSE);
}
//====================
void DataMonitor::ReadPosition() {
  std::ifstream fin("./Position_Data/Last.log");
  fin >> fClosestCell;
  fin >> fPosX;
  fin >> fPosY;
  fin.close();
  std::cout << " * Position Read: " << fPosX << " " << fPosY << " " << fClosestCell.Data() << std::endl;
}
//====================
void DataMonitor::ReadConfig() {
  std::ifstream fin("./config/constants.dat");
  std::cout << " * reading ./config/constants.dat" << std::endl;
  TString var, val;
  for(;;) {
    fin >> var >> val;
    if(!fin.good()) break;
    if(var.Contains("kNumberOfBoards")) {
      kNumberOfBoards = val.Atoi();
    } else if(var.Contains("kTotalNumberOfChannels")) {
      kTotalNumberOfChannels = val.Atoi();
    } else if(var.Contains("kNumberOfChannels")) {
      kNumberOfChannels = val.Atoi();
    } else if(var.Contains("kNumberOfSamples")) {
      kNumberOfSamples = val.Atoi();
    } else if(var.Contains("kNumberOfEventsLearning")) {
      kNumberOfEventsLearning = val.Atoi();
    } else if(var.Contains("kMergeRefresh")) {
      kMergeRefresh = val.Atoi();
    } else if (var.Contains("kNewEventRefresh")) {
      kNewEventRefresh = val.Atoi();
    }
  }
  fin.close();
  std::cout << "  -kNumberOfBoards " << kNumberOfBoards << std::endl;
  std::cout << "  -kTotalNumberOfChannels " << kTotalNumberOfChannels << std::endl;
  std::cout << "  -kNumberOfChannels " << kNumberOfChannels << std::endl;
  std::cout << "  -kNumberOfSamples " << kNumberOfSamples << std::endl;
  std::cout << "  -kNumberOfEventsLearning " << kNumberOfEventsLearning << std::endl;
  std::cout << "  -kMergeRefresh " << kMergeRefresh << std::endl;
  std::cout << "  -kNewEventRefresh " << kNewEventRefresh << std::endl;

  fin.open( "./config/boards.dat" );
  std::cout << " * reading ./config/boards.dat" << std::endl << " ";
  for(int nlines=0;nlines!=kMaxNumberOfBoards;++nlines) {
    Int_t bdid;
    TString bdtype, bdtec;
    fin >> bdid >> bdtype >> bdtec;
    if(!fin.good()) break;
    if(bdid<0||bdid>kNumberOfBoards) continue;
    fBoardCode[bdid] = bdtype;
    fBoardTech[bdid] = bdtec;
  }
  for(int i=0; i!=kNumberOfBoards; ++i) {
    std::cout << " | " << i << " " << fBoardCode[i] << " " << fBoardTech[i];
  }
  fin.close();
  std::cout << std::endl;

  for(int i=0; i!=kNumberOfBoards; ++i) {
    fin.open( Form("./ped/board%d.ped",i) );
    std::cout << Form(" * reading ./ped/board%d.ped",i) << "... ";
    int nread;
   for(nread=0;nread!=kTotalNumberOfChannels;++nread) {
      fin >> fPedestals[i][nread];
      if(!fin.good()) break;
    }
    fin.close();
    std::cout << nread << " values found" << std::endl;
  }

}
//====================
void DataMonitor::ConfigureChannels() {
  // read geometrical configuration of cell
  // from board database
  // needs: GetLocalCell(), GetMappingPlane()
  // returns: fDREAMChannel index, fPitchX
  for(int bd=0; bd!=kNumberOfBoards; ++bd) {
    TString cell = GetLocalCell(bd);
    fBoardCELL[bd] = cell;
    if(fDiagrams[bd])
      fDiagrams[bd]->SetTitle( Form("[%d]%s  %s;X [mm];Y [mm]",bd,fBoardTech[bd].Data(),cell.Data()) );
    const MappingPlane *mplane = GetMappingPlane(bd, cell);
    if(!mplane) continue;
    fDREAMChannels[bd] = mplane->GetStripCount();
    if(fDREAMChannels[bd]>kNumberOfChannels) {
      std::cout << "NO NO " << fDREAMChannels[bd] << " " << kNumberOfChannels << std::endl;
      std::cout << "at board " << bd << std::endl;
      exit(0);
    }
    //std::cout << " Channels in board number " << bd << " : " << fDREAMChannels[bd] << std::endl;
    for(int ch=0; ch!=fDREAMChannels[bd]; ++ch) {
      fDREAMChannel[bd][ch] = mplane->GetDreamChannel(ch);
      fPitchX[bd][ch] = mplane->GetPitch();
      fPeriodY[bd][ch] = 1;
      fStretch[bd][ch] = 1;
      //std::cout << fDREAMChannel[bd][ch] << " ";
    }
    //std::cout << std::endl;
  }
}
//====================
const MappingPlane* DataMonitor::GetMappingPlane(Int_t bd, TString sp, Int_t pl) {
  MappingTable *mtable = NULL;
  mtable = fMapCollection->GetMappingTable( fBoardCode[bd].Data()  );
  if(!mtable) {
    //std::cout << " Could not find board code " << bd << " " << fBoardCode[bd].Data() << " in collection." << std::endl;
    fNotInstalled[bd] = true;
  }
  if(fNotInstalled[bd]) return NULL;
  fNotInstalled[bd] = false;
  MappingSpot *mspot = NULL;
  mspot = mtable->GetMappingSpot(sp.Data()); assert(mspot);
  if(!mspot) {
    std::cout << " Could not find cell code " << sp.Data() << " in table. abort" << std::endl;
    exit(0);
  }
  return mspot->GetMappingPlane(pl);
}
//====================
void DataMonitor::ModelPads(Int_t bd) {
  double x[100];
  double y[100];
  double gx = -fPitchX[bd][0]*fDREAMChannels[bd]/2; //0;
  double gy = -fPitchX[bd][0]*fDREAMChannels[bd]/2; //0;
  int widY = fPitchX[bd][0]*fDREAMChannels[bd];
  for(int str=0; str!=fDREAMChannels[bd]; ++str) {
    double xp = fPitchX[bd][str];
    double wd = 0.8*xp;
    double yp = fPeriodY[bd][str];
    double st = fStretch[bd][str];
    int nvert = (widY/yp*2+1)*2;
    for(int i=0; i!=nvert; ++i) {
      double sx, sy;
      if(i<(nvert/2)) {
	sx = (i%2)*st;
	sy = i*yp/2;
      } else {
	sx = ((i+1)%2)*st;
	sy = (nvert-i-1)*yp/2;
      }
      sx += (i/(nvert/2))*wd;
      x[i] = gx + sx;
      y[i] = sy + gy;
      //cout << x[i] << " | " << y[i] << endl;
    }
    fDiagrams[bd]->AddBin(nvert,x,y);
    gx += xp;
  }
}
//====================
void DataMonitor::GetXYFromCell(TString mycell, Double_t &x, Double_t &y) {
  char cellstr[10] = {'A','B','C','D','E',
		      'F','G','H','I','J'};
  const char tmpC = mycell[0];
  TString tmpR = mycell[1];
  int row = tmpR.Atoi();
  int col = 0;
  for(int i=0; i!=10; ++i) {
    if(tmpC==cellstr[i])
      col = i;
  }
  x = -45 + 10*col;
  y = -45 + 10*row;
  return;
}
//====================
TString DataMonitor::GetCellName(Double_t x, Double_t y) {
  TString cell;
  char cellstr[10] = {'A','B','C','D','E',
		      'F','G','H','I','J'};
  int xx = x + 50;
  int yy = y + 50;
  if(xx<0) xx=0;
  if(yy<0) yy=0;
  xx = xx/10;
  yy = yy/10;
  if(xx>9) xx=9;
  if(yy>9) yy=9;
  cell = Form("%c%d",cellstr[xx],yy);
  return cell;
}
//====================
void DataMonitor::CreateDisplayConfiguration(TGCompositeFrame *mf) {
  TGCompositeFrame *mfR1 = new TGCompositeFrame(mf, 170, 20, kVerticalFrame);
  mf->AddFrame(mfR1, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
  TGTextEdit *fConfigurationFile = new TGTextEdit(mfR1,0.95*__window_wd__, 0.2*__window_ht__);
  mfR1->AddFrame(fConfigurationFile, new TGLayoutHints(kLHintsCenterX|kLHintsExpandX, 5, 5, 2, 2));
  TString sConfFile = "./config/boards.dat";
  fConfigurationFile->LoadFile(sConfFile.Data());
  fConfigurationFile->GetText();
  fConfigurationFile->Goto( fConfigurationFile->GetText()->RowCount() , 0 );

  TGCompositeFrame *mfR2 = new TGCompositeFrame(mf, 170, 20, kVerticalFrame);
  mf->AddFrame(mfR2, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
  TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("config file",mfR2,0.95*__window_wd__, 0.4*__window_ht__,kSunkenFrame);
  Int_t cId = embeddedCanvas->GetCanvasWindowId();
  fCanvasMapCF = new TCanvas("CanvasMapCF", 10, 10, cId);
  embeddedCanvas->AdoptCanvas(fCanvasMapCF);
  mfR2->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasMapCF->Divide(kNumberOfBoards,1,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    TVirtualPad *tmp = fCanvasMapCF->cd(i+1);
    ModelPads(i);
    fDiagrams[i]->Draw("COL");
  }

  TGCompositeFrame *mfR3 = new TGCompositeFrame(mf, 170, 20, kVerticalFrame);
  mf->AddFrame(mfR3, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
  CreateHistory(mfR3);
}
//====================
void DataMonitor::CreateScans(TGCompositeFrame *mf) {
  TGCompositeFrame *mf2 = new TGCompositeFrame(mf, 170, 20, kVerticalFrame);
  mf->AddFrame(mf2, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
  TRootEmbeddedCanvas *embeddedCanvas2 = new TRootEmbeddedCanvas("overall scan",mf2,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId2 = embeddedCanvas2->GetCanvasWindowId();
  fCanvasScans = new TCanvas("CanvasScans", 10, 10, cId2);
  embeddedCanvas2->AdoptCanvas(fCanvasScans);
  mf2->AddFrame(embeddedCanvas2, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasScans->Divide(kNumberOfBoards,1,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    TVirtualPad *tmp = fCanvasScans->cd(i+1);
    fScanAdc[i]->Draw("COL");
  }
}
//====================
void DataMonitor::CreateScansWav(TGCompositeFrame *mf) {
  TGCompositeFrame *mf2 = new TGCompositeFrame(mf, 170, 20, kVerticalFrame);
  mf->AddFrame(mf2, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
  TRootEmbeddedCanvas *embeddedCanvas2 = new TRootEmbeddedCanvas("average waveform",mf2,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId2 = embeddedCanvas2->GetCanvasWindowId();
  fCanvasScansWav = new TCanvas("CanvasScansWav", 10, 10, cId2);
  embeddedCanvas2->AdoptCanvas(fCanvasScansWav);
  mf2->AddFrame(embeddedCanvas2, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasScansWav->Divide(kNumberOfBoards,1,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    TVirtualPad *tmp = fCanvasScansWav->cd(i+1);
    fScanWav[i]->Draw("hist");
  }
}
//====================
void DataMonitor::CreateScansPed(TGCompositeFrame *mf) {
  TGCompositeFrame *mf2 = new TGCompositeFrame(mf, 170, 20, kVerticalFrame);
  mf->AddFrame(mf2, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
  TRootEmbeddedCanvas *embeddedCanvas2 = new TRootEmbeddedCanvas("overall scan pedestal",mf2,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId2 = embeddedCanvas2->GetCanvasWindowId();
  fCanvasScansPed = new TCanvas("CanvasScansPed", 10, 10, cId2);
  embeddedCanvas2->AdoptCanvas(fCanvasScansPed);
  mf2->AddFrame(embeddedCanvas2, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasScansPed->Divide(kNumberOfBoards,1,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    TVirtualPad *tmp = fCanvasScansPed->cd(i+1);
    fScanPed[i]->Draw("COL");
  }
}
//====================
void DataMonitor::CreateScansSgn(TGCompositeFrame *mf) {
  TGCompositeFrame *mf2 = new TGCompositeFrame(mf, 170, 20, kVerticalFrame);
  mf->AddFrame(mf2, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
  TRootEmbeddedCanvas *embeddedCanvas2 = new TRootEmbeddedCanvas("overall signals",mf2,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId2 = embeddedCanvas2->GetCanvasWindowId();
  fCanvasScansSgn = new TCanvas("CanvasScansSgn", 10, 10, cId2);
  embeddedCanvas2->AdoptCanvas(fCanvasScansSgn);
  mf2->AddFrame(embeddedCanvas2, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasScansSgn->Divide(kNumberOfBoards,1,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    TVirtualPad *tmp = fCanvasScansSgn->cd(i+1);
    fScanSgn[i]->Draw("hist");
  }
}
//====================
void DataMonitor::CreateLearning(TGCompositeFrame *mf) {
  TGCompositeFrame *mfR3 = new TGCompositeFrame(mf, 170, 20, kVerticalFrame);
  mf->AddFrame(mfR3, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
  TRootEmbeddedCanvas *embeddedCanvas2 = new TRootEmbeddedCanvas("learning",mfR3,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId2 = embeddedCanvas2->GetCanvasWindowId();
  fCanvasLearning = new TCanvas("CanvasLearning", 10, 10, cId2);
  embeddedCanvas2->AdoptCanvas(fCanvasLearning);
  mfR3->AddFrame(embeddedCanvas2, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasScans->Divide(kNumberOfBoards,1,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    TVirtualPad *tmp = fCanvasLearning->cd(i+1);
    fHitLearning[i]->Draw();
  }
}
//====================
void DataMonitor::CreateSignals(TGCompositeFrame *mf) {
  TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("signalsplot",mf,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId = embeddedCanvas->GetCanvasWindowId();
  fCanvasMapS = new TCanvas("CanvasMapS", 10, 10, cId);
  embeddedCanvas->AdoptCanvas(fCanvasMapS);
  mf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasMapS->Divide(kNumberOfChannels,kNumberOfBoards,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    for(int j=0; j!=kNumberOfChannels; ++j) {
      TVirtualPad *tmp = fCanvasMapS->cd(i*kNumberOfChannels+j+1);
      tmp->SetTopMargin(0.1);
      tmp->SetBottomMargin(0.1);
      tmp->SetLeftMargin(0.2);
      tmp->SetRightMargin(0.1);
      fSignal[i][j]->Draw("HIST");
    }
  }
  fCanvasMapS->SetEditable(kFALSE);
}
//====================
void DataMonitor::CreateHeights(TGCompositeFrame *mf) {
  TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("heightsplot",mf,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId = embeddedCanvas->GetCanvasWindowId();
  fCanvasMapH = new TCanvas("CanvasMapH", 10, 10, cId);
  embeddedCanvas->AdoptCanvas(fCanvasMapH);
  mf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasMapH->Divide(kNumberOfChannels,kNumberOfBoards,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    for(int j=0; j!=kNumberOfChannels; ++j) {
      TVirtualPad *tmp = fCanvasMapH->cd(i*kNumberOfChannels+j+1);
      tmp->SetTopMargin(0.1);
      tmp->SetBottomMargin(0.1);
      tmp->SetLeftMargin(0.2);
      tmp->SetRightMargin(0.1);
      fHeight[i][j]->Draw("HIST");
    }
  }
  fCanvasMapH->SetEditable(kFALSE);
}
//====================
void DataMonitor::CreateWidths(TGCompositeFrame *mf) {
  TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("heightsplot",mf,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId = embeddedCanvas->GetCanvasWindowId();
  fCanvasMapW = new TCanvas("CanvasMapW", 10, 10, cId);
  embeddedCanvas->AdoptCanvas(fCanvasMapW);
  mf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasMapW->Divide(kNumberOfChannels,kNumberOfBoards,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    for(int j=0; j!=kNumberOfChannels; ++j) {
      TVirtualPad *tmp = fCanvasMapW->cd(i*kNumberOfChannels+j+1);
      tmp->SetTopMargin(0.1);
      tmp->SetBottomMargin(0.1);
      tmp->SetLeftMargin(0.2);
      tmp->SetRightMargin(0.1);
      fWidth[i][j]->Draw("HIST");
    }
  }
  fCanvasMapW->SetEditable(kFALSE);
}
//====================
void DataMonitor::CreateHitSummary(TGCompositeFrame *mf) {
  TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("hitsummaryplot",mf,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId = embeddedCanvas->GetCanvasWindowId();
  fCanvasMapHS = new TCanvas("CanvasMapHS", 10, 10, cId);
  embeddedCanvas->AdoptCanvas(fCanvasMapHS);
  mf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasMapHS->Divide(2,kNumberOfBoards/2,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    TVirtualPad *tmp = fCanvasMapHS->cd(i+1);
    tmp->SetTopMargin(0.1);
    tmp->SetBottomMargin(0.1);
    tmp->SetLeftMargin(0.2);
    tmp->SetRightMargin(0.1);
    //fHitSummary[i]->Draw("HIST");
    fHitSummary[i]->Draw("E");
  }
  fCanvasMapHS->SetEditable(kFALSE);
}
//====================
void DataMonitor::CreateTimeSummary(TGCompositeFrame *mf) {
  TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("hitsummaryplot",mf,0.95*__window_wd__, 0.9*__window_ht__,kSunkenFrame);
  Int_t cId = embeddedCanvas->GetCanvasWindowId();
  fCanvasMapTS = new TCanvas("CanvasMapTS", 10, 10, cId);
  embeddedCanvas->AdoptCanvas(fCanvasMapTS);
  mf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasMapTS->Divide(2,kNumberOfBoards/2,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    TVirtualPad *tmp = fCanvasMapTS->cd(i+1);
    tmp->SetTopMargin(0.1);
    tmp->SetBottomMargin(0.1);
    tmp->SetLeftMargin(0.2);
    tmp->SetRightMargin(0.1);
    //fTimeSummary[i]->Draw("HIST");
    fTimeSummary[i]->Draw("E");
  }
  fCanvasMapTS->SetEditable(kFALSE);
}
//====================
void DataMonitor::CreateHistory(TGCompositeFrame *mf) {
  TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas("historygraph",mf, 0.95*__window_wd__, 0.2*__window_ht__, kSunkenFrame);
  Int_t cId = embeddedCanvas->GetCanvasWindowId();
  fCanvasMapHI = new TCanvas("CanvasMapHI", 10, 10, cId);
  embeddedCanvas->AdoptCanvas(fCanvasMapHI);
  mf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,2,2));
  fCanvasMapHI->Divide(kNumberOfBoards,0,0,0);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    TVirtualPad *tmp = fCanvasMapHI->cd(i+1);
    tmp->SetTopMargin(0.1);
    tmp->SetBottomMargin(0.1);
    tmp->SetLeftMargin(0.2);
    tmp->SetRightMargin(0.1);
    fHistory[i]->Draw("APL");
    fHistory[i]->SetMarkerStyle(20);
    fHistory[i]->SetMarkerColor(kBlue-3);
    fHistory[i]->SetTitle( Form("Charge collected vs time B%d",i) );
  }
  fCanvasMapHI->SetEditable(kFALSE);
}
//====================
void DataMonitor::RefreshAll() {
  //std::cout << "  >> RefreshAll called" << std::endl;

  // refresh TabContainerOne
  if(fTabContainerOne) {
    if(fTabContainerOne->GetCurrent()==0&&fCanvasScans) {
      fCanvasScans->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfBoards; ++i) {
	if(fNotInstalled[i]) continue;
	fCanvasScans->cd(i+1)->Modified();
	fCanvasScans->cd(i+1)->Update();
      }
      fCanvasScans->Modified();
      fCanvasScans->Update();
      fCanvasScans->SetEditable(kFALSE);
    }
    if(fTabContainerOne->GetCurrent()==1&&fCanvasScansSgn) {
      fCanvasScansWav->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfBoards; ++i) {
	if(fNotInstalled[i]) continue;
	fCanvasScansWav->cd(i+1)->Modified();
	fCanvasScansWav->cd(i+1)->Update();
      }
      fCanvasScansWav->Modified();
      fCanvasScansWav->Update();
      fCanvasScansWav->SetEditable(kFALSE);
    }
    if(fTabContainerOne->GetCurrent()==2&&fCanvasScansPed) {
      fCanvasScansPed->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfBoards; ++i) {
	if(fNotInstalled[i]) continue;
	fCanvasScansPed->cd(i+1)->Modified();
	fCanvasScansPed->cd(i+1)->Update();
      }
      fCanvasScansPed->Modified();
      fCanvasScansPed->Update();
      fCanvasScansPed->SetEditable(kFALSE);
    }
    if(fTabContainerOne->GetCurrent()==3&&fCanvasScansSgn) {
      fCanvasScansSgn->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfBoards; ++i) {
	if(fNotInstalled[i]) continue;
	fCanvasScansSgn->cd(i+1)->Modified();
	fCanvasScansSgn->cd(i+1)->Update();
      }
      fCanvasScansSgn->Modified();
      fCanvasScansSgn->Update();
      fCanvasScansSgn->SetEditable(kFALSE);
    }
  }

  if(fTabContainer) {
    fCanvasMapHI->SetEditable(kTRUE);
    for(int i=0; i!=kNumberOfBoards; ++i) {
      if(fNotInstalled[i]) continue;
      fCanvasMapHI->cd(i+1)->Modified();
      fCanvasMapHI->cd(i+1)->Update();
    }
    fCanvasMapHI->Modified();
    fCanvasMapHI->Update();
    fCanvasMapHI->SetEditable(kFALSE);
  }
  if(fTabContainer) {
    if(fTabContainer->GetCurrent()==1&&fCanvasMap) {
      //std::cout << "   current canvas 1" << std::endl;
      fCanvasMap->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfChannels*kNumberOfBoards; ++i) {
	if(fNotInstalled[i/kNumberOfChannels]) continue;
	fCanvasMap->cd(i+1)->Modified();
	fCanvasMap->cd(i+1)->Update();
      }
      fCanvasMap->Modified();
      fCanvasMap->Update();
      fCanvasMap->SetEditable(kFALSE);
    }
    if(fTabContainer->GetCurrent()==2&&fCanvasMapS) {
      //std::cout << "   current canvas 2" << std::endl;
      fCanvasMapS->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfChannels*kNumberOfBoards; ++i) {
	if(fNotInstalled[i/kNumberOfChannels]) continue;
	fCanvasMapS->cd(i+1)->Modified();
	fCanvasMapS->cd(i+1)->Update();
      }
      fCanvasMapS->Modified();
      fCanvasMapS->Update();
      fCanvasMapS->SetEditable(kFALSE);
    }
    if(fTabContainer->GetCurrent()==3&&fCanvasMapH) {
      //std::cout << "   current canvas 3" << std::endl;
      fCanvasMapH->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfChannels*kNumberOfBoards; ++i) {
	if(fNotInstalled[i/kNumberOfChannels]) continue;
	fCanvasMapH->cd(i+1)->Modified();
	fCanvasMapH->cd(i+1)->Update();
      }
      fCanvasMapH->Modified();
      fCanvasMapH->Update();
      fCanvasMapH->SetEditable(kFALSE);
    }
    if(fTabContainer->GetCurrent()==4&&fCanvasMapW) {
      //std::cout << "   current canvas 4" << std::endl;
      fCanvasMapW->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfChannels*kNumberOfBoards; ++i) {
	if(fNotInstalled[i/kNumberOfChannels]) continue;
	fCanvasMapW->cd(i+1)->Modified();
	fCanvasMapW->cd(i+1)->Update();
      }
      fCanvasMapW->Modified();
      fCanvasMapW->Update();
      fCanvasMapW->SetEditable(kFALSE);
    }
  }
  if(fTabContainer2) {
    if(fTabContainer2->GetCurrent()==1&&fCanvasMapCL) {
      fCanvasMapCL->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfBoards*4; ++i) {
	if(fNotInstalled[i]) continue;
	fCanvasMapCL->cd(i+1)->Modified();
	fCanvasMapCL->cd(i+1)->Update();
      }
      fCanvasMapCL->Modified();
      fCanvasMapCL->Update();
      fCanvasMapCL->SetEditable(kFALSE);
    }
    if(fTabContainer2->GetCurrent()==0&&fCanvasMapHS) {
      fCanvasMapHS->SetEditable(kTRUE);
      for(int i=0; i!=kNumberOfBoards; ++i) {
	if(fNotInstalled[i]) continue;
	fCanvasMapHS->cd(i+1)->Modified();
	fCanvasMapHS->cd(i+1)->Update();
      }
      fCanvasMapHS->Modified();
      fCanvasMapHS->Update();
      fCanvasMapHS->SetEditable(kFALSE);
    }
    
    fCanvasMapTS->SetEditable(kTRUE);
    for(int i=0; i!=kNumberOfBoards; ++i) {
      if(fNotInstalled[i]) continue;
      fCanvasMapTS->cd(i+1)->Modified();
      fCanvasMapTS->cd(i+1)->Update();
    }
    fCanvasMapTS->Modified();
    fCanvasMapTS->Update();
    fCanvasMapTS->SetEditable(kFALSE);
    
    fCanvasMapCF->SetEditable(kTRUE);
    for(int i=0; i!=kNumberOfBoards; ++i) {
      if(fNotInstalled[i]) continue;
      fCanvasMapCF->cd(i+1)->Modified();
      fCanvasMapCF->cd(i+1)->Update();
    }
    fCanvasMapCF->Modified();
    fCanvasMapCF->Update();
    fCanvasMapCF->SetEditable(kFALSE);
  }
  //std::cout << "  << RefreshAll called" << std::endl;
}
//====================
void DataMonitor::StyleH1(TH1D *tmp) {
  tmp->SetFillColor(kBlue-3);
  tmp->SetLineColor(kBlue-3);
  tmp->SetStats(1);
  tmp->GetYaxis()->SetLabelSize(0.1);
  tmp->GetXaxis()->SetLabelSize(0.1);
  tmp->GetXaxis()->SetNdivisions(504);
  tmp->GetYaxis()->SetNdivisions(508);
}
//====================
void DataMonitor::ReCreateHistograms() {
  // creates new histogram objects
  // reads kNumberOfBoards, kNumberOfChannels, kNumberOfSamples
  std::cout << "ReCreateHistograms INIT" << std::endl;
  //====== HISTOGRAMS
  fTimeSummaryTemp = new TH1D( "TST","TST",kNumberOfSamples,-0.5,kNumberOfSamples-0.5);
  for(int i=0; i!=kNumberOfBoards; ++i) {
    fDiagrams[i] = new TH2Poly( Form("BD%d",i), Form("[%d]%s  %s;X [mm];Y [mm]",i,
						     fBoardTech[i].Data(),fBoardCELL[i].Data()),-8,+8, -8,+8);
    fDiagrams[i]->SetStats(0);
    int st = -kNumberOfChannels/2;
    fHistory[i] = new TGraph();
    fHitSummary[i] = new TProfile( Form("HS%d",i),Form("HS%d;chn  idx;counts",i),kNumberOfChannels,
				   st-0.5,st+kNumberOfChannels-0.5,"S"); 
    fClusters_Num[i] = new TH1D( Form("CLUN%d",i),Form("CLUN%d",i),5,-0.5,4.5);     StyleH1(fClusters_Num[i]);
    fClusters_Wid[i] = new TH1D( Form("CLUW%d",i),Form("CLUW%d",i),30, -0.5, 29.5); StyleH1(fClusters_Wid[i]);
    fClusters_Amp[i] = new TH1D( Form("CLUA%d",i),Form("CLUA%d",i),100, 0, 10000);  StyleH1(fClusters_Amp[i]);
    fClusters_xCE[i] = new TH1D( Form("CLUX%d",i),Form("CLUX%d",i),100,-5.1,+5.1);  StyleH1(fClusters_xCE[i]);
    fTimeSummary[i] = new TProfile( Form("TS%d",i),Form("TS%d;sample  idx;counts",i),
				    kNumberOfSamples,-0.5,kNumberOfSamples-0.5,"S");
    fScan[i] = new TH2D( Form("SCAN%d",i), Form("SCAN%d;Channel;Sample",i),
			 kTotalNumberOfChannels, -0.5, kTotalNumberOfChannels-0.5,
			 kNumberOfSamples, -0.5, kNumberOfSamples-0.5 );
    fScanAdc[i] = new TProfile2D( Form("SCANADC%d",i), Form("SCANADC%d;Channel;Sample;< ADC >",i),
				  kTotalNumberOfChannels, -0.5, kTotalNumberOfChannels-0.5,
				  kNumberOfSamples, -0.5, kNumberOfSamples-0.5 );
    fScanWav[i] = new TProfile( Form("SCANWAV%d",i), Form("SCANWAV%d;Sample;< ADC >",i),
				kNumberOfSamples, -0.5, kNumberOfSamples-0.5 );
    fScanPed[i] = new TProfile2D( Form("SCANPED%d",i), Form("SCANPED%d;Channel;Sample; < ADC-PED >",i),
				  kTotalNumberOfChannels, -0.5, kTotalNumberOfChannels-0.5,
				  kNumberOfSamples, -0.5, kNumberOfSamples-0.5 );
    fScanSgn[i] = new TProfile( Form("SCANSGN%d",i), Form("SCANSGN%d;Channel;< sgn >",i),
				kTotalNumberOfChannels, -0.5, kTotalNumberOfChannels-0.5 );
    for(int j=0; j!=kNumberOfChannels; ++j) {
      fChannel[i][j] = NULL;
      fSignal[i][j] = NULL;
      fHeight[i][j] = NULL;
      fWidth[i][j] = NULL;
      fChannel[i][j] = new TH1D( Form("C%dP%d",i,j),Form("C%dP%d",i,j),
				 kNumberOfSamples,-0.5,kNumberOfSamples-0.5);
      fSignal[i][j] = new TProfile( Form("S%dP%d",i,j),Form("S%dP%d",i,j),
				    kNumberOfSamples,-0.5,kNumberOfSamples-0.5);
      fHeight[i][j] = new TH1D( Form("H%dP%d",i,j),Form("H%dP%d",i,j),100,0,800);
      fWidth[i][j] = new TH1D( Form("W%dP%d",i,j),Form("W%dP%d",i,j),100,0,1000);
      StyleH1(fChannel[i][j]);
      StyleH1(fSignal[i][j]);
      StyleH1(fHeight[i][j]);
      StyleH1(fWidth[i][j]);
    }
  }
  std::cout << "ReCreateHistograms DONE" << std::endl;
}
//====================
DataMonitor::DataMonitor(TApplication *app, UInt_t w, UInt_t h) {
  std::cout << "DM: I am awaking" << std::endl;
  fLearning = kTRUE;
  fReady = kFALSE;
  fApp = app;
  gStyle->SetTitleFontSize(0.1);

  // INITIALIZATION
  fNoEventsSampled = 0;
  for(int bd=0; bd!=kNumberOfBoards; ++bd) {
    fDiagrams[bd] = NULL;
    fIntMean[bd] = 8;
    fIntHWindow[bd] = 5;
  }
  fThisRun = NULL;
  fEventsSampled = NULL;
  fThisEvent = -1;
  fEventsReaded = 0;
  fEventsCataloged = 0;

  //====== MAPCOLLECTION
  fMapCollection = new MappingTableCollection();
  std::cout << " * loading map collection for... ";
  {
    const char *mfiles[] = {"./Run_Data/B00e.root",
			    "./Run_Data/B01a.root",
			    "./Run_Data/Z03k.A.root",
			    "./Run_Data/Z03k.B.root",
			    "./Run_Data/Z03k.D.root",
			    "./Run_Data/V00a.root"};
    for(unsigned i=0; i<sizeof(mfiles)/sizeof(mfiles[0]); i++) {
      TFile *ff = new TFile(mfiles[i]);
      MappingTable *mtable = (MappingTable*)ff->Get("MappingTable"); assert(mtable);
      std::cout << mtable->GetTag().Data() << " ";
      fMapCollection->AddMappingTable(mtable);
    }
  }
  std::cout << std::endl;
  fMapCollection->CreateLookupTables();

  ReadPosition();
  ReadConfig();
  ConfigureChannels();
  ReCreateHistograms();


  ///////////////
  // WINDOW ONE
  fWindowOne = new TGMainFrame(gClient->GetRoot(), __window_wd__, __window_ht__);
  fWindowOne->SetWindowName("Quick Scan");
  fWindowOne->SetWMPosition(0,0);
  fTabContainerOne = new TGTab(fWindowOne,96,26);
  fWindowOne->AddFrame(fTabContainerOne,
		       new TGLayoutHints(kLHintsTop | kLHintsExpandX,2,2,2,2));
  TGCompositeFrame *tab_1_0 = fTabContainerOne->AddTab("Last Event");
  CreateScans(tab_1_0);
  TGCompositeFrame *tab_1_1 = fTabContainerOne->AddTab("Average Waveform");
  CreateScansWav(tab_1_1);
  TGCompositeFrame *tab_1_2 = fTabContainerOne->AddTab("Last Event Reduced");
  CreateScansPed(tab_1_2);
  TGCompositeFrame *tab_1_3 = fTabContainerOne->AddTab("Average Signals");
  CreateScansSgn(tab_1_3);
  fWindowOne->MapSubwindows();
  fWindowOne->Layout();
  fWindowOne->MapWindow();
  fWindowOne->Move(0,0);
  ///////////////


  
  /*
  //====== WINDOWS
  fWindowHitSummary = new TGMainFrame(gClient->GetRoot(), __window_wd__, __window_ht__);//950, 1050);
  fWindowHitSummary->SetWindowName("On-the-fly Cluster Reconstruction");
  fWindowHitSummary->SetWMPosition(2000,0);
  fTabContainer2 = new TGTab(fWindowHitSummary,96,26);
  fWindowHitSummary->AddFrame(fTabContainer2, new TGLayoutHints(kLHintsTop | kLHintsExpandX,2,2,2,2));
  TGCompositeFrame *tab0_2 = fTabContainer2->AddTab("Weighted Average of Hit Distribution");
  CreateHitSummary(tab0_2); ///
  TGCompositeFrame *tab1_2 = fTabContainer2->AddTab("Clusters");
  CreateClusters(tab1_2);   ///

  fWindowTimeSummary = new TGMainFrame(gClient->GetRoot(), __window_wd__, __window_ht__);// 950, 1050);
  fWindowTimeSummary->SetWindowName("Weighted Average of Wave Profile");
  fWindowTimeSummary->SetWMPosition(100,0);
  CreateTimeSummary(fWindowTimeSummary); ///

  //========= SIGNAL INSPECTOR
  fWindowDetails = new TGMainFrame(gClient->GetRoot(), __window_wd__, __window_ht__); //1900, 1000);
  fWindowDetails->SetWindowName("Basic Information");
  fWindowDetails->SetWMPosition(0,0);
  TGCompositeFrame *mfR1 = new TGCompositeFrame(fWindowDetails, 170, 20, kHorizontalFrame);
  fWindowDetails->AddFrame(mfR1, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,2,2));
  StampRun(mfR1); ///
  fTabContainer = new TGTab(fWindowDetails,96,26);
  fWindowDetails->AddFrame(fTabContainer, new TGLayoutHints(kLHintsTop | kLHintsExpandX,2,2,2,2));
  TGCompositeFrame *tab0 = fTabContainer->AddTab("Configuration");
  CreateDisplayConfiguration(tab0);
  TGCompositeFrame *tab1 = fTabContainer->AddTab("Last Event");
  CreateChannels(tab1);
  TGCompositeFrame *tab2 = fTabContainer->AddTab("Signal");
  CreateSignals(tab2);
  TGCompositeFrame *tab3 = fTabContainer->AddTab("Height");
  CreateHeights(tab3);
  TGCompositeFrame *tab4 = fTabContainer->AddTab("Integral");
  CreateWidths(tab4);
  //TGCompositeFrame *tab6 = fTabContainer->AddTab("History");
  //CreateHistory(tab6);

  
  fWindowDetails->MapSubwindows();
  fWindowDetails->Layout();
  fWindowDetails->MapWindow();
  fWindowDetails->Move(0,0);

  fWindowHitSummary->MapSubwindows();
  fWindowHitSummary->Layout();
  fWindowHitSummary->MapWindow();
  fWindowHitSummary->Move(3000,0);

  fWindowTimeSummary->MapSubwindows();
  fWindowTimeSummary->Layout();
  fWindowTimeSummary->MapWindow();
  fWindowTimeSummary->Move(1500,0);

  fTabContainer->SetTab(0);
  */

  fReady = kTRUE;
  RefreshAll();

  std::cout << " DM: Please wait few seconds. Starting engine..." << std::endl;
}
//====================
DataMonitor::~DataMonitor() {
  fWindowDetails->Cleanup();
  fWindowHitSummary->Cleanup();
  fWindowTimeSummary->Cleanup();
  fApp->Terminate();
}

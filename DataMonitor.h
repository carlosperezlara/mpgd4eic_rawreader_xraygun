#ifndef DATAMONITOR_H
#define DATAMONITOR_H

#include <TApplication.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGTab.h>
#include <TGTextEdit.h>
#include <TGNumberEntry.h>
#include <TGIcon.h>

#include <MappingTableCollection.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TTimer.h>
#include <TImage.h>

const Int_t kMaxNumberOfBoards=9;
const Int_t kMaxNumberOfChannels=128;

class DataMonitor {
 private:
  void CreateScansAdc(TGCompositeFrame *mf);
  void CreateScansPed(TGCompositeFrame *mf);
  void CreateScansSgn(TGCompositeFrame *mf);
  void CreateScansWav(TGCompositeFrame *mf);
  void Learn();
  void CreateAnalysisMax(TGCompositeFrame *mf);
  void CreateAnalysisSgn(TGCompositeFrame *mf);
  void CreateAnalysisWid(TGCompositeFrame *mf);
  void CreateAnalysisAmp(TGCompositeFrame *mf);
  void CreateAnalysisCen(TGCompositeFrame *mf);


  void CreateDisplayConfiguration(TGCompositeFrame *mf);
  void CreateClusters(TGCompositeFrame *mf);
  void CreateChannels(TGCompositeFrame *mf);
  void CreateSignals(TGCompositeFrame *mf);
  void CreateHeights(TGCompositeFrame *mf);
  void CreateWidths(TGCompositeFrame *mf);
  void CreateHitSummary(TGCompositeFrame *mf);
  void CreateTimeSummary(TGCompositeFrame *mf);
  void CreateHistory(TGCompositeFrame *mf);
  const MappingPlane* GetMappingPlane(Int_t bd, TString sp, Int_t pl=0);
  void ReadPosition();
  void ReadConfig();
  void ConfigureChannels();
  void StampRun(TGCompositeFrame *mf);
  void ModelPads(Int_t bd);
  TString GetLocalCell(Int_t bd);
  void StyleH1(TH1D*);
  void GetXYFromCell(TString, Double_t&,Double_t&);
  TString GetCellName(Double_t,Double_t);
  void ReCreateHistograms();
  TApplication *fApp;

  TProfile *fSignal[kMaxNumberOfBoards][kMaxNumberOfChannels]; // display kMaxNumberOfChannels strips
  TProfile *fHitSummary[kMaxNumberOfBoards];
  TProfile *fTimeSummary[kMaxNumberOfBoards];
  TH1D *fHeight[kMaxNumberOfBoards][kMaxNumberOfChannels]; // display kMaxNumberOfChannels strips
  TH1D *fWidth[kMaxNumberOfBoards][kMaxNumberOfChannels]; // display kMaxNumberOfChannels strips
  TH1D *fClusters_Num[kMaxNumberOfBoards];
  TH1D *fClusters_Wid[kMaxNumberOfBoards];
  TH1D *fClusters_Amp[kMaxNumberOfBoards];
  TH1D *fClusters_xCE[kMaxNumberOfBoards];
  TH1D *fTimeSummaryTemp;
  TH2Poly *fDiagrams[kMaxNumberOfBoards];


  TH1D *fHitLearning[kMaxNumberOfBoards];
  
  TH2D *fScan[kMaxNumberOfBoards];
  TH1D *fChannel[kMaxNumberOfBoards][kMaxNumberOfChannels];

  TGMainFrame *fWindowOne;
  TCanvas *fCanvasScansAdc;
  TCanvas *fCanvasScansWav;
  TCanvas *fCanvasScansPed;
  TCanvas *fCanvasScansSgn;
  TCanvas *fCanvasLearning;
  TProfile2D *fScanAdc[kMaxNumberOfBoards];
  TProfile   *fScanWav[kMaxNumberOfBoards];
  TProfile2D *fScanPed[kMaxNumberOfBoards];
  TProfile   *fScanSgn[kMaxNumberOfBoards];

  TGMainFrame *fWindowTwo;
  TCanvas *fCanvasAnalysisMax;
  TCanvas *fCanvasAnalysisSgn;
  TCanvas *fCanvasAnalysisWid;
  TCanvas *fCanvasAnalysisAmp;
  TCanvas *fCanvasAnalysisCen;
  TH2D *fAnalysisMax[kMaxNumberOfBoards];
  TH2D *fAnalysisSgn[kMaxNumberOfBoards];
  TH1D *fAnalysisAmp[kMaxNumberOfBoards];
  TH1D *fAnalysisWid[kMaxNumberOfBoards];
  TH1D *fAnalysisCen[kMaxNumberOfBoards];
  
  TCanvas *fCanvasMap;
  TCanvas *fCanvasMapS;
  TCanvas *fCanvasMapH;
  TCanvas *fCanvasMapW;
  TCanvas *fCanvasMapHS;
  TCanvas *fCanvasMapTS;
  TCanvas *fCanvasMapCF;
  TCanvas *fCanvasMapHI;
  TCanvas *fCanvasMapCL;
  TGraph *fHistory[kMaxNumberOfBoards];


  TGMainFrame *fWindowHitSummary;
  TGMainFrame *fWindowTimeSummary;
  TGMainFrame *fWindowDetails;
  TGTextEdit *fConfigurationFile;
  MappingTableCollection *fMapCollection;

  
  TTimer *fRefresh;
  Double_t fPedestals[kMaxNumberOfBoards][kMaxNumberOfChannels];
  Int_t fThisEvent;
  Int_t fNoEventsSampled;
  Double_t fEventsReaded;
  Double_t fEventsCataloged;
  TString fClosestCell;
  Double_t fPosX;
  Double_t fPosY;
  TString fBoardCode[kMaxNumberOfBoards];
  TString fBoardTech[kMaxNumberOfBoards];
  TString fBoardCELL[kMaxNumberOfBoards];
  Double_t fPitchX[kMaxNumberOfBoards][kMaxNumberOfChannels];
  Double_t fPeriodY[kMaxNumberOfBoards][kMaxNumberOfChannels];
  Double_t fStretch[kMaxNumberOfBoards][kMaxNumberOfChannels];
  Int_t fDREAMChannels[kMaxNumberOfBoards];
  Int_t fDREAMChannel[kMaxNumberOfBoards][kMaxNumberOfChannels];
  Bool_t fNotInstalled[kMaxNumberOfBoards];
  Double_t fIntMean[kMaxNumberOfBoards];
  Double_t fIntHWindow[kMaxNumberOfBoards];
  Bool_t fReady;
  Bool_t fLearning;
  Bool_t fEBEPedestals;
  Int_t fCellChannels[kMaxNumberOfBoards][20][2];
  Double_t fCellPitch[kMaxNumberOfBoards][20];
  
 public:
  DataMonitor(TApplication *app, UInt_t w, UInt_t h);
  virtual ~DataMonitor();
  void RefreshAll();
  TH1D* GetChannel(int b, int s) {return fChannel[b][s];}
  TH2D* GetScan(int b) {return fScan[b];}
  void NewRun(Int_t run);
  void NewEvent(Int_t evr);
  void Merge();
  void SetDREAMChannel(Int_t bd,Int_t ch, Int_t idx) {fDREAMChannel[bd][ch] = idx;}
  Int_t GetDREAMChannels(Int_t bd) {return fDREAMChannels[bd];}
  Int_t GetDREAMChannel(Int_t bd,Int_t ch) {return fDREAMChannel[bd][ch];}
  Int_t IsBoardNotInstalled(Int_t bd) {return fNotInstalled[bd];}
  void SetXY(Double_t x, Double_t y) {fPosX=x;fPosY=y;fClosestCell=GetCellName(x,y);}
  Bool_t IsReady() {return fReady;}
  Bool_t IsLearning() {return fLearning;}
  
  TGLabel *fThisRun;
  TGLabel *fEventsSampled;
  TGLabel *fSamplingFraction;

  TGTab *fTabContainerOne;
  TGTab *fTabContainerTwo;

  TGTab *fTabContainer;
  TGTab *fTabContainer2;

  // these values are not constant anymore
  // in fact they will be overwritten by config.dat
  // before reserving memory for the histogram
  // in order to reduce the load for some applications
  Int_t kNumberOfBoards=9;
  Int_t kTotalNumberOfChannels=128;
  Int_t kNumberOfChannels=26;
  Int_t kNumberOfSamples=16;
  Int_t kNumberOfEventsLearning=500;
  Int_t kMergeRefresh=1000; // noevents
  Int_t kNewEventRefresh=100; // noevents

  ClassDef(DataMonitor, 0)
};

#endif

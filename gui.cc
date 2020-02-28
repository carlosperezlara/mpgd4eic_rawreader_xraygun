#include <TApplication.h>
#include <pmonitor/pmonitor.h>
#include "DataMonitor.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>

DataMonitor *gDM;
TH1D *gHist1D;
TH2D *gHist2D;
Int_t gCurrentRun;

//======================
int pinit() {
  //pidentify(0);
  return 0;
}
//======================
int process_event (Event * e) {
  //std::cout << " [pmonitor::process_event called" << std::endl;
  if(!gDM->IsReady()) return 0;
  
  //REAL DATA
  Int_t run = e->getRunNumber();
  Int_t evnr = e->getEvtSequence();

  double xpos, ypos;
  Packet *p926 = e->getPacket(926);
  if(p926) {
    xpos = p926->iValue(0);
    ypos = p926->iValue(1);
    delete p926;
    xpos = -xpos/100;
    ypos = -ypos/100-50;
    gDM->SetXY(xpos,ypos);
    std::cout << " * positions " << xpos << "  " << ypos << std::endl;
  }

  if(gCurrentRun!=run) {
    gDM->NewRun(run);
    gCurrentRun = run;
  }
  
  Packet *p = e->getPacket(3000);
  if(p) {
    gDM->NewEvent( evnr );
    int nfeu = p->iValue(0, "NR_FEU");
    //std::cout << " Number of FEUs " << nfeu << std::endl;
    if(nfeu>9) {
      std::cout << " [pmonitor] Ev" << evnr;
      std::cout << " NR_FEU reports to be greater than 9. Unreasonable. Skipping event." << std::endl;
    } else {
      for(int feu=0; feu<nfeu; ++feu) {
	if(gDM->IsBoardNotInstalled(feu)) continue;
	Int_t feuid = p->iValue(feu, "FEU_ID");
	Int_t channels =  p->iValue(feuid, "CHANNELS");
	Int_t samples =  p->iValue(feuid, "SAMPLES");
	Int_t nchips = p->iValue(feuid, "NR_DREAM"); //max 2
	//std::cout << " Reading FEU ID " << feuid << " which has " << nchips*64 << " channels and " << samples << " samples." <<  std::endl;
	if(samples>256) {
	  std::cout << " [pmonitor] Ev" << evnr <<  " FEU ID " << feuid;
	  std::cout << " reports number of samples greater than 256. Skipping FEU." << std::endl;
	  continue;
	}
	if(nchips>8) {
	  std::cout << " [pmonitor] Ev" << evnr << " FEU ID " << feuid;
	  std::cout << " reports number of chips greater than 8. Skipping FEU." << std::endl;
	  continue;
	}

	// Performs a quick scan of the data being recorded
	// determines which channels are really exposed to beam
	// after learning it will know which cell to watch
	if(gDM->IsLearning()) {
	  for(int chip=0; chip<nchips; ++chip) {
	    if( !p->iValue(feuid,chip,"DREAM_ENABLED") ) continue;
	    for(int ch=0; ch!=64; ++ch) {
	      gHist2D = gDM->GetScan(feu);
	      for(int sa=0; sa!=samples; ++sa) {
		Int_t dreamch = ch+64*chip;
		Int_t adc = p->iValue( feuid, dreamch, sa);
		gHist2D->Fill( double(dreamch), double(sa),  double(adc) );
	      }	  
	    }
	  }
	}
	// Checks only channels of interest for
	// a quick analysis of the data
	if(!gDM->IsLearning()) {
	  int nch = gDM->GetDREAMChannels(feu);
	  int minch = gDM->GetDREAMChannel(feu,0);
	  int maxch = minch + nch;
	  //std::cout << " Listening to " << minch << "--" << maxch << std::endl;
	  for(int ch=minch; ch!=maxch; ++ch) {
	    Int_t chipD = ch/64;
	    if( !p->iValue(feuid,chipD,"DREAM_ENABLED") ) continue;
	    gHist1D = gDM->GetChannel(feu,ch-minch);
	    for(int sa=0; sa!=samples; ++sa) {
	      Int_t adc = p->iValue( feuid, ch, sa);
	      gHist1D->SetBinContent( sa+1, double(adc) );
	      //std::cout << adc << "|";
	    }
	  }
	}
	//std::cout << std::endl;
      }
    }
    delete p;
  }
  //std::cout << " pmonitor::process_event built ]. Launching merge..." << std::endl;
  gDM->Merge();
  //std::cout << " DONE" << std::endl << std::endl;
  return 0;
}
//======================
//======================
//======================
int main(int nn, char** arg) {
  if(nn>1) {
    //TString filename = "/data/purschke/mrwell/mrwell-00000363-0000.evt";
    TString filename = arg[1];
    std::cout << filename.Data() << std::endl;
    pfileopen( filename.Data() );
  } else {
    //rcdaqopen();
  }
  TApplication *app = new TApplication("gui",0,0);
  gDM = new DataMonitor(app, 1600, 850);
  
  pstart();
  app->Run();
  return 0;
}

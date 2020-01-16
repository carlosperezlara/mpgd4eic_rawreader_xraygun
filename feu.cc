
#include <iostream>
#include <arpa/inet.h>

#include <string.h>


#include "feu.h"

#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
int init_done = 0;

using namespace std;

TH1F *hmax; 
TH1F *hped; 
TH1F *hsignal; 
TH1F *hsignalg; 

TH1F *h3; 
TH1F *h4; 

TH2F *sed_100;
TH2F *sed_119;

TH1F *sed_proj_100;
TH1F *sed_proj_119;

TH1F *hcent;
TH2F *hcorr;
TH2F *hres_corr;


int pinit()
{

  if (init_done) return 1;
  init_done = 1;

  hmax = new TH1F ( "hmax","max signal distribution", 128, 0, 4200); 
  hped = new TH1F ( "hped","pedestal distribution", 128, 0, 500); 
  hsignal = new TH1F ( "hsignal","signal distribution per channel", 128, 0, 3000); 
  hsignalg = new TH1F ( "hsignalg","signal distribution per event", 256, 0, 10000); 

  //  h3 = new TH1F ( "h3","Signal distribution per event", 128, 0, 1000); 
  //h4 = new TH1F ( "h4","Signal distribution", 128, 0, 1000); 


  // not the single event displays
  sed_100 = new TH2F ( "sed_100","FEU 100 SED Run = 0000 Evt 000000", 2*64, -0.5, 2*64-0.5, 64, -0.5, 63.5);
  sed_100->SetStats(0);
  sed_100->GetXaxis()->SetTitle("channel");
  sed_100->GetYaxis()->SetTitle("sample");
  sed_100->GetYaxis()->SetTitleOffset(1.6);
  sed_100->GetYaxis()->SetTitleOffset(1.7);
  
  sed_119 = new TH2F ( "sed_119","FEU 119 SED Run = 0000 Evt 000000", 8*64, -0.5, 8*64-0.5, 31, -0.5, 31.5);  
  sed_119->SetStats(0);
  sed_119->GetXaxis()->SetTitle("channel");
  sed_119->GetYaxis()->SetTitle("sample");
  sed_119->GetYaxis()->SetTitleOffset(1.6);
  sed_119->GetYaxis()->SetTitleOffset(1.7);


  sed_proj_100 = new TH1F ( "sed_proj_100","FEU PROJ 100 SED Run = 0000 Evt 000000", 2*64, -0.5, 2*64-0.5);
  sed_proj_100->GetXaxis()->SetTitle("channel");
  sed_proj_100->GetYaxis()->SetTitle("ADC counts");

  hcent = new TH1F("hcent", "Centroid Distr.", 10000, -50, 50);
  hcent->GetXaxis()->SetTitle("Centroid (mm)");
  hcent->GetYaxis()->SetTitle("Counts");
  
  hcorr = new TH2F("hcorr", "Motor Position Vs Centoid", 1000, -50, 50, 10000, -50, 50);
  hcorr->GetXaxis()->SetTitle("Motor Position (mm)");
  hcorr->GetYaxis()->SetTitle("Centroid (mm)");
  
  hres_corr = new TH2F("hres_corr", "Motor Position Vs Centoid", 1000, -50, 50, 10000, -50, 50);
  hres_corr->GetXaxis()->SetTitle("Motor Position (mm)");
  hres_corr->GetYaxis()->SetTitle("Centroid-MP (mm)");
  
  return 0;

}
int reset_hists()
{

  // add something to preserve them as needed
  
  //  h2->Reset();
  sed_100->Reset();
  sed_proj_100->Reset();
  hcent->Reset();
}

int process_event (Event * e)
{
  static int old_runnumber = -1;

  static float  x_position = -10000;
  static float  y_position = -10000;


  
  if ( e->getEvtType() == BEGRUNEVENT || old_runnumber != e->getRunNumber() )
    {
      old_runnumber = e->getRunNumber();
      reset_hists();
      Packet *p926 = e->getPacket(926);
      if ( p926)
	{
	  x_position = p926->iValue(0);
	  y_position = p926->iValue(1);
	  delete p926;
	  cout << __FILE__ << " " << __LINE__ << " positions " << x_position << "  " << y_position << endl;
	}
      else
	{
	  x_position = -10000;
	  y_position = -10000;
	}
      
      hcent->Reset();
    }

  
  Packet *p = e->getPacket(3000);
  if (p)
    {
      
      sed_100->Reset();
      sed_119->Reset();
      sed_proj_100->Reset();
      
      
      char title[512];
      sprintf (title,"FEU 100 SED Run = %4d Evt %6d", e->getRunNumber(), e->getEvtSequence());
      sed_100->SetTitle(title);
      sprintf (title,"FEU 119 SED Run = %4d Evt %6d", e->getRunNumber(), e->getEvtSequence());
      sed_119->SetTitle(title);
      
      sprintf (title,"FEU PROJ 100 SED Run = %4d Evt %6d", e->getRunNumber(), e->getEvtSequence());
      sed_proj_100->SetTitle(title);
      
      
      
      int max_signal = 0;
      float global_signal = 0;
      
      for ( int ich = 0; ich < 2*64; ich++)
	{
	  // FEU id 100 first, pedestal ch 20..31
	  float pedestal = 0;
	  for ( int isample = 55; isample < 63;  isample++)
	    {
	      pedestal +=  p->iValue( 100, ich, isample);
	    }
	  
	  pedestal /= 9;
	  hped->Fill(pedestal);
	  float signal = 0;
	  
	  //      	  for ( int isample = 0; isample < p->iValue(100,"SAMPLES");  isample++)
      	  for ( int isample = 0; isample < 64;  isample++)
      	    {
	      int v =  pedestal - p->iValue( 100, ich, isample);
	      //int v = p->iValue( 100, ich, isample);
	      sed_100->Fill(ich, isample, v );
	      sed_proj_100->Fill(ich, v );
      	    }
	  
	  for ( int isample = 0; isample < 16;  isample++)
	    {
	      //	      int v =  pedestal - p->iValue( 100, ich, isample);
	      int v =  p->iValue( 100, ich, isample);
	      if (v > max_signal) max_signal = v; 
	      signal += v;
	    }
	  signal /= 16;
	  hsignal->Fill(signal);
	  
	  if ( signal > 10) global_signal += signal;
	}	  
      hmax->Fill(max_signal);
      hsignalg->Fill(global_signal);
      
      // for ( int ich = 0; ich < 8*64; ich++)
      // 	{
      
      // 	  // fill in the FEU 100 single event display
      // 	  for ( int isample = 0; isample < p->iValue(100,"SAMPLES");  isample++)
      // 	    {
      // 	      int v =  p->iValue( 100, ich, isample);
      // 	      sed_100->Fill(ich, isample, v );
      // 	    }
      
      //



      delete p;
    }


  
  for (int i=0; i<sed_proj_100->GetNbinsX(); i++)
    {
      sed_proj_100->SetBinError(i, 200); //100*.05*sed_proj_100->GetBinContent(ch_max));
    }
  
  sed_proj_100->Fit("pol1");
  TF1 *f1 = new TF1("f1", "pol1", 0, 1000);
  f1 = sed_proj_100->GetFunction("pol1");
  float yint = f1->GetParameter(0);
  float slope = f1->GetParameter(1);
  
  
  //float bin_max = sed_proj_100->GetMaximumBin();
  //float bin_min = sed_proj_100->GetMinimumBin();
  
  float peak_max = -100000000.0;
  float peak_max_raw = -100000000.0;
  float bin_amp=0.0;
  float ch_max = 0.0;
  int bin_max = 0;
  for (int i=2; i<sed_proj_100->GetNbinsX()-1; i++)
    {
      bin_amp = abs(sed_proj_100->GetBinContent(i) - ((slope*sed_proj_100->GetBinCenter(i))+yint));
      if (peak_max < bin_amp)
	{
	  peak_max = bin_amp;
	  ch_max = sed_proj_100->GetBinCenter(i);
	  bin_max = i;
	  peak_max_raw = sed_proj_100->GetBinContent(i) - ((slope*sed_proj_100->GetBinCenter(i))+yint);
	}
      
    }
  
  float baseline = (slope*ch_max)+yint;
  cout<<"Baseline: "<<baseline<<", Peak_Max: "<<peak_max<<", CH_Max: "<<ch_max<<endl;

  /*  
  TF1 *ff0 = new TF1("ff0", "pol0");
  ff0->FixParameter(0, baseline);
  TF1 *ff1 = new TF1("ff1", "gaus");
  ff1->SetParLimits(0, 0.5*peak_max_raw-baseline, peak_max_raw-baseline);
  ff1->SetParLimits(1, ch_max-3, ch_max+3);
  ff1->SetParLimits(2,0,5);
  
  TF1 *ff2 = new TF1("ff2", "ff0 + ff1");
  ff2->FixParameter(0, baseline);
  //ff2->SetParLimits(1, 0.75*peak_max_raw-baseline, peak_max_raw-baseline);
  ff2->SetParLimits(1, -200,-100);
  ff2->SetParLimits(2, ch_max-3, ch_max+3);
  ff2->SetParLimits(3,0.5,5);
  //ff2->SetParameter(0, (slope*ch_max)+yint);
  //ff2->FixParameter(1, ch_max);
  //sed_proj_100->GetXaxis()->SetRangeUser(ch_max-20, ch_max+20);
  sed_proj_100->Fit("ff2", "", "", ch_max-20, ch_max+20);
  //ff2 = sed_proj_100->GetFunction("gaus");
  float mean = ff2->GetParameter(1);

  cout<<"Centroid: "<<mean<<endl;  
  */

  float sum = 0.0;
  float cent = 0.0;
  for (int i=bin_max-2; i<bin_max+3; i++) //!! units=bins !!
    {
      bin_amp = abs(sed_proj_100->GetBinContent(i) - ((slope*sed_proj_100->GetBinCenter(i))+yint));
      cent += bin_amp*sed_proj_100->GetBinCenter(i);
      sum += bin_amp;
    }

  cent = cent/sum;
  
  cout<<"XmotorSteps: "<<x_position<<", YmotorSteps: "<<y_position<<endl;  
  cout<<"Centroid: "<<cent<<", XmotorPos: "<<x_position/160<<", YmotorPos: "<<(y_position/160)+50.<<endl;  


  cent = ((128 - cent) - 25)*1; //mm
  cent = 50.0-cent;
  float mp = (y_position/160)+64.0;

  cout<<"cent_scaled: "<<cent<<", mp_scaled: "<<mp<<endl;
  
  hcent->Fill(cent);
  hcorr->Fill(mp, cent, 1); //mm
  hres_corr->Fill(mp, cent-mp, 1); //mm
  
  return 0;
}


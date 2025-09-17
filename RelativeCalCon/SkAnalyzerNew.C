#define SkAnalyzer_cxx
#include "SkAnalyzerNew.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace std;

double CalCon(5e-3);
TString tName;
double nBin=500;
double minADC = 1000.;
double maxADC = 20000.;
double relGain[25],dummyVec[25],PeakGain[25];
const int CryRef=12;
TH1D *hSumADC,*hPkCls;
TH2D *hPkFrc, *hPkADC, *hTotal, *hCluster, *hScale, *hADC, *hCal, *hLaser;
// ==============================
void SkAnalyzer::Book()
{
  cout << " Booking histograms ... " << endl;
  
  hCluster = new TH2D("hCluster","Cluster Energy by Crystal",nCry,0.,(double)nCry,nBin,minADC,maxADC); th2Stack.push_back(hCluster);
  hPkFrc   = new TH2D("hPkFrc","Fraction Peak/Overall",nCry,0.,(double)nCry,101,0.,1.01); th2Stack.push_back( hPkFrc );
  hPkCls   = new TH1D("hPkCls","Spread Peak Counts",100,0.95,1.05); tObjStack.push_back( hPkCls );
  hPkADC   = new TH2D("hPkADC","ADC counts, peak",nCry,0.,(double)nCry,nBin,minADC,maxADC); th2Stack.push_back( hPkADC );
  hTotal   = new TH2D("hTotal","Total Energy by Crystal",nCry,0.,(double)nCry,nBin,minADC, maxADC); th2Stack.push_back(hTotal);
  hScale   = new TH2D("hScale","Scale Factors",5,-0.5,4.5,5,-0.5,4.5);th2Stack.push_back(hScale);
  hADC     = new TH2D("hADC","ADC counts by counter",nCry,0.,(double)nCry,nBin,minADC, maxADC); th2Stack.push_back( hADC );
  hLaser   = new TH2D("hLaser","ADC counts by counter, laser events",nCry,0.,(double)nCry,800,minADC,maxADC); th2Stack.push_back(hLaser);
  hCal     = new TH2D("hCal","ADC counts by counter,scaled",nCry,0.,(double)nCry,nBin,minADC, maxADC); th2Stack.push_back( hCal );
  hSumADC  = new TH1D("hSumADC","ADC sum",1000,0.,1e5); th1Stack.push_back( hSumADC );
}

// ==============================
void SkAnalyzer::iteration(double *gain, int *nEv, int *off , int *stp, TString *name){};

// ==============================
// questa funzione esengue nIter-volte:
// Loop( nEvt, iter+1, step,"") e Fit2Dmap(hPkADC,iter+1,PeakGain)

void SkAnalyzer::iterate( int nIter=2, int eBeam=-100, int step=10, TString OutputFileName="SaveResults.hist" )
{

  Book();
  SetFitStyle();
  for( int i=0; i<nCry;i++ ){
    relGain[i]=dummyVec[i]=1;
  }
  
  // int nEvt = -1*abs(eBeam);
  int nEvt = eBeam;
  for( int iter=0; iter<nIter; iter++ )
    {
      printf(" *** Fraction at  Iteration number %d \n",iter);
      Loop( nEvt, iter+1, step,"");
      
      Fit2Dmap(hPkADC,iter+1,PeakGain);
      
      for( int i=0; i<nCry; i++ )
	{
	  relGain[i] *= PeakGain[i]/PeakGain[CryRef];
	  if( relGain[i]<0.1 )
	    relGain[i]=1;
	  printf("  %5.3f", relGain[i] );
	  if( i%5==4 )
	    printf("\n");
	}
    }
  
  Loop( nEvt,0,1,"");
  Fit2Dmap(hPkADC,0,PeakGain);
  Fit2Dmap(hLaser,0,PeakGain);
  gStyle->SetStatX(0.45);
  Fit2Dmap(hCluster,-1,PeakGain);
  for( int iCry=0; iCry<nCry; iCry++ )
    {
      hPkCls->Fill( PeakGain[iCry]/PeakGain[CryRef] );      
      if( (iCry%5)== 0 ) cout << endl;
      double sf = ( relGain[iCry]==0 ? 1: relGain[iCry]/relGain[CryRef]);
      printf(" %5.3f,",sf);
      hScale->Fill( iCry%5,4-iCry/5, sf );
    }
  //DrawBig2Dmap(hADC,eBeam, "ADC/");
  //DrawBig2Dmap(hCal,eBeam, "Cal/");

  //Draw2Dmap(hADC,eBeam);
  //Draw2Dmap(hCal,eBeam);
  
  TCanvas tFit("tFit","Counts in the 9 Clusters",1024,1024);
  TH1D *hfit = hCluster->ProjectionY();
  double pval[5], perr[5];
  FitHisto( hfit, pval, perr );
  tFit.SaveAs(OutputDir+"ClusterCounts.png");
  
  SaveResults( OutputFileName );
}

// ==============================
void SkAnalyzer::Loop() { Loop(-1,"SaveResults.hist"); }

// ==============================
void SkAnalyzer::Loop( int nEvt=-1, TString OutputFileName="SaveResults.hist" ){ Loop( nEvt,100, 1, OutputFileName) ;}

// ==============================
void SkAnalyzer::Loop( int nEvt=-1, int offs=100, int step=0, TString OutputFileName="SaveResults.hist" )
{
  cout << " ***  Analyze " << nEvt/step << "  data at " << eBeam << " beam energy; iter# " << offs << endl;
 
  eBeam = ( nEvt<0 ? (double) abs(nEvt) : 100 );
  Long64_t nbytes = 0, nb = 0;  
  Long64_t nentries = fChain->GetEntriesFast();  
  if( nEvt < 0 || nEvt>nentries ) nEvt = nentries;
   
  cout << " ***  Analyze " << nEvt/step << "  data with " << eBeam << " beam energy; iter# " << offs << endl;
  for( int i=0; i<nCry; i++ )
    {
      printf("  %5.3f", relGain[i] );
      if( i%5==4 ) printf("\n");
    }
  
  hLaser->Reset();
  hCal->Reset();
  hADC->Reset();
  hPkADC->Reset();
  hPkFrc->Reset();
  hTotal->Reset();
  hCluster->Reset();
  hSumADC->Reset();
    
  if (fChain == 0) return;
  for( Long64_t jentry=0; jentry<nEvt;jentry+=1 )
    {      
      Long64_t ientry = LoadTree(jentry);
      if( ientry < 0) break;
      if( ientry%10000 == offs ) printf(" +++   at event %5d \n", (int) ientry ) ;
      if( ientry%step != offs ) continue;
      
      nb = fChain->GetEntry(jentry);
      nbytes += nb;      

      isLaser = TagLaser();
      
      double sumADC=0;
      for( int i=0; i<nCry; i++)
	sumADC+=FitCry[i];
      hSumADC->Fill( sumADC );

      // per i run di laser fa il fit sui singoli cristalli ed esce dal loop
      if( isLaser )
	for( int i=0; i<nCry; i++ )
	  hLaser->Fill(i,FitCry[i]);

      if( isLaser )
	continue;

      //per i run di fascio fa il fit dei cristalli e calcola/calibra l'energia
      bool good=true;
      for( int i=0 ; i<nCry; i++ )
	{	  
	  if( FitTim[i]<0 || FitTim[i]>127 ) {
	    good = false;
	    break;
	  } 
	  eCry[i] = FitCry[i]/relGain[i];
	  hADC->Fill(i,FitCry[i]) ;
	  hCal->Fill(i,eCry[i]) ;
	}
      if( !good )
	continue;
	  
      double ePeak(0.);           
      int iPeak = GetPeakCry( ePeak ) ;
      if( iPeak<0 || iPeak>nCry )
	continue ;
      
      //      if( offs==2 ) cout << "step 0" << endl;
      
      double xPeak(0.), yPeak(0.);
      double eCls = GetClusterEnergy( iPeak, xPeak, yPeak ) ;
      double eTot = GetTotalEnergy() ;
      if( ePeak/eTot>0.95 )
	continue;

      //      if( offs==2 ) cout << "step " << iPeak << " " << ePeak << "  " << eTot << endl;
      
      hTotal->Fill( iPeak, eTot );
      hCluster->Fill( iPeak, eCls );      
      hPkFrc->Fill( (double) iPeak, ePeak/eTot );
      if( ePeak/eTot>0.7 && ePeak/eTot<0.95 )
	hPkADC->Fill( (double) iPeak, ePeak );
    }

}

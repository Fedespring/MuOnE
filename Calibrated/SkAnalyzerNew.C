#define SkAnalyzer_cxx
#include "../SkAnalyzerNew.h"
#include <TH2.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace std;

double relGain[25];
double CorFactorLowEnergy[25] ={ 1.013, 1.048, 1.038, 0.998, 1.030,
				 1.010, 0.891, 1.152, 1.072, 1.042,
				 1.172, 1.217, 1.226, 1.376, 1.038,
				 1.036, 1.042, 1.018, 1.123, 1.008,
				 1.013, 1.136, 0.997, 1.087, 1.024 };

double CorFactor50GeV[25] ={ 0.996, 0.990, 1.001, 0.996, 0.992,
			     0.998, 1.001, 1.000, 0.998, 0.999,
			     1.003, 0.999, 1.000, 1.001, 0.999,
			     1.001, 1.002, 1.001, 1.001, 1.000,
			     0.997, 0.995, 0.999, 0.999, 0.996};

double CorFactor75GeV[25] ={ 0.986, 0.963, 0.989, 0.975, 0.980,
			     0.987, 0.990, 0.986, 0.979, 0.975,
			     0.982, 0.981, 0.979, 0.985, 0.989,
			     0.993, 0.995, 0.993, 0.990, 0.989,
			     0.987, 0.975, 0.987, 0.985, 0.985};

double CorFactor150GeV[25] ={ 1.001, 1.016, 0.999, 1.009, 1.004,
			      1.002, 1.002, 1.006, 1.007, 1.004,
			      1.000, 1.006, 1.004, 1.002, 1.002,
			      1.002, 1.002, 1.000, 1.001, 1.003,
			      1.004, 1.008, 1.004, 1.005, 1.004} ;

double CorFactor2025[25] ={  1.000, 1.000, 1.000, 1.000, 1.000,
			     1.000, 1.021, 1.019, 1.021, 1.000,
			     1.000, 1.016, 0.990, 0.996, 1.000,
			     1.000, 1.014, 1.015, 1.025, 1.000,
			     1.000, 1.000, 1.000, 1.000, 1.000};
			     
double calcon =5.460e-3; // GeV/ADC count PEAK

TH2D *ClusEnergy,*ClusFrac,*ClusTot;
TH2D *hBeam, *hLaser, *hPkEnf, *hPkClf, *hRMSMeV, *hRMSADC;
TH2D *hX, *hY, *hXY;
TNtuple *tnt;
TFile *fOut;

// ====================================================================================
void SkAnalyzer::iterate( int nIter, int eBeam, int step, TString name ){};
// ====================================================================================
void SkAnalyzer::Book()
{
  double nBinL=1000;
  double nBin=250;
  double minAdc = eBeam*0.75 ;
  double maxAdc = eBeam*1.10 ;
  double maxFrc =1.15;
  double minFrc =0.75;
  
  if( eBeam<20 ) maxAdc = eBeam*1.25;
  if( eBeam==1  ) { nBin=200; minAdc=0.25 ; maxAdc=10.25; maxFrc=3.; }
  else if( eBeam ==2 ) { nBin=200; minAdc=0.5 ; maxAdc=10.5; maxFrc=3.; }
  else if( eBeam == 5 || eBeam == 6 ) { nBin=200; minAdc=2.5 ; maxAdc=22.5; maxFrc=3.; }

  //  fOut = new TFile( "Results.root", "recreate" );  
  hX = new TH2D("hX","E frac vs X shower",150,-7.5,7.5,100,0.9,1.1); tObjStack.push_back(hX);
  hY = new TH2D("hY","E frac vs X shower",150,-7.5,7.5,100,0.9,1.1); tObjStack.push_back(hY);
  hXY = new TH2D("hXY","Y vs X shower",75,-7.5,7.5,75,-7.5,7.5); tObjStack.push_back(hXY);
  tnt = new TNtuple("tnt","tnt","f_Crys:f_xPeak:f_yPeak:f_ePeak:f_eClus"); tObjStack.push_back(tnt);
  
  hBeam = new TH2D("hBeam","Beam Signal ADC",nCry,0,nCry,nBinL,0.,25000); tObjStack.push_back(hBeam);
  hLaser = new TH2D("hLaser","Laser Signal",nCry,0,nCry,nBinL,500.,10500); tObjStack.push_back(hLaser);
  hPkEnf = new TH2D("hPkEnf","Peak ADC Energy/E beam ",nCry,0.,nCry,nBin,0.5,maxFrc); tObjStack.push_back(hPkEnf);
  hPkClf = new TH2D("hPkClf","Peak ADC Energy/E cluster ",nCry,0.,nCry,nBin,0.91,1.01); tObjStack.push_back(hPkClf);
  hRMSMeV = new TH2D("hRMSMeV","Channel Noise Beam (MeV)",nCry,0,nCry,100,0.,100.); tObjStack.push_back(hRMSMeV);
  hRMSADC = new TH2D("hRMSADC","Channel Noise Laser (ADC)",nCry,0,nCry,nBin,0.,100.); tObjStack.push_back(hRMSADC);
  ClusEnergy = new TH2D("ClusEnergy","Energy 3x3 Cluster in Calo",nCry,0.,nCry,nBin,minAdc,maxAdc); tObjStack.push_back(ClusEnergy);  
  ClusFrac = new TH2D("ClusFrac","Energy 3x3 Cluster/E beam",nCry,0.,nCry,nBin,minFrc,maxFrc); tObjStack.push_back(ClusFrac);
  ClusTot = new TH2D("ClusTot","Emore ../Relnergy 3x3 Cluster/E Calo",nCry,0.,nCry,nBin,0.55,1.05); tObjStack.push_back(ClusTot);

  
}
// ====================================================================================
void SkAnalyzer::iteration (double *gain, int *nEv, int *off, int *stp, TString *name){}
void SkAnalyzer::Loop() { Loop(-1,"SaveResults.hist"); }
void SkAnalyzer::Loop( int nEvt=-1, TString OutputFileName="SaveResults.hist" ){ Loop( nEvt,2024, 1, OutputFileName) ;}
void SkAnalyzer::Loop( int nEvt=-1, int year=2024, int step=0, TString OutputFileName="SaveResults.hist" )
{
  eBeam = ( nEvt<0 ? (double) abs(nEvt) : 100 );
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntriesFast();  
  if( nEvt < 0 || nEvt>nentries ) nEvt = nentries;

  fOut = new TFile( OutputFileName.Data(), "recreate" );
  Book();  
  calcon = CalCon( year, relGain );

  cout << " ***  Analyze " << nEvt << " data in year" << year << ", E beam= "<< eBeam  << " Cal Constant = " << calcon << endl;

  for( int i=0; i<nCry; i++ )
    { printf("%5.4f ",relGain[i]); if( i%5==4 ) cout <<endl;}
  
  if (fChain == 0) return;

  double dummyVec[25];
  for( int i=0; i<25;i++ ) {dummyVec[i]=1; }

  if( year==2023 )
    {
      if( eBeam == 75 ) { for( int i=0; i<nCry; i++ ) relGain[i] /= CorFactor75GeV[i]; }
      else if( eBeam < 10 ) { for( int i=0; i<nCry; i++ ) relGain[i] /= CorFactorLowEnergy[i]; }
      else if( eBeam == 50 ) { for( int i=0; i<nCry; i++ ) relGain[i] /= CorFactor50GeV[i]; }
      else if( eBeam ==150 ) { for( int i=0; i<nCry; i++ ) relGain[i] /= CorFactor150GeV[i]; }
    }
  else if( year==2025 ) { for( int i=0; i<nCry; i++ ) relGain[i] /= CorFactor2025[i]; }
  
  for( Long64_t jentry=0; jentry<nEvt;jentry+=1 )
    {
      
      Long64_t ientry = LoadTree(jentry);
      if( ientry < 0) break;
      if( ientry%100000 == 0 ) printf(" +++   at event %5d \n", (int) ientry ) ;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      isLaser = TagLaser();
      if( isLaser == 0 )
	{
	  for( int i=0 ; i<nCry; i++ )
	    {
	      hBeam->Fill((double) i,FitCry[i]);
	      eCry[i] = calcon*FitCry[i]/relGain[i];
	      hRMSMeV->Fill((double) i,RmsCry[i]*relGain[i]*calcon*1000);
	    }
	  
	  double ePeak(0.);           
	  int iPeak = GetPeakCry( ePeak ) ;
	  if( iPeak<0 || iPeak>nCry ) continue ;
    
	  double xPeak(0.), yPeak(0.);
	  double eCls = GetClusterEnergy( iPeak, xPeak, yPeak ) ;
	  double eTot = GetTotalEnergy() ;

	  hPkEnf->Fill( double(iPeak), ePeak/eBeam ) ;
	  ClusEnergy->Fill( double(iPeak), eCls ) ;
	  ClusFrac->Fill( double(iPeak), eCls/eBeam ) ;
	  if( eCls>0 ) hPkClf->Fill( double(iPeak), ePeak/eCls ) ;
          if( eTot>0 ) ClusTot->Fill( double(iPeak), eCls/eTot ) ;
	  if( eCls/eBeam>0.75 ) hXY->Fill( xPeak, yPeak );
	  hX->Fill( xPeak,eCls/eBeam );
	  hY->Fill( yPeak,eCls/eBeam );

	  if( eCls/eBeam>1.5) continue;
	  if( eCls/eBeam<0.5) continue;
	  if( OnBorder[iPeak] ) continue;
	  tnt->Fill( (float) iPeak, (float) xPeak, (float) yPeak, (float) ePeak, (float) eCls ); 
	}
      else
	{
	  for( int i=0 ; i<nCry; i++ )
	    {
	      hLaser->Fill((double) i,FitCry[i]);
	      hRMSADC->Fill((double) i,RmsCry[i]);
	    }      
	}
    }


  Plot2D( ClusEnergy, (int) eBeam  );
  Plot2D( hPkEnf, (int) eBeam  );
  Plot2D( hRMSMeV, (int) eBeam  );
  Plot2D( hRMSADC, (int) eBeam  );
  Plot2D( hLaser, (int) eBeam  );
  Plot2D( hBeam, (int) eBeam  );

  SetFitStyle();
  Fit2Dmap( ClusEnergy,-eBeam,dummyVec ) ;
  Fit2Dmap( hPkEnf,eBeam,dummyVec ) ;
  Fit2Dmap( ClusFrac,-eBeam,dummyVec ) ;

  double fpar[5],ferr[5];
  TString tName;
  tName.Form("%s_Sum.jpg",ClusEnergy->GetName());
  TCanvas tCls("tCls","Total Ecal energy",2058,1024);
  FitHisto( ClusEnergy->ProjectionY(), fpar, ferr );
  tCls.SaveAs( tName.Data() );

  cout << " Histos " << tObjStack.size() << endl; 
  cout << " End of Loop, Save results to " << OutputFileName.Data() << endl ;

  for( uint i=0; i<tObjStack.size(); i++ ) tObjStack.at(i)->Write();
  
  // SaveResults(OutputFileName);
}

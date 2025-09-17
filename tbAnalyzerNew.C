#define tbAnalyzer_cxx
#include "tbAnalyzerNew.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStopwatch.h>

//TFile *fOut;
TTree *Skim;
TH1D *hBeam1D, *hLaser1D;
TH2D *hBeam, *hLaser;
int eBeam;
const double Alpha=1.75;
const double Beta=1.25;

/*  ========================================================================
    analyze a chain using just a fraction (=1/step) of the events, 
    with an arbitrary offset (offs)
    ==========================================================================*/
void tbAnalyzer::Loop() { Loop(-1,"SaveResults.hist"); }
void tbAnalyzer::Loop( int nEvt=-1, TString OutputFileName="SaveResults.hist" ){ Loop( nEvt, 0, 1, OutputFileName) ;}
void tbAnalyzer::Loop( int nEvt=-1, int offs=0, int step=1, TString OutputFileName="SaveResults.hist" )
{

  if (fChain == 0) return;
  
  // fOut = new TFile(OutputFileName.Data(),"recreate");

  TFile fOut(OutputFileName.Data(),"recreate");
  
  bool isLaser(0);
  float MaxCry[25], MaxTim[25], PedCry[25], RmsCry[25];
  float FitCry[25], FitTim[25], FitWlf[25], FitWrg[25];
  float PedMin[25], FitPed[25];
  
 
  Skim = new TTree("Skim","Skimmed Ntuple") ;
  Skim->Branch("isLaser",&isLaser,"isLaser/B");
  Skim->Branch("FitCry",FitCry,"FitCry[25]/F");
  Skim->Branch("FitTim",FitTim,"FitTim[25]/F");
  Skim->Branch("FitWlf",FitWlf,"FitWlf[25]/F");
  Skim->Branch("FitWrg",FitWrg,"FitWrg[25]/F");
  Skim->Branch("MaxCry",MaxCry,"MaxCry[25]/F");
  Skim->Branch("MaxTim",MaxTim,"MaxTim[25]/F");
  Skim->Branch("PedCry",PedCry,"PedCry[25]/F");
  Skim->Branch("RmsCry",RmsCry,"RmsCry[25]/F");
  Skim->Branch("PedMin",PedMin,"PedMin[25]/F");
  Skim->Branch("FitPed",FitPed,"FitPed[25]/F");
  Skim->Print();
  cout <<  "  ***** Test New Setup ******  " << endl;
  
  hBeam = new TH2D("hBeam","Sum ADC counts beam events",128,0.,128.,nCry,0.,(double) nCry);
  hLaser = new TH2D("hLaser","Sum ADC counts laser events",128,0.,128.,nCry,0.,(double) nCry); 

  // load gain factors
  double relGain[nCry] ;
  std::fill( relGain,relGain+nCry,1. ) ;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int nFiles = fChain->GetTreeNumber();
  int current_file(-1);
  int npr(0), nLaser(0);
  
  if( nEvt<0 ) { eBeam = abs(nEvt); nEvt=nentries ;} 
  else nEvt+=offs;

  cout << " ========================================================== " << endl;
  cout << " Start analysis on  " << nEvt-offs << " events from " << nFiles <<  "files " << endl;
  cout << " fraction analized is 1/" << step << ", offset is "<< offs    << endl;
  cout << " ========================================================== " << endl;

  // first fast loop, find peak delay
  for ( Long64_t jentry=1; jentry<5000 ;jentry+=1 )
    {
      Long64_t ientry = LoadTree(jentry);
      if( ientry < 0 ) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      FillEventHistogram( ientry );
      if( TagLaser() ) { hLaser->Add( hCrySgn ) ;}
      else             { hBeam->Add( hCrySgn ) ;}
    }

  hLaser1D = (TH1D*) hLaser->ProjectionX("hLaser1D");
  hBeam1D = (TH1D*) hBeam->ProjectionX("hBeam1D");

  double TsLaser = hLaser1D->GetMaximumBin()-1 ;
  double TsBeam  = hBeam1D->GetMaximumBin()-1 ;

  gStyle->SetOptStat(10);
  gStyle->SetOptFit(1);

  // Fit Laser Events
  TString tcName;
  const int npar=5;
  double fpar[npar]={hLaser1D->GetMaximum()-hLaser1D->GetBinContent(5), TsLaser,Alpha,Beta,hLaser1D->GetBinContent(5)};
  double ferr[npar]={100, TsLaser*0.1,0.,0.,100.};

  cout << " Start parameters:" ;
  for( int i=0; i<npar; i++ ) cout << "  " << fpar[i] ;
  cout << endl;
  
  TCanvas *tc1 = new TCanvas("tc1","Laser",100,100,512,256);
  //  FitHisto( fpar, ferr, hLaser1D, true );
  FitPedes( fpar, ferr, hLaser1D, true );
  tcName.Form("PlotLaser_%dGeV.png",eBeam);
  tc1->SaveAs(tcName.Data());

  cout << " Laser ";
  for( int j=0; j<5; j++) cout << fpar[j] << " ";
  cout << endl;
  
  // Fit data events; shape fixed to Laser result
  fpar[1] = TsBeam;
  ferr[2] = 0.;
  ferr[3] = 0.;
  fpar[4] = hBeam1D->GetBinContent(5);
  fpar[0] = hBeam1D->GetMaximum()-fpar[0];
  
  // ferr[2] = fpar[2]*0.1;   ferr[3] = fpar[3]*0.1;
  TCanvas *tc2 = new TCanvas("tc2","Beam Fix",100,300,512,256);
  //  FitHisto( fpar, ferr, hBeam1D, true );
  FitPedes( fpar, ferr, hBeam1D, true );
  tcName.Form("PlotBeam_%dGeV.png",eBeam);
  tc2->SaveAs(tcName.Data());

  cout << " Beam ";
  for( int j=0; j<5; j++) cout << fpar[j] << " ";
  cout << endl;

  TStopwatch timer;
  timer.Start();
  
  TCanvas *tc25 = new TCanvas("tc25","Crystal Gain",1600,900);
  tc25->Divide(5,5,0,0);
  
  TH1D *hFit = new TH1D();

  double tReal(0), tCPU(0);
  // ready for event analysis
  for ( Long64_t jentry=offs; jentry<nEvt;jentry+=1 ) {

    Long64_t ientry = LoadTree(jentry);
    if( ientry < 0 ) break;

    
    npr++ ;
    bool DisplayFit = false;
    if( npr%10000 == 1 )
      {
	tReal += timer.RealTime();
	tCPU += timer.CpuTime();
	timer.Stop();
	cout << " at event " << npr << " entry " << jentry << " @file:" << fChain->GetCurrentFile()->GetName() << " time " << tReal << "  " << tCPU << endl ;
	timer.Print();
	timer.Start();
	DisplayFit = true;
      }
	    
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Load Pedestals, then fill event histogram  
    FillEventHistogram( ientry );
    isLaser = (TagLaser())? 1: 0 ;

    ComputePedestalMin(1,30,PedMin);
    
    for( int i=0; i<nADC; i++ )
      {
	int iCry = CrystalMap[i];
	if( iCry==-1 ) continue;
	PedCry[iCry] = pedestal->at(iCry) ;
	RmsCry[iCry] = pedRMS->at(iCry) ;

	// fit signal
	fpar[4]=0;
	ferr[4]=1;
	
	if( isLaser )
	  {
	    fpar[1] = TsLaser;
	    ferr[1] = fpar[1]*0.1;
	    fpar[2] = Alpha;
	    fpar[3] = Beta;
	    ferr[2] = 0.;
	    ferr[3] = 0.;
	    fpar[4] = PedCry[iCry];
	    ferr[4] = RmsCry[iCry];     
	  }
	else
	  {
	    fpar[1] = TsBeam;
	    ferr[1] = fpar[1]*0.1;
	    fpar[2] = Alpha;
	    fpar[3] = Beta;
	    ferr[2] = 0.;
	    ferr[3] = 0.;
	    fpar[4] = PedCry[iCry];
	    ferr[4] = RmsCry[iCry];     

	  }
	//	FitCry[iCry] = SigCry( fpar, ferr, iCry, iCry);
        hFit = (TH1D*) hCrySgn->ProjectionX("hFit",iCry+1,iCry+1);
	fpar[0] = hFit->GetMaximum();
	ferr[0] = fpar[0]*0.1;
	
	if( DisplayFit )
	  {
	    tc25->cd(iCry+1);
	    FitPedes( fpar, ferr, (TH1D*)hFit->Clone(Form("Fit_%d",iCry)), true );
	    // cout << iCry << "   ";
	    // for( int j=0; j<5; j++) cout << fpar[j] << " ";
	    // cout << endl;
	  }
	else
	  { FitPedes( fpar, ferr, hFit , false );}
	
	FitCry[iCry] = fpar[0];
	FitTim[iCry] = fpar[1];
	FitWlf[iCry] = fpar[2];
	FitWrg[iCry] = fpar[3];
	FitPed[iCry] = fpar[4];
      }
    if( DisplayFit )
      {
	TString tcName;
	tcName.Form("Display_%d.C",(int)ientry);
	tc25->SaveAs(tcName);
      }
    Skim->Fill();
  }

  hLaser->Write();
  hLaser1D->Write();
  hBeam1D->Write();
  hBeam->Write();
  Skim->Write();
  Skim->Print();

  cout << " Close Output File " << fOut.GetName() << endl;
  
  fOut.Close();
  // delete fOut;
}

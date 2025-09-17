//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 25 19:07:35 2023 by ROOT version 6.24/06
// from TTree Skim/Skimmed Ntuple
// found on file: ../data/skim_230531/Chain100Center.root
//////////////////////////////////////////////////////////

#ifndef SkAnalyzer_h
#define SkAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.

#include "vector"
#include "iostream"
#include "cmath"
#include "math.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TString.h>
#include <TF1.h>
#include <TSystem.h>

class SkAnalyzer {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  const static int nCry=25;
  // Declaration of leaf types
  int isLaser;
  Float_t         FitCry[nCry];
  Float_t         FitTim[nCry];
  Float_t         FitWlf[nCry];
  Float_t         FitWrg[nCry];
  Float_t         MaxCry[nCry];
  Float_t         MaxTim[nCry];
  Float_t         PedCry[nCry];
  Float_t         RmsCry[nCry];
  //  Float_t         SigCry[nCry];

  Float_t        xCluster;   //
  Float_t        yCluster;   //
  Float_t        TrackXecal;   //
  Float_t        TrackYecal;   //
  
  // List of branches
  TBranch        *b_isLaser;   //!
  TBranch        *b_FitCry;   //!
  TBranch        *b_FitTim;   //!
  TBranch        *b_FitWlf;   //!
  TBranch        *b_FitWrg;   //!
  TBranch        *b_MaxCry;   //!
  TBranch        *b_MaxTim;   //!
  TBranch        *b_PedCry;   //!
  TBranch        *b_RmsCry;   //!

  // new branches

  TBranch        *b_xCluster;   //
  TBranch        *b_yCluster;   //
  TBranch        *b_TrackXecal;   //
  TBranch        *b_TrackYecal;   //

  
  
  // Added features
  bool OnBorder[nCry]= { 1,1,1,1,1 ,1,0,0,0,1 ,1,0,0,0,1 ,1,0,0,0,1, 1,1,1,1,1 };
  double alpha[nCry ]= { 1,1,1,1,1 ,1,1,1,1,1 ,1,1,1,1,1 ,1,1,1,1,1, 1,1,1,1,1 }; 
  double eCry[nCry];
  double CrySize=2.95;
  double xCenter[5]={-5.9,-2.95,0.0, 2.95, 5.9};
  double yCenter[5]={ 5.9, 2.95,0.0,-2.95,-5.9};
   
  SkAnalyzer(TTree *tree=0, TString TreeName="Skim" );
  SkAnalyzer(TString FileName, TString TreeName="Skim" );
  virtual ~SkAnalyzer();
  
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual void     Loop( int N, TString T );
  virtual void     Loop( int nevt, int ofs, int step, TString T );
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  virtual void     Book();
  virtual void     iterate( int nIter, int eBeam, int step, TString name );
  virtual void     iteration(double *gain, int *nEv, int *off , int *stp, TString *name);

  TString OutputDir = "/eos/user/f/federica/www/MuOne/test/";
  
  vector <TObject*> tObjStack;
  vector <TH1D*> th1Stack;
  vector <TH2D*> th2Stack;

  // ====================================================
  // Laser event if all the ACDs have at least 1000 counts
  bool TagLaser()
  {
    // double minADC=1000;
    //    for( int i=0; i<nCry; i++ ) if( FitCry[i]<minADC ) return false;
    double sum(0);
    for( int i=0; i<nCry; i++ ) sum+=FitCry[i]; 
    if( sum>19000 ) return true;
    return false;
  }
  // ====================================================
  double CalCon( int year, double *ptr )
  {
    double calcon(1);
    if( year==2023 ) // 2023 hE run not adjusted for laser
      {
	double sf[nCry]={0.995, 1.053, 1.035, 1.038, 1.012,
			 0.979, 0.951, 1.021, 0.982, 1.024,
			 1.009, 1.016, 1.000, 0.992, 1.066,
			 1.011, 1.019, 0.953, 0.989, 1.016,
			 0.978, 1.049, 1.020, 0.996, 1.000};
	calcon = 5.460e-3;
	for( int i=0; i<nCry; i++) *ptr++ = sf[i];
      }
    else if( year==2024 || year==20240 ) // 2024 July run
      {
	double sf[nCry] ={0.855, 1.000, 0.789, 0.735, 0.959,
			  0.988, 0.636, 1.715, 0.852, 1.000,
			  0.842, 0.893, 1.328, 1.015, 0.929,
			  0.979, 1.000, 0.837, 1.000, 0.744,
			  1.000, 0.785, 0.809, 0.954, 0.936};
	calcon = 9.05e-3;
	for( int i=0; i<nCry; i++) *ptr++ = sf[i];
      }
    else if( year==202506 ) // 2025 June run
      {
	double sf[nCry]={1.014, 0.977, 1.015, 0.988, 1.002,
			 1.016, 0.987, 0.998, 1.008, 1.024,
			 1.032, 0.978, 1.000, 0.988, 1.026,
			 1.025, 0.951, 1.008, 0.993, 0.958,
			 0.995, 1.018, 0.963, 1.042, 0.988};
	calcon =13.11e-3;
	for( int i=0; i<nCry; i++) *ptr++ = sf[i];				
      }
    else if( year==202507 ) // 2025 July run
      {
	double sf[nCry]={0.972, 1.000, 0.976, 0.979, 0.983,
			 0.965, 0.964, 0.941, 1.001, 1.007,
			 0.986, 1.019, 1.000, 1.002, 0.982,
			 0.969, 0.977, 0.954, 0.972, 0.978,
			 0.985, 0.993, 0.993, 0.983, 0.989};
	calcon =11.53e-3;
	for( int i=0; i<nCry; i++) *ptr++ = sf[i];	
      }
    else if( year==20230 ) // 2025 July run
      {
	double sf[nCry]={1.000, 0.950, 0.964, 0.966, 0.989,
			 1.015, 1.045, 0.976, 1.013, 0.977,
			 0.987, 0.980, 1.000, 1.004, 0.938,
			 0.987, 0.978, 1.043, 1.011, 0.984,
			 1.023, 0.948, 0.973, 1.002, 0.999};
	calcon = 5.46e-3/0.988;
	for( int i=0; i<nCry; i++) *ptr++ = sf[i];	
      }
    
    else
      {
	calcon = 1.;
	for( int i=0; i<nCry; i++) *ptr++ = 1.;
      }    
    return calcon;
  }
  // ====================================================
  double TimeCorrection( double x )
  {
    const int np=5;

    double t = x - (int) x;
    t -= 0.5 ;
    double p[np]={1.000,0.0116,-0.3065,-0.1174,0.3667};
    double val(0);
    for( int i=0; i<np; i++ ) val += p[i]*pow(t,i);
    return val;
  }
  // ======================================================
  inline int GetPeakCry( double &Epeak ) {
    int maxCry(-1);
    Epeak = 0. ;
    for( int iCry=0; iCry<nCry; iCry++ ) {
      if( eCry[iCry]< Epeak ) continue;
      maxCry=iCry;
      Epeak = eCry[iCry] ;
    }
    return maxCry ;
  }
  // ======================================================
  inline double GetTotalEnergy(){
    double s(0); 
    for( int iCry=0; iCry<nCry; iCry++ ) s+=eCry[iCry];
    return s;
  }
  // ======================================================
  inline double GetClusterNoise( int iPeak )
  {
    if( OnBorder[iPeak] ) return -1;
    double s(0);
    int xPeak = iPeak%5 ; 
    int yPeak = iPeak/5 ;
    for( int ix=xPeak-1; ix<=xPeak+1; ix++ ){
      for( int iy=yPeak-1; iy<=yPeak+1; iy++ ){ s+=pow(RmsCry[ix+5*iy],2) ;}}
    return sqrt(s);
  }
  // ======================================================
  inline double GetIntegratedClusterNoise( int iPeak )
  {
    if( OnBorder[iPeak] ) return -1;
    double s(0);
    int xPeak = iPeak%5 ; 
    int yPeak = iPeak/5 ;
    for( int ix=xPeak-1; ix<=xPeak+1; ix++ ){
      for( int iy=yPeak-1; iy<=yPeak+1; iy++ )
	{
	  double CryNoise = MaxTim[ix+5*iy]*pow(RmsCry[ix+5*iy],2) ; 
	  s+=CryNoise;
	}
    }
    return sqrt(s);
  }
  // ======================================================
  inline void Plot2D( TH2D *h, int eBeam )
  {
    TCanvas tc("tc",h->GetTitle(),1600,900);
    h->Draw("box");
    TString tcName;
    tcName.Form("%d_GeV_2D%s.png",eBeam,h->GetName());
    tc.SaveAs( OutputDir+tcName.Data() );
  }
  // ======================================================
  inline void Plot2Dmap( TH2D *h, bool zoom=false )
  {
    auto TC25 = (TCanvas*) gROOT->FindObject("TC25");
    if (TC25) delete TC25;
    TC25 = new TCanvas("TC25",h->GetTitle(),2048,1024);
    if( zoom )     TC25->Divide(3,3,0,0) ;
    else TC25->Divide(5,5,0,0) ;

    for( int ipad=1; ipad<nCry+1; ipad++ )
      {
	if( zoom && OnBorder[ipad-1] ) continue;
	TC25->cd( ipad );
	gPad->SetTickx(2);
	gPad->SetTicky(2);
	TH1D* h1 = h->ProjectionY("h1",ipad,ipad);
	h1->SetLineColor( kRed+2 );
	h1->SetFillColor( kRed+2 );
	h1->DrawCopy();
      }
    TString tcName;
    tcName.Form("%s.png",h->GetTitle());
    cout << "Save Canvas as " << tcName.Data() << endl;
    TC25->SaveAs(OutputDir+tcName.Data());
  }
  // ======================================================
  inline bool IsInCluster( int iCrys, int iPeak )
  {
    if( OnBorder[iPeak] ) return false;
    int xCrys = iCrys%5;
    int yCrys = iCrys/5;
    int xPeak = iPeak%5;
    int yPeak = iPeak/5;

    if( abs(xCrys-xPeak)>1 ) return false;
    if( abs(yCrys-yPeak)>1 ) return false;

    return true ;
  }
  // ======================================================  
  inline double GetArrayEnergy( vector<int> *CryArray, double &xClus, double &yClus, double offs=4.9 )
  {
    xClus = yClus = 0.;
    int *ptr ;
    
    //for( uint i=0 ; i<CryArray->size(); i++ ) cout << CryArray->at(i) << " " ;
    // return 0;
    
    double esum(0);
    for( uint i=0 ; i<CryArray->size(); i++ ) esum += eCry[CryArray->at(i)] ;
    if( esum<=0 ) return esum;

    double wSum(0);
    for( uint a=0 ; a<CryArray->size(); a++ ) 
      {
	int i = CryArray->at(a);
	double x = xCenter[i%5];
	double y = yCenter[i/5];
	double e = eCry[i];
	double w =(e<=0? 0:max(0.,offs+log(e/esum))); 
	xClus += w*x ;
	yClus += w*y ;
	wSum += w;
      }
    if( wSum <= 0 )
      { xClus=yClus=999 ;}
    else
      { xClus/=wSum; yClus/=wSum;}

    return esum; 
  }
  // ======================================================  
  inline double GetClusterEnergy( int iPeak, double &xClus, double &yClus, double offs=4.9 )
  {
    if( OnBorder[iPeak] ) return -1;
    double s(0);
    double wCry[nCry];
    int xPeak = iPeak%5 ; 
    int yPeak = iPeak/5 ;
    for( int ix=xPeak-1; ix<=xPeak+1; ix++ ){
      for( int iy=yPeak-1; iy<=yPeak+1; iy++ ){ s+=eCry[ix+5*iy] ;}}

    xClus =0;
    yClus =0;
    double wsum=0;
    for( int ix=-1; ix<=+1; ix++ ){
      for( int iy=-1; iy<=+1; iy++ )
	{
	  int x = xPeak+ix;
	  int y = yPeak+iy;
	  double w = offs+log(eCry[x+5*y]/s) ;
	  if( w<0 ) w=0;
	  // xClus += w * (double) ix ;
	  // yClus += w * (double) iy ;
	  xClus += w * xCenter[x] ;
	  yClus += w * yCenter[y] ;
	  wsum += w;
	}
    }
    
    if( wsum<=0 ) { xClus = yClus = 999 ;}
    else { xClus /= wsum; yClus /= wsum; }

    return s;
   }
   // ======================================================
   inline double GetTotalEnergyMax(double *gn){
     double s(0); 
     for( int iCry=0; iCry<nCry; iCry++ ) s+=*gn*MaxCry[iCry];
     gn++;
     return s;
   }
  // ====================================================
  virtual void SaveResults( TString OutputFileName)
  {

	cout << " ====> Save results on " << OutputFileName.Data() << endl ;
	TFile fOut(OutputFileName.Data(),"RECREATE") ;
	// TFile fOut("SaveResults.hist","RECREATE") ;
	for( uint i=0; i<th1Stack.size(); i++ ) {th1Stack.at(i)->Write();} 
	for( uint i=0; i<th2Stack.size(); i++ ) {th2Stack.at(i)->Write();} 
	for( uint i=0; i<tObjStack.size(); i++ )
	  {
	    tObjStack.at(i)->GetName();
	    tObjStack.at(i)->Write();
	  }
	fOut.Print();
	fOut.Close();

  }

  // ===================================================================
  inline void SetFitStyle()
  {
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(10);
    //gStyle->SetStatX(0.45);
    gStyle->SetStatX(0.95);
    gStyle->SetStatY(0.94);
    gStyle->SetStatH(0.2);
    gStyle->SetStatW(0.2);     
    gStyle->SetStatFontSize(0.05);
  }
  // ===================================================================
  inline virtual void FitHisto( TH1D* hFit, double *fpar, double *ferr )
  {
    SetFitStyle() ;    
    hFit->SetMarkerColor(39);
    hFit->SetMarkerSize(0.5);
    hFit->SetMarkerStyle(20);
    hFit->SetLineColor(39);
    hFit->SetLineWidth(3);
    
    // define fit range & starting values with a gaussian fit
    TF1 GFit("GFit","gaus",hFit->GetBinLowEdge(1),hFit->GetBinLowEdge(hFit->GetNbinsX())+1);
    GFit.SetParameter(0,hFit->GetMaximum());
    GFit.SetParameter(1,hFit->GetBinLowEdge( hFit->GetMaximumBin()));
    GFit.SetParameter(2,hFit->GetRMS());

    double xmin = GFit.GetParameter(1)-1.5*GFit.GetParameter(2);
    double xmax = GFit.GetParameter(1)+2.5*GFit.GetParameter(2);

    hFit->Fit( &GFit,"N0","",xmin,xmax );

    // Cristal Ball
    TF1 cbFit("cbFit","crystalball",hFit->GetBinLowEdge(1),hFit->GetBinLowEdge(hFit->GetNbinsX())+1);
    cbFit.SetLineColor(46);
    cbFit.SetLineWidth(2);
    for( int i=0; i<3 ;i++)
      cbFit.SetParameter(i,GFit.GetParameter(i)) ;
    

    cbFit.SetParameter(3,-hFit->GetSkewness());
    cbFit.SetParameter(4,2);

    hFit->SetOption("EML");
    hFit->Fit(&cbFit,"","",xmin,xmax);

    // display Gauss and tail of the Xball
    for( int i=0; i<3 ;i++)
      GFit.SetParameter(i,cbFit.GetParameter(i)) ;
    
    auto xb = (TF1*) gROOT->FindObject("xb"); if( xb ) delete xb;
    auto gs = (TF1*) gROOT->FindObject("gs"); if( gs ) delete gs;
    
    xb = (TF1*) cbFit.Clone("xb");
    gs = (TF1*) GFit.Clone("gs");
    
    auto ft = (TF1*) gROOT->FindObject("ft");
    if( ft ) ft->Clear() ;
	    
    ft = new TF1("ft","xb-gs",hFit->GetBinLowEdge(1),hFit->GetBinLowEdge(hFit->GetNbinsX())+1);  

    GFit.SetLineStyle(9) ;
    GFit.SetLineColor(44);
    GFit.SetLineWidth(3);
    ft->SetLineStyle(9) ;
    ft->SetLineWidth(3) ;
    ft->SetLineColor(32);
    
    // hFit->SetMaximum(cbFit.GetParameter(0)*1.25);
    // hFit->GetXaxis()->SetRangeUser(cbFit.GetParameter(1)-15*cbFit.GetParameter(2) , cbFit.GetParameter(1)+3*cbFit.GetParameter(2));
    hFit->SetMaximum( cbFit.GetParameter(0)*1.25 );
    
    hFit->DrawCopy();
    GFit.Draw("same");
    ft->Draw("same");
      
    for( int i=0; i<5; i++ ) {
      fpar[i]=cbFit.GetParameter(i);
      ferr[i]=cbFit.GetParError(i);
    }
    
  }
  // ============================================================================
  virtual void DrawBig2Dmap( TH2D* h, double eBeam=1000, TString subDir = "./" )
  {
 
    for( int iCry=0 ; iCry<25; iCry++ )
      {
	auto TC25 = (TCanvas*) gROOT->FindObject("TC25");
	if (TC25) delete TC25;
	TC25 = new TCanvas("TC25",h->GetTitle(),2048,1024);
        
	TC25->cd();
	TH1D *hDraw = h->ProjectionY("hDraw",iCry+1,iCry+1);
	hDraw->SetTitle("");
	hDraw->SetLineWidth(2);
	hDraw->DrawCopy();
      
	TString s;
	s.Form("%s_%d%s",h->GetName(),(int)eBeam,".png");
	TC25->SaveAs(OutputDir+subDir+s.Data());
      }
  }
  // ============================================================================
  virtual void Draw2Dmap( TH2D* h, double eBeam=1000 )
  {
    auto TC25 = (TCanvas*) gROOT->FindObject("TC25");
    if (TC25) delete TC25;
    TC25 = new TCanvas("TC25",h->GetTitle(),2048,1024);
    TC25->Divide(5,5,0,0);

    for( int iCry=0 ; iCry<25; iCry++ )
      {
	TC25->cd( iCry+1 );
	TH1D *hDraw = h->ProjectionY("hDraw",iCry+1,iCry+1);
	hDraw->SetTitle("");
	hDraw->SetLineWidth(2);
	hDraw->DrawCopy();
      }
    TString s;
    s.Form("%s_%d%s",h->GetName(),(int)eBeam,".png");
    TC25->SaveAs(OutputDir+s.Data());
  }
    // ============================================================================
 
  virtual void Profile3D( TH3D* h, double eBeam=1000 ) 
  {
    auto TC25 = (TCanvas*) gROOT->FindObject("TC25");
    if (TC25) delete TC25;
    TC25 = new TCanvas("TC25",h->GetTitle(),2048,1024);
    TC25->Divide(5,5,0,0);

    TH2D *hDraw = (TH2D*) h->Project3D("zy");
    hDraw->Reset();
    cout << " axis " << hDraw->GetYaxis()->GetBinLowEdge(1) << endl;
    
    for( int iCry=0 ; iCry<25; iCry++ )
      {
	TC25->cd( iCry+1 );
	for( int j=1; j<=hDraw->GetNbinsX(); j++ )
	  {
	    for( int k=1; k<=hDraw->GetNbinsY(); k++ )
	      {
		hDraw->SetBinContent(j,k,h->GetBinContent( iCry+1,j,k ));
	      }
	  }
      
	hDraw->SetTitle("");
	hDraw->ProfileX()->DrawCopy();
	hDraw->Reset();
      }
    TString s;
    s.Form("%s_%d%s",h->GetName(),(int)eBeam,".png");
    TC25->SaveAs(OutputDir+s.Data());
    return;
    
    }
  
  // ============================================================================
  virtual void Fit2DLan( TH2D* h, double eBeam=1000 )
  {
    auto TC25 = (TCanvas*) gROOT->FindObject("TC25");
    if (TC25) delete TC25;
    TC25 = new TCanvas("TC25",h->GetTitle(),2048,1024);
    TC25->Divide(5,5,0,0);

    for( int iCry=0 ; iCry<25; iCry++ )
      {
	TC25->cd( iCry+1 );
	TH1D *hDraw = h->ProjectionY("hDraw",iCry+1,iCry+1);
	hDraw->SetTitle("");
	hDraw->SetLineWidth(2);
	if( hDraw->GetMean()>100 ) hDraw->Fit("landau","","",50,1000);
	hDraw->DrawCopy();
      }
    TString s;
    s.Form("%s_%d%s",h->GetName(),(int)eBeam,".png");
    TC25->SaveAs(OutputDir+s.Data());
  }
  // ============================================================================
  virtual void Fit2Dmap( TH2D* h , int offs, double *gn) 
  {
    const int npar=5;
    double fpar[npar],ferr[npar];
    const int MinCounts=100;
    
    //auto TC25 = (TCanvas*) gROOT->FindObject("TC25");
    //if (TC25) delete TC25;
    //TC25 = new TCanvas("TC25",h->GetTitle(),2048,1024);
    //if( offs<0 ) {TC25->Divide(3,3,0,0) ;}
    //else {TC25->Divide(5,5,0,0) ;}

    double ipad(0);
    double nOk=0;
    double PeakAverage=0.;
    
    TString tsFit,tsAvg,tsSig;
    tsFit.Form("%s_%d.png",h->GetName(),abs(offs));
    tsAvg.Form("%s_%d_Avg",h->GetName(),abs(offs));
    tsSig.Form("%s_%d_RMS",h->GetName(),abs(offs));
    TH2D *hAvg = new TH2D( tsAvg.Data(),tsAvg.Data(),25,0.5,25.5,50,0.95,1.05);
    TH2D *hSig = new TH2D( tsSig.Data(),tsSig.Data(),25,0.5,25.5,50,0.,0.1);
    tObjStack.push_back(hAvg);
    tObjStack.push_back(hSig);
    
    for( uint iCry=0; iCry< nCry ; iCry++ ) 
      {
	if( offs<0 && OnBorder[iCry] ) continue;

	auto TC25 = (TCanvas*) gROOT->FindObject("TC25");
	if (TC25) delete TC25;
	TC25 = new TCanvas("TC25",h->GetTitle(),2048,1024);
	TC25->cd();
	//TC25->cd(ipad+1);
	//ipad++;
	
	gPad->SetTickx(2);
	gPad->SetTicky(2);
	TH1D *hFit = h->ProjectionY("hFit",iCry+1,iCry+1);
	hFit->SetNameTitle("hFit","");

	gn[iCry]=0;
	if( hFit->GetEntries() < 10 )
	  {
	    hFit->DrawCopy("e") ;
	    continue ;
	  }
	if( hFit->GetEntries() < MinCounts  )
	  {
	    hFit->DrawCopy("e") ;
	  }
	else
	  {
	    nOk++;
	    FitHisto( hFit, fpar, ferr);
	    cout << " *** " << iCry << " " << fpar[1] << " " << ferr[1] << endl ;
	    hAvg->Fill( iCry+1, fpar[1]);
	    hAvg->SetMarkerStyle(20);
	    hAvg->SetMarkerSize(0.8);
	    // hSig->SetBinContent( iCry+1, 100*fpar[2]/fpar[1]);
	    // hSig->SetBinError( iCry+1, 100*ferr[2]/fpar[1]);
	    hSig->Fill( iCry+1,fpar[2]/fpar[1] );
	    hSig->SetMarkerStyle(21);
	    hSig->SetMarkerSize(0.8);
	    gn[iCry] = fpar[1];
	  }

 	
	TString s;
	s.Form("%s_%d.png",h->GetName(),int(iCry));
	TC25->SaveAs(OutputDir+s) ;
      }
    // TC25->SaveAs(OutputDir+tsFit) ;
  }
  // ============================================================================

private:

  double eBeam;
  bool NewBranch;
};

#endif

#ifdef SkAnalyzer_cxx

SkAnalyzer::SkAnalyzer(TString FileName, TString TreeName) : fChain(0) 
{
  cout << " Open File " << FileName.Data() << endl << endl ;
  TFile *f = new TFile(FileName.Data()) ;
  TTree *tree =new TTree(); 
  f->GetObject(TreeName.Data(),tree);
  //  f->GetObject(Skim,tree);
  Init(tree) ;
}
SkAnalyzer::SkAnalyzer(TTree *tree, TString TreeName) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0)
     {
       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/skim_230531/Chain100Center.root");
       if (!f || !f->IsOpen()) { f = new TFile("../data/skim_230531/Chain100Center.root"); }
       f->GetObject(TreeName.Data(),tree);
       //  f->GetObject(Skim,tree);
     }
   Init(tree);
}

SkAnalyzer::~SkAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t SkAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t SkAnalyzer::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void SkAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);


   fChain->SetBranchAddress("isLaser", &isLaser, &b_isLaser);
   fChain->SetBranchAddress("FitCry", FitCry, &b_FitCry);
   fChain->SetBranchAddress("FitTim", FitTim, &b_FitTim);
   //   fChain->SetBranchAddress("FitWlf", FitWlf, &b_FitWlf);
   //   fChain->SetBranchAddress("FitWrg", FitWrg, &b_FitWrg);
   fChain->SetBranchAddress("MaxCry", MaxCry, &b_MaxCry);
   fChain->SetBranchAddress("MaxTim", MaxTim, &b_MaxTim);
   fChain->SetBranchAddress("PedCry", PedCry, &b_PedCry);
   fChain->SetBranchAddress("RmsCry", RmsCry, &b_RmsCry);
   //   fChain->SetBranchAddress("SigCry", SigCry, &b_SigCry);

   static TString xClusterBranch("xCluster");
   cout << "Test new Branches " << endl;
   
   TBranch *b = (TBranch *) tree->GetListOfBranches()->FindObject(xClusterBranch);
   // if( b==NULL ) cout << " gheminga " << endl;
   // else cout << "ghe " << endl ;
   if( b )
     {
       fChain->SetBranchAddress("xCluster", &xCluster, &b_xCluster);
       fChain->SetBranchAddress("yCluster", &yCluster, &b_yCluster);
       fChain->SetBranchAddress("TrackXecal",&TrackXecal,&b_TrackXecal);
       fChain->SetBranchAddress("TrackYecal",&TrackYecal,&b_TrackYecal);
       NewBranch = true ;
       
       cout << " New Branches OK" << endl;
     }
   else
     {
       NewBranch = false;
       cout << " Miss new Branches " << endl;
     }
     
   
   Notify();
}

Bool_t SkAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SkAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SkAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SkAnalyzer_cxx

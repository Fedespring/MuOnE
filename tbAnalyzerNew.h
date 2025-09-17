/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 24 16:36:51 2022 by ROOT version 6.24/06
// from TTree simpleNtuple/simpleNtuple
// found on file: LocalDataFiles/6GeV_A2_20220725_151400.root
//////////////////////////////////////////////////////////

#ifndef tbAnalyzer_h
#define tbAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "iostream"
#include "cmath"
#include "math.h"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TF1.h>
#include <TSystem.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>

using namespace std;

class tbAnalyzer {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types
  UInt_t          iEvt;
  UShort_t        nSig;
  vector<unsigned short> *signal=nullptr;
  vector<int>     *signalOffsets=nullptr;
  vector<int>     *signalCounts=nullptr;
  vector<double>  *pedestal=nullptr;
  vector<double>  *pedRMS=nullptr;
  vector<double>  *signOverRMS=nullptr;
  
  // List of branches
  TBranch        *b_iEvt;   //!
  TBranch        *b_nSig;   //!
  TBranch        *b_signal;   //!
  TBranch        *b_signalOffsets;   //!
  TBranch        *b_signalCounts;   //!
  TBranch        *b_pedestal;   //!
  TBranch        *b_pedRMS;   //!
  TBranch        *b_signOverRMS;   //!


  
  //  tbAnalyzer(TChain *tChain=0 ) ;
  tbAnalyzer(TTree *tree=0, bool run2022=false );
  tbAnalyzer(TChain *aChain=0, bool run2022=false );
  tbAnalyzer(TString FileName, bool run2022=false );
  virtual ~tbAnalyzer();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual void     Loop( int N, TString T );
  virtual void     Loop( int nevt, int ofs, int step, TString T );
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  const double LaserThreshold=10000 ;
  const double SingleThreshold=500 ;
  //  static const int nADC=32, nCry=25, n2Cry=625 ;
  static const int nADC=32, nCry=27; // two more channels to account for laser lamp info.
  vector<int> *colMap, *rowMap;
  vector<int> cryMap, adcMap;

  bool OldRun = false;
  
  //  Maps for 2022 data do not contain laser lamp information
  int ADCMap2022[32]={1, 0,29,30,31, 3, 2,26,27,28, 6, 5, 4,24,25, 9, 8, 7,22,23,12,11,10,20,21,-1,-1,-1,-1,-1,-1,-1};
  int CrystalMap2022[32]={1, 0, 6, 5,12,11,10,17,16,15,22,21,20,-1,-1,-1,-1,-1,-1,-1,23,24,18,19,13,14, 7, 8, 9, 2, 3, 4};

  // MAPs for 2023 data: ADC 18 & 19 contain laser lamp information; they are mapped on Crystals 25 & 26 (!)
  int ADCMap[32]={1, 0,29,30,31, 3, 2,26,27,28, 6, 5, 4,24,25, 9, 8, 7,22,23,12,11,10,20,21,18,19,-1,-1,-1,-1,-1};
  int CrystalMap[32]={1, 0, 6, 5,12,11,10,17,16,15,22,21,20,-1,-1,-1,-1,-1,25,26,23,24,18,19,13,14, 7, 8, 9, 2, 3, 4};

  double sCry[nCry], mCry[nCry], aCry[nCry], bCry[nCry];
  
  int bySide[nCry]={ 1,1,1,1,1,
		     1,0,0,0,1,
		     1,0,0,0,1,
		     1,0,0,0,1,
		     1,1,1,1,1 }; //  1,1 };
  
  double CalCon[nCry], t_vec[nCry] ;
  double r_mat[nCry][nCry] ; 
  
  vector <TH1D*> th1Stack;
  vector <TH2D*> th2Stack;

  // ====================================================
  inline int nEntries() { return (int) fChain->GetEntriesFast(); }
  inline double GetMaxCry() { return CryMax;}
  inline double GetMaxTlc() { return TlcMax;}
  inline double GetMaxTgb() { return TgbMax;}
  inline double GetLocalMaximum() { return LocMax;}
  inline double GetGlobalMaximum() { return GlbMax;}
  // ====================================================
  virtual TH1D* DrawSingleCry(  int iCry=-1, int ievt=-1 )
  {
    if( iCry<0 || iCry>nCry )
      {
	cout << " Non existing crystal" << iCry << endl;
	return 0;
      }
    if( ievt<0 ) FillEventHistogram( ievt );
    auto hCry = (TH1D*) gROOT->FindObject("hCry"); if (hCry) delete hCry;
    hCry = hCrySgn->ProjectionX("hCry",iCry,iCry);
    hCry->DrawCopy();
    return hCry;
  }
  // ====================================================
  inline void FitArray( int ievt=-1 )
  {
    if( ievt<0 ) FillEventHistogram(abs(ievt));

    gStyle->SetOptFit(1);
    gStyle->SetOptStat(10);
    double fpar[4],ferr[4];
    
    fpar[0] = GlbMax;
    fpar[1] = TgbMax;
    fpar[2] = 1.00;
    fpar[3] = 1.50;
    for( int i=0; i<4; i++ ) ferr[i]=0.1*fpar[i];

    auto tc1 = (TCanvas*) gROOT->FindObject("tc1");
    if( tc1 ) delete tc1;
    tc1 = new TCanvas("tc1","Sum ADC counts",600,600);
    double s = SigCry(fpar,ferr,0,25,0);

    // fix the widths and the peak position
    ferr[1] = ferr[2] = ferr[3]=0;
    
    auto tc25 = (TCanvas*) gROOT->FindObject("tc25");
    if( tc25 ) delete tc25;
    tc25 = new TCanvas("tc25","Sum ADC counts",1600,900);
    tc25->Divide(5,5,0,0);
    
    for( int ip=2; ip<4; ip++) ferr[ip]=0;
    for( int i=0; i<nCry; i++ )
      {
	tc25->cd(i+1);
	double s = SigCry(fpar,ferr,i,i,0);
      }
  }
  // ====================================================
  inline void FitHisto( double* fpar, double *ferr, TH1D* hFit, bool plot=false )
  {
    int npar=4;
    auto fFit = (TF1*) gROOT->FindObject("fFit"); if (fFit) delete fFit;
    fFit=new TF1("fFit",fCMS,0.,128.,npar);    
    for( int i=0; i<npar; i++ )
      {
	fFit->SetParameter(i,fpar[i]);
	fFit->SetParError(i,ferr[i]);
	if( ferr[i]==0 ) { fFit->FixParameter(i,fpar[i]);}
      }
    TString sFitOpt="LQN";
    if( plot ) sFitOpt="LQ";

    int ixmin = (int) fpar[1]-1;
    int ixmax = (int) fpar[1]+3;
    double xmin = (double) ixmin;
    double xmax = (double) ixmax;
    double xmean = (xmax+xmin)/2;
    fFit->SetParameter(0,hFit->GetBinContent( (int) xmean) );
    TFitResultPtr r = hFit->Fit( fFit,sFitOpt,"",xmin,xmax ) ;
    if( plot ) { gSystem->ProcessEvents();}
    for( int i=0; i<npar; i++) { fpar[i]=fFit->GetParameter(i); ferr[i]=fFit->GetParError(i); }
    //if( plot ) delete hFit;
    
  }  // ====================================================
  inline void FitPedes( double* fpar, double *ferr, TH1D* hFit, bool plot=false )
  {
    int npar=5;
    auto fFit = (TF1*) gROOT->FindObject("fFit"); if (fFit) delete fFit;
    fFit=new TF1("fFit",fCMS0,0.,128.,npar);    
    for( int i=0; i<npar; i++ )
      {
	fFit->SetParameter(i,fpar[i]);
	fFit->SetParError(i,ferr[i]);
	if( ferr[i]==0 ) { fFit->FixParameter(i,fpar[i]);}
      }
    TString sFitOpt="LQN";
    if( plot ) sFitOpt="LQ";

    double xmin = fpar[1]-5;
    double xmax = fpar[1]+5;
    //  double xmean = (xmax+xmin)/2;
    //  fFit->SetParameter(0,hFit->GetBinContent( (int) xmean) );
    TFitResultPtr r = hFit->Fit( fFit,sFitOpt,"",xmin,xmax ) ;
    if( plot ) { gSystem->ProcessEvents();}
    for( int i=0; i<npar; i++) { fpar[i]=fFit->GetParameter(i); ferr[i]=fFit->GetParError(i); }

  }
  // ====================================================
  inline double SigCry( int imin=0, int imax=25, int ievt=0, bool plot=false )
  {
    if( ievt<0 ) FillEventHistogram(ievt);
    double fpar[4]={0.0,0.0,0.5,1.0};
    double ferr[4]={1.0,1.0,0.1,0.1};
    return SigCry( fpar, ferr, imin, imax, abs(ievt), plot );
  }
  // ====================================================
  inline double SigCry( double* fpar, double *ferr, int imin=0, int imax=25, int ievt=0, bool plot=false )
  {
    if( ievt<0 ) FillEventHistogram(ievt);
    auto hFit = (TH1D*) gROOT->FindObject("hFit"); if (hFit) delete hFit;
    hFit = hCrySgn->ProjectionX("hFit",imin,imax);
    
    if( fpar[0] == 0 ) fpar[0] = hFit->GetMaximum();
    if( fpar[1] == 0 ) fpar[1] = hFit->GetMaximumBin()-1;
 
    FitHisto( fpar, ferr, hFit, plot );
    return fpar[0];
  }
  // ====================================================
  inline TH2D* GetLampHistogram() { return hLmpSgn; }
  // ====================================================
  inline TH2D* GetEventHistogram(){ return hCrySgn; }
  // ====================================================
  inline void FillEventHistogram( int ievt=-1 )
  {
    int nevt = fChain->GetEntriesFast();

    if( abs(ievt) >= nevt ) { cout << " Event " << ievt << " does not exist \n"; return; }
    if( ievt<0 )
      {
	int kevt = LoadTree( abs(ievt) ) ; 
	int nb = fChain->GetEntry(kevt) ;
      }
    hCrySgn->Reset();
    
    double pval[100], prms[100];
    double pmin[25];
    ComputePedestal( 0, 30, pval, prms );
    //ComputePedestalMin(0, 127, pmin);

    //  Lamp[0] = signal->at(signalOffsets->at(16)+iofs)-pedestal->at(16) ;
    //  Lamp[1] = signal->at(signalOffsets->at(17)+iofs)-pedestal->at(17) ;
    
    for( uint ichan=0; ichan<pedestal->size(); ichan++ )
      {
	// if( CrystalMap[ichan]==-1 ) continue;
	// int j = CrystalMap[ichan] ;

	if( cryMap.at(ichan) == -1 ) continue;
	
	int j = cryMap.at(ichan);	
	pedestal->at(j) = pval[ichan];
	pedRMS->at(j) = prms[ichan];
	for( int iofs=0; iofs<128; iofs++ )
	  {
	    // double val = signal->at(signalOffsets->at(ichan)+iofs)-pval[ichan] ;
	    double val = signal->at(signalOffsets->at(ichan)+iofs) ;
	    hCrySgn->SetBinContent( iofs+1,j+1, val) ;
	  }
      }

    
    int b,c;
    LocMax = hCrySgn->GetMaximum();
    GlbMax = hCrySgn->ProjectionX()->GetMaximum();
    hCrySgn->GetMaximumBin(TlcMax,CryMax,c);
    hCrySgn->ProjectionX()->GetMaximumBin(TgbMax,b,c) ;

    if( ievt<0 ) hCrySgn->DrawCopy("lego");
  }
  // ====================================================
  inline void ComputePedestal( int imin, int imax, double *pedval, double *pedrms )
  {
    imin = max( imin,0 ) ;
    imax = min( imax,127 ) ;
    for( uint ichan=0; ichan<pedestal->size(); ichan++ ) {

      //      if( CrystalMap[ichan] == -1 ) continue ;
      if( cryMap.at(ichan) == -1 ) continue;
      int nval(0);
      double p_ave(0);
      double p_val(0);
      double p_sqr(0);
      for( int iofs=imin; iofs<imax; iofs++ ) {
	p_val = signal->at(signalOffsets->at(ichan)+iofs) ;
	p_ave += p_val; 
	p_sqr += p_val*p_val;
	nval++;
      }
      p_ave /= nval;
      p_sqr /= nval;
      p_sqr -= p_ave*p_ave;
      p_sqr = sqrt(p_sqr);
      //printf(" pedestal %d %d %f %f \n",ichan,nval,p_ave,p_sqr); 
      pedval[ichan] = p_ave ;
      pedrms[ichan] = p_sqr ;
    }
  }  

  // Computes the pedestal following the same structure but 
  inline void ComputePedestalMin(int imin, int imax, float *pedmin) {
    imin = std::max(imin, 0);
    imax = std::min(imax, 127);

    for (uint ichan = 0; ichan < pedestal->size(); ichan++) {
        if (cryMap.at(ichan) == -1) continue;

        double p_min = std::numeric_limits<double>::max(); // Initialize to a large value

        for (int iofs = imin; iofs < imax; iofs++) {
            double p_val = signal->at(signalOffsets->at(ichan) + iofs);
            p_min = std::min(p_min, p_val); // Update if a smaller value is found
        }
	int j = CrystalMap[ichan];
        pedmin[j] = p_min;  // Store the minimum value
    }
  }
  // ====================================================
  virtual void LoadRelativeGain( double *gain )
  {
    TString fName = "CrystalGain.hist" ;
    if( !gSystem->AccessPathName(fName.Data()) ) 
      {
	TFile *fGain=0; 
	fGain = TFile::Open(fName.Data());
	TH1D *hGain = (TH1D*) fGain->Get( "hRes1" ) ;
	if( hGain!=0 ) { 
	  for( int i=0; i<hGain->GetNbinsX(); i++) {
	    double r =(hGain->GetBinContent(i+1)==0 ? 1 : hGain->GetBinContent(13)/hGain->GetBinContent(i+1));
	    *gain = r*(*gain);
	    gain++;
	  }
	}
	else cout << " histogram not found " << endl ;
	fGain->Close() ;
      }
    else
      { cout << " Gain file not found " << fName.Data() << endl ;}
  }
  
  // ====================================================
  static double biGaus( double *t, double *par )
  {
    double x= t[0] ;
    double f(0) ;
    if( x<par[1] )
      { f=par[0]*TMath::Gaus(x,par[1],par[2]) ;}
    else
      { f=par[0]*TMath::Gaus(x,par[1],par[3]) ;}
    return f;
  }
  // ====================================================
  inline void FillTimeSlot( TH2D *h )
  {
    for( uint iCry=0 ; iCry<nCry; iCry++ ) 
      {
	int iADC = adcMap.at(iCry) ;
	for( int i=0; i<signalCounts->at(iADC) ; i++ )
	  {h->Fill( i+0.5, (double) iCry, signal->at(signalOffsets->at(iADC)+i)-pedestal->at(iADC));}
      }
  }
  // ====================================================
  virtual bool TagLaser() // { return GlbMax>LaserThreshold ;}
  {
    if( GlbMax<LaserThreshold ) return false;
    for( int i=1; i<nCry+1; i++ )
      { if( hCrySgn->ProjectionX("",i,i)->GetMaximum()< SingleThreshold ) return false; }
	//	cout << " " << i << " " <<  hCrySgn->ProjectionX("",i,i)->GetMaximum() <<endl;
    return true;
  }
  // ====================================================
  static double fCMS( double *t, double *par )
  {
    double x = t[0]-par[1];
    double a = par[2];
    double b = par[3];
    //    cout << x << " " << a << "  " << b << endl ;
    if( a==0 || b == 0 ) return 0;
    double arg = 1+x/a/b;
    if( arg<0 ) return 0;
    return par[0]*pow(arg,a)/exp(x/b);
  }
  // ====================================================
  static double fCMS0( double *t, double *par )
  {
    double x = t[0]-par[1];
    double a = par[2];
    double b = par[3];
    //    cout << x << " " << a << "  " << b << endl ;
    if( a==0 || b == 0 ) return 0;
    double arg = 1+x/a/b;
    if( arg<0 ) return par[4];
    return par[0]*pow(arg,a)/exp(x/b)+par[4];
  }
  // ====================================================

  
private:

  // May 23 set up
  int TimeSlotMin=39;
  int TimeSlotMax=43;

  // int TimeSlotMin=76;
  // int TimeSlotMax=83;
  TH2D *hCrySgn, *hLmpSgn;
  int CryMax=-1;
  int TlcMax=-1;
  int TgbMax=-1;
  double LocMax=0, GlbMax=0;
};

#endif

#ifdef tbAnalyzer_cxx

tbAnalyzer::tbAnalyzer(TString FileName, bool run2022 ) : fChain(0) 
{
  cout << " Open File " << FileName.Data() << endl << endl ;
  TFile *f = new TFile(FileName.Data()) ;
  TTree *tree =new TTree(); 
  f->GetObject("simpleNtuple",tree);
  colMap = f->Get<vector<int>>("colMap") ;
  rowMap = f->Get<vector<int>>("rowMap") ;
  OldRun = run2022;
  
  Init(tree) ;
}

tbAnalyzer::tbAnalyzer( TTree *tree, bool run2022 ) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  OldRun = run2022;
  if (tree == 0)
    { cout << "   ====   Tree not Found ==== " << endl ;}
  else
    { Init(tree); }
}
tbAnalyzer::tbAnalyzer(TChain *aChain, bool run2022 ) : fChain(0) 
{
  OldRun = run2022;
  TTree *tree = new TTree();
   if (tree == 0)
     { cout << "   ====   Tree not Found ==== " << endl ;}
   else
     {
       tree = aChain;
       tree->Print();
       Init( tree ) ;		       
     }
}

tbAnalyzer::~tbAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tbAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tbAnalyzer::LoadTree(Long64_t entry)
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

void tbAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
  // tree->Print();

   // Set object pointer
   signal = 0;
   signalOffsets = 0;
   signalCounts = 0;
   pedestal = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("iEvt", &iEvt, &b_iEvt);
   fChain->SetBranchAddress("nSig", &nSig, &b_nSig);
   fChain->SetBranchAddress("signal", &signal, &b_signal);
   fChain->SetBranchAddress("signalOffsets", &signalOffsets, &b_signalOffsets);
   fChain->SetBranchAddress("signalCounts", &signalCounts, &b_signalCounts);
   fChain->SetBranchAddress("pedestal", &pedestal, &b_pedestal);
   fChain->SetBranchAddress("pedRMS", &pedRMS, &b_pedRMS);
   fChain->SetBranchAddress("signOverRMS", &signOverRMS, &b_signOverRMS);
   Notify();

   if( OldRun )
     {
       for( uint iADC=0; iADC<nADC; iADC++ ) { cryMap.push_back( CrystalMap2022[iADC] ) ; adcMap.push_back( ADCMap2022[iADC] ) ;}
     }
   else 
     {
       for( uint iADC=0; iADC<nADC; iADC++ ) { cryMap.push_back( CrystalMap[iADC] ) ; adcMap.push_back( ADCMap[iADC] ) ;}
     }

   cout << " ***** Total entries ===== " << fChain->GetEntries() << endl ;
   cout << " ***** Total files ===== " << fChain->GetTreeNumber() << endl ;
   cout << " ***** ADC Map: " << endl ;
   cout << endl << " Cry " ;
   for( uint i=0; i<nADC; i++ ) printf("%2d,",cryMap.at(i)) ;
   cout << endl << " ADC " ;
   for( uint i=0; i<nADC; i++ ) printf("%2d,",adcMap.at(i)) ;
   cout << endl << endl ;   
   cout << " *******  Signal ADC time slut : " << TimeSlotMin << "   ****   " << TimeSlotMax << endl << endl ;
   // Create the event histogram
   hLmpSgn = new TH2D("hLmpSgn","Lamp signal",128,0.,128,4,0.,4.);
   hCrySgn = new TH2D("hCrySgn","ADC counts per Crystal",128,0.,128,nCry,0.,(float) nCry);
   cout << " *** Created Histogram " << hCrySgn->GetTitle() << endl ;
}

Bool_t tbAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tbAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tbAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tbAnalyzer_cxx

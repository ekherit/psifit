#define _HESSE_ 1 // calculate errors with HESSE
#define _MINOs_ 1 // use minos errors
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <complex>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <fstream>
#include <unistd.h>
#include <error.h>
#include <argp.h>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include <TROOT.h>
#include <TH1.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TAxis.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TPad.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TFormula.h>
#include <TString.h>
#include <TVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <TString.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TVirtualFitter.h>
#include <TCanvas.h>
#include <TPostScript.h>
#include <TApplication.h>

R__EXTERN TSystem *gSystem;
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("ROOTPROGRAM","ROOTPROGRAM", initfuncs);
using namespace std;

#include "FitOniumRCompactLib.hh"

#define  NumMaxP 120
Int_t    NpPP=0;
Double_t CrossBhabha=265;
Double_t CrossGG=50;
Double_t LUM_CROSS_SECTION=0;

Double_t CrossSInScan[NumMaxP];
Double_t CrossSErrInScan[NumMaxP];
Double_t CrossSBBInScan[NumMaxP];
Double_t CrossSBBErrInScan[NumMaxP];
Double_t CrossSGGInScan[NumMaxP];
Double_t CrossSGGErrInScan[NumMaxP];
Double_t LumInScan[NumMaxP];
Double_t LumInScanGG[NumMaxP];
Double_t LumInScanEE[NumMaxP];
Double_t LumInScanEEGG[NumMaxP];
Double_t NmhInScan[NumMaxP];
Double_t LumLgammaInScan[NumMaxP];
Double_t NbbInScan[NumMaxP];
Double_t NggInScan[NumMaxP];
Double_t NLum[NumMaxP];
Double_t EInScan[NumMaxP];
Double_t SigmaWInScan[NumMaxP];
Double_t dSigmaWInScan[NumMaxP];
Double_t WInScan[NumMaxP];
Double_t EErrInScan[NumMaxP];
Double_t WErrInScan[NumMaxP];
Double_t SignalDiscrepancy[NumMaxP];
Double_t SignalDiscrepancyError[NumMaxP];
Double_t BBLumCorr[NumMaxP];

Int_t    NumEpoints=0;
Double_t MinChi2=1e+7;
void fcnResMult    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnResMult2    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnResMult_both2    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnResMult3    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

bool USE_CBS_SIGMAW = false; 
bool USE_CHI2 = false;
bool LUM_COR = true;
bool FREE_ENERGY_FIT=false;
unsigned BOTH_FIT=0;
double DEDIFF=0.02;

enum LuminosityType
{
  BESLUM, //BES online luminosity monitor
  NEELUM, //Bhabha scattering
  NGGLUM, //Gamma gamma production luminosity
  NEEGGLUM,  //both bhabha and gamma gamma luminosity
  NMHEELUM,  //Resonance and bhabha together fit
  NMHGGLUM,  //Resonance and gamma gamma together fit
  NMHEEGGLUM //Resonance, bhabha and gg together fit
};

LuminosityType LUMINOSITY;
double PDGMASS=0;

enum ResonanceType
{
  JPSIRES,
  PSI2SRES
};

unsigned RESONANCE;


double CHI2_TOTAL;
double CHI2_LUM;
double CHI2_ENERGY;
double CHI2_SIGNAL;




inline double sq(double x) { return x*x;}
inline double cb(double x) { return x*x*x;}

int main(int argc, char **argv)
{
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  opt_desc.add_options()
    ("help,h","Print this help")
    ("verbose,v","Verbose output")
    ("scan,s",po::value<std::string>()->default_value("scan.txt"),"data file")
    ("lum,l", po::value<std::string>()->default_value("ee"),"luminosity used: bes, ee (default), gg, eegg, mhee, mhgg, mheegg")
    ("resonance,R", po::value<unsigned>(&RESONANCE)->default_value(JPSIRES),"resonance: 0 - JPSI, 1 - PSI2S")
    ("free-energy,E","allow energy points to be variated")
    ("both", po::value<unsigned>(&BOTH_FIT)->default_value(0),"Fit mode 0,1,2")
    ("lum-cor","Use MC luminosity correction")
    ("chi2", "Use chi2 fit")
    ("cross-section-ee",po::value<double>(&CrossBhabha)->default_value(200), "Bhabha cross section, nb")
    ("cross-section-gg", po::value<double>(&CrossGG)->default_value(20), "Gamma-gamma cross section, nb")
    ("dE",po::value<double>(&DEDIFF)->default_value(0.02),"Combine run with energy in dE interval, MeV")
    ("skip",po::value< std::vector<unsigned> >(),"list of skipped points")
    ("ems-error",po::value< double >(),"ems energy measurement error for each point")
    ;
  po::positional_options_description pos;
  pos.add("scan",-1);
  po::variables_map opt; //options container
  try
  {
    po::store(po::command_line_parser(argc, argv).options(opt_desc).positional(pos).run(), opt);
    std::ifstream config_file("fit.cfg");
    po::store(po::parse_config_file(config_file,opt_desc,true), opt);
    po::notify(opt);
  } 
  catch (boost::program_options::error & po_error)
  {
    cerr << "WARGNING: configuration: "<< po_error.what() << endl;
  }

  if(opt.count("help"))
  {
    std::clog << opt_desc;
    return 0;
  }
  std::cout << "Fit resonance: ";
  switch(RESONANCE)
  {
    case JPSIRES:
      PDGMASS=_MJPsi;
      std::cout << "JPSI";
      break;
    case PSI2SRES:
      PDGMASS=_MPsiPrime;
      std::cout << "PSI2S";
      break;
  }
  std::cout << std::endl;


  USE_CHI2 = opt.count("chi2");
  std::cout << boolalpha;

  std::cout << "Chi2 fit: " << USE_CHI2 << std::endl;

  LUM_COR = opt.count("lum-cor") && LUMINOSITY==NEELUM;
  std::cout << "Luminosity corrections: " << LUM_COR << std::endl;

  FREE_ENERGY_FIT = opt.count("free-energy");
  std::cout << "Use free energy fit: " << FREE_ENERGY_FIT << std::endl;

  std::cout << "Average cross section: Bhabha: " << CrossBhabha << " nb, Gamma-gamma: " << CrossGG << " nb" << std::endl;

	cout << "Luminosity used: ";
  std::string lumstr=opt["lum"].as<string>();
  if(lumstr=="bes")
  {
      LUMINOSITY = BESLUM;
      LUM_CROSS_SECTION=CrossBhabha;
      std::cout << "BES online luminosity monitor"<< std::endl;
  }
  if(lumstr=="ee")
  {
      LUMINOSITY = NEELUM;
      LUM_CROSS_SECTION=CrossBhabha;
      std::cout << "Bhabha"<< std::endl;
  }
  if(lumstr=="gg")
  {
      LUMINOSITY = NGGLUM;
      LUM_CROSS_SECTION=CrossGG;
      std::cout << "gamma-gamma"<< std::endl;
  }
  if(lumstr=="eegg")
  {
      LUMINOSITY = NEEGGLUM;
      LUM_CROSS_SECTION=CrossGG+CrossBhabha;
      std::cout << "Bhabha + gamma-gamma"<< std::endl;
  }
  if(lumstr=="mhee")
  {
      LUMINOSITY = NMHEELUM;
      LUM_CROSS_SECTION=CrossBhabha;
      std::cout << "MH + Bhabha "<< std::endl;
  }
  if(lumstr=="mhgg")
  {
      LUMINOSITY = NMHGGLUM;
      LUM_CROSS_SECTION=CrossGG;
      std::cout << "MH + gamma-gamma "<< std::endl;
  }
  if(lumstr=="mheegg")
  {
      LUMINOSITY = NMHEEGGLUM;
      LUM_CROSS_SECTION=CrossGG+CrossBhabha;
      std::cout << "MH + Bhabha + gamma-gamma "<< std::endl;
  }
  cout << "Average cross section of luminosity measurement process: " << LUM_CROSS_SECTION << " nb" << endl;
  TApplication* theApp=0;
  theApp=new  TApplication("App", &argc, argv);


  TMinuit*  MinuitRes=0; 
  TF1* FitRes=0;
  TF1* FitResBG=0;
  TGraphErrors* GrRes=0;
  TLine* LineRes=0;
  int dimMHFile=10;
  int MHRun=0;
  int MHLum=1;
  int MHEnergy=2;
  int MHEnergyErr=3;
  int MHSigmaW=4;
  int MHdSigmaW=5;
  int MHNmh=6; 
  int MHNee=7;
  int MHNgg=8; 
  int MHLumCor=9; 
  double** AllMH=0;
  int npMHFile=0; 

  /* ***** Reading data from file ****************** */
  AllMH=new double* [npMHFile];
  std::string data_file_name=opt["scan"].as<std::string>();
  std::cout << "Reading data from file " << data_file_name;
  npMHFile=GetNumRows(data_file_name.c_str(),dimMHFile);       
  FillArrayFromFile(data_file_name.c_str(),AllMH,dimMHFile,npMHFile);  
  std::cout << "Read " << npMHFile << " points." << std::endl;

  int dimAP= 19;
  int ARun=0;
  int AEnergy =1;
  int AEnergyErr=2;
  int AMHEv  =3;
  int AMHEvH =4;
  int AEE    =5;
  int AGG    =6;
  int ALe    =7;
  int ALp    =8;
  int ACr    =9;
  int ACrErr =10;
  int ALRatE  =11;
  int ALRatP  =12;
  int ALRatEErr  =13;
  int ALRatPErr  =14;
  int ASE  =15;
  int ASigmaW = 16;
  int AdSigmaW = 17;
	int ALumCor = 18;
  int npAP=0;    
  npAP=npMHFile;
  double** AP=new double* [npAP];
  int Aind=0;

  for(int i=0;i<npMHFile;i++ )
  {      
    AP[Aind]=new double [dimAP];    
    AP[Aind][ARun]=AllMH[i][MHRun];                    
    AP[Aind][AEnergy]=AllMH[i][MHEnergy]*0.5;        
    AP[Aind][AEnergyErr]=AllMH[i][MHEnergyErr]*0.5;       
    AP[Aind][ALe]=AllMH[i][MHLum];//Fill luminocity from bes online lum
    AP[Aind][ALp]=AllMH[i][MHLum];//This is the same for both electron and positron
    AP[Aind][AMHEv]=AllMH[i][MHNmh];
    AP[Aind][ASigmaW] = AllMH[i][MHSigmaW];
    AP[Aind][AdSigmaW] = AllMH[i][MHdSigmaW];
    AP[Aind][AEE]=AllMH[i][MHNee];
    AP[Aind][AGG]=AllMH[i][MHNgg];
    AP[Aind][ALumCor]=AllMH[i][MHLumCor];
    Aind++;
  }

  if(Aind!=npAP) cout<<"PROBLEM !!!!!!!!!"<<endl;
  int   NumParForC[3];  
  int   CUpDown[3];
  NumParForC[0] = AEnergy;
  CUpDown[0]  = 1;
  NumParForC[1] = ARun;
  CUpDown[1]  = 1;

  /*  SKIP SORTING OF POINTS IF separate scan*/
  if(BOTH_FIT==0)
  {
    for(Int_t i=0;i<npAP;i++) // simplest bubble's method for each condition
    {
      for(Int_t j=0;j<i;j++)
      {
        for(Int_t ic=1;ic>=0;ic--)
        {

          if(comparDRows(CUpDown[ic],AP[i],AP[j],NumParForC[ic])<0)
          {
            swapDRows(AP[i],AP[j],dimAP);

          }
        }
      }
    }
  }

  char ouputstring[140];
  ofstream AA("test.txt", ios::out);
  for(int i=0;i<npAP;i++)
  {
    sprintf(ouputstring,"  %7.4f %.0f %.0f %.0f",AP[i][AEnergy],AP[i][ARun],AP[i][AMHEv],AP[i][AEE]); 
    AA<<ouputstring<<"\n";	            
  }
  AA.close();        
  Double_t EMin=AP[0][AEnergy]-1 ;      

  Double_t EMax=AP[npAP-1][AEnergy]+1;
  int np=npAP ;   

  Double_t *En_=new Double_t[np];
  Double_t *Eerr_=new Double_t[np];
  Double_t *SigmaW_=new Double_t[np];
  Double_t *dSigmaW_=new Double_t[np];
  Double_t *Nmh_=new Double_t[np];
  Double_t *Nbb_=new Double_t[np];
  Double_t *Ngg_=new Double_t[np];
  Double_t *Le_=new Double_t[np];
  Double_t *Lp_=new Double_t[np];    
  Double_t *Lcor_=new Double_t[np];    
  Int_t*        Euse=new Int_t[np];
  np=0;

  Double_t SigmaW=1.0;
  for(int i=0;i<npAP;i++)
  {   
    En_[np]=AP[i][AEnergy];
    Eerr_[np]=AP[i][AEnergyErr];
    Nmh_[np]=AP[i][AMHEv];
    Nbb_[np]=AP[i][AEE];  //here will be Nee or Ngg
    Ngg_[np]=AP[i][AGG];  //here will be Nee or Ngg
    Le_[np]=AP[i][ALe];   //here will be BES online lum
    Lp_[np]=AP[i][ALp];  
    SigmaW_[np]=AP[i][ASigmaW];
    dSigmaW_[np]=AP[i][AdSigmaW];
		Lcor_[np] = AP[i][ALumCor];
    np++;      
  }
  int  NEp=0;
  SeparatePointsPartNew(np,Nbb_,&NEp,Euse,En_,DEDIFF); 
  Double_t *En  =new Double_t[NEp];
  Double_t *Eerr=new Double_t[NEp];
  Double_t *SW  =new Double_t[NEp];
  Double_t *dSW =new Double_t[NEp];
  Double_t *Nmh=new  Double_t[NEp];
  Double_t *Nbb=new  Double_t[NEp];
  Double_t *Ngg=new  Double_t[NEp];
  Double_t *Le=new  Double_t[NEp];
  Double_t *Lp=new  Double_t[NEp];
  Double_t *Lcor=new  Double_t[NEp];

  SumPointsByQuant(np,NEp,Euse,En_,Nbb_,Eerr_,En,Eerr,false);
  SumPointsSimple(np,NEp,Euse,Nmh_,Nmh);
  SumPointsSimple(np,NEp,Euse,Nbb_,Nbb); 
  SumPointsSimple(np,NEp,Euse,Ngg_,Ngg); 
  SumPointsSimple(np,NEp,Euse,Le_,Le); 
  SumPointsSimple(np,NEp,Euse,Lp_,Lp);     
  SumPointsSimple(np,NEp,Euse,SigmaW_,SW);     
  SumPointsSimple(np,NEp,Euse,dSigmaW_,dSW);     
  SumPointsSimple(np,NEp,Euse,Lcor_,Lcor);     
  NumEpoints=NEp;
  int numpar=4;
  if(FREE_ENERGY_FIT) numpar+=NEp;
  if(BOTH_FIT==1) numpar+=2;
  if(BOTH_FIT==2) numpar+=3;
  Double_t ECorrBB=0;
  Double_t ECorrGG=0;
  Double_t LG=0,Lee=0,Lgg=0; //integrated luminosity measured different way
  //fill global arrayes which needed to access from the FCN function
      double DEE[7] = {-2.27063e-04,
            1.12744e-01   ,
           -9.56351e-02   ,
            1.04181e-02   ,
            1.16408e-01   ,
           -9.23315e-02   ,
           -8.21663e-03,
      };
      double DEE2[7] = {
          -2.27063e-04,   
           1.12744e-01,   
          -9.56351e-02,   
           1.04181e-02,   
           1.16408e-01,   
          -9.23315e-02,   
          -8.21663e-03
      };
  double ems_error=0;
  if(opt.count("ems-error")) 
  {
    ems_error = opt["ems-error"].as<double>();
    cout << "Set all ems errors to : " << ems_error << endl;
  }
  TRandom random;
  for(int is=0;is<NEp;is++)
  {
    double rn = random.Gaus(0,ems_error);
    double dE=DEE2[is]+ rn;
    dE=0;
    cout << is << " " << dE <<  " " << rn << endl;
    EInScan[is]=En[is]+dE;
    WInScan[is]=2.*EInScan[is];
    if(opt.count("ems-error")) EErrInScan[is]=ems_error;
    else EErrInScan[is]=Eerr[is];
    WErrInScan[is]=EErrInScan[is]*2;
    SigmaWInScan[is]=SW[is];
    dSigmaWInScan[is]=dSW[is];
    NmhInScan[is]=Nmh[is];   
    NbbInScan[is]=Nbb[is];       
    NggInScan[is]=Ngg[is];       
    ECorrBB=1./CrossSBhabhaPP(En[is],&CrossBhabha);
    ECorrGG=1./CrossSBhabhaPP(En[is],&CrossGG);
    LumInScanEE[is]=NbbInScan[is]*ECorrBB;
    LumInScanGG[is]=NggInScan[is]*ECorrGG;
    LumInScanEEGG[is]=(NggInScan[is] + NbbInScan[is])/(1./ECorrBB + 1./ECorrGG);
    LG+=TMath::Max(Le[is],Lp[is]);
    Lee+= LumInScanEE[is];
    Lgg+= LumInScanGG[is];
    LumLgammaInScan[is]=TMath::Max(Le[is],Lp[is]); //BES lum
		BBLumCorr[is] = Lcor[is];
    switch(LUMINOSITY)
    {
      case BESLUM:
        LumInScan[is]=LumLgammaInScan[is];
        NLum[is]=1e100;
        break;
      case NEELUM:
        LumInScan[is]=LumInScanEE[is];
        NLum[is] = NbbInScan[is];
        break;
      case NGGLUM:
        LumInScan[is]=LumInScanGG[is];
        NLum[is]=NggInScan[is];
        break;
      case NEEGGLUM:
        LumInScan[is]=LumInScanEEGG[is];
        NLum[is]=NggInScan[is]+NbbInScan[is];
        break;
      case NMHEELUM:
        LumInScan[is]=LumInScanEE[is];
        NLum[is]=NbbInScan[is];
        break;
      case NMHGGLUM:
        LumInScan[is]=LumInScanGG[is];
        NLum[is]=NggInScan[is];
        break;
      case NMHEEGGLUM:
        LumInScan[is]=LumInScanGG[is];
        NLum[is]=NggInScan[is]+NbbInScan[is];
        break;
    }
    //cout<<"BESLUM:"<<LumLgammaInScan[is]<<" LumInScan:"<<LumInScan[is]<< ", BBLumCor=" << BBLumCorr[is] << endl;;
    //calculate mhadr cross section
    CrossSInScan[is]=Nmh[is]/LumInScan[is];
    if(Nmh[is]>4) CrossSErrInScan[is]=sqrt(Nmh[is]*(1.+Nmh[is]/NLum[is]))/LumInScan[is];        
    else CrossSErrInScan[is]=sqrt(Nmh[is]+4)/LumInScan[is];     
    //calculate bhabha or gg cross section
    CrossSBBInScan[is]=Nbb[is]/LumInScan[is];
    CrossSGGInScan[is]=Ngg[is]/LumInScan[is];
    if(Nbb[is]>4) CrossSBBErrInScan[is]=sqrt(Nbb[is])/LumInScan[is];        
    else CrossSBBErrInScan[is]=sqrt(Nbb[is]+4)/LumInScan[is];         
    if(Ngg[is]>4) CrossSGGErrInScan[is]=sqrt(Ngg[is])/LumInScan[is];        
    else CrossSGGErrInScan[is]=sqrt(Ngg[is]+4)/LumInScan[is];         
  }
  std::cout << "Total luminosity: Lbes=" << LG << "/nb, Lee=" << Lee << "/nb, Lgg="<<Lgg<<"/nb"<<std::endl;
  if(opt.count("verbose"))
  {
    for(int is=0;is<NEp;is++)
    {
      cout<<"point:"<<is<<" Energy:"<<En[is]<<"NMH:"<<Nmh[is]<<endl;	    
    }
  }    

  MinuitRes = new TMinuit(numpar);
  MinuitRes->SetFCN(fcnResMult);  

  Double_t arglistRes[numpar*2];

  Int_t ierflgRes = 0;
  if(opt.count("verbose"))
  {
    arglistRes[0]=-1;
    MinuitRes->mnexcm("SET PRINT", arglistRes,1,ierflgRes);
    arglistRes[0] =0;
    MinuitRes->mnexcm("SET NOW", arglistRes,0,ierflgRes);
  }
  else
  {
    arglistRes[0]=-1;
    MinuitRes->mnexcm("SET PRINT", arglistRes,1,ierflgRes);
    arglistRes[0] = 0;
    MinuitRes->mnexcm("SET NOW", arglistRes,0,ierflgRes);
  }
  MinuitRes->mnexcm("SET ERR", arglistRes,1,ierflgRes);
  arglistRes[0] = 2;
  MinuitRes->mnexcm("SET STRATEGY", arglistRes,1,ierflgRes);

  Double_t vstartRes[5]= {5,0.34,0.5,1.59,LUM_CROSS_SECTION};   

  Double_t stepRes[5] =  {1,0.1,0.01,0.05,0.0};



  MinuitRes->SetMaxIterations(100000000);                         
  MinuitRes->DefineParameter(0,"bg",vstartRes[0],stepRes[0],-100,+100);
  MinuitRes->DefineParameter(1,"eff",vstartRes[1],stepRes[1],0.0,1.0);      
  MinuitRes->DefineParameter(2,"dM/2.",vstartRes[2],stepRes[2],-0.5,0.5);      
  MinuitRes->DefineParameter(3,"SigmaW",vstartRes[3],stepRes[3],0.5,2.0);
  //if(USE_CBS_SIGMAW) MinuitRes->FixParameter(3);
  if(FREE_ENERGY_FIT)
  {
    for(int j=0;j<NEp;j++)
    {
      char  NameP[10];
      sprintf(NameP,"dE%d",j);         
      MinuitRes->DefineParameter(j+4,NameP,0,0.1,-0.5,+0.5);        
    }
  }

  int nep = FREE_ENERGY_FIT==false ? 0 : NEp;
  /* 
   *  Fit both scan togegher
   *  This code must be refubrished to use with boost program options.
  if(arguments.scan==3)
  {
    cout << "Fit scan1 and scan2"<< endl;
    if(BOTH_FIT==1)
    {
      MinuitRes->DefineParameter(0+4+nep,"bg2",vstartRes[0],stepRes[0],-150,150.0);
      MinuitRes->DefineParameter(1+4+nep,"eff2",vstartRes[1],stepRes[1],0.1,1.0);      
      cout << "Variant with own bg and eff" << endl;
    }
    if(BOTH_FIT==2)
    {
      MinuitRes->DefineParameter(0+4+nep,"bg2",vstartRes[0],stepRes[0],-150,150.0);
      MinuitRes->DefineParameter(1+4+nep,"eff2",vstartRes[1],stepRes[1],0.1,1.0);      
      MinuitRes->DefineParameter(2+4+nep,"SigmaW2",vstartRes[3],stepRes[3],0.5,1.8);      
      cout << "Variant with own bg, eff and sigmaW" << endl;
    }
  }
  */

  MinuitRes->mnexcm("MIGRAD", arglistRes,numpar,ierflgRes);
  MinuitRes->mnimpr();
  MinuitRes->mnexcm("HESSE", arglistRes,0,ierflgRes);
  MinuitRes->mnexcm("MINOs 10000000 3 3", arglistRes,0,ierflgRes);
  // Print results
  Double_t aminRes,edmRes,errdefRes;
  Int_t nvparRes,nparxRes,icstatRes;
  Double_t * parRes= new Double_t [numpar] ;
  Double_t*       parErrRes= new Double_t [numpar] ;    
  for(Int_t i=0;i<numpar;i++)
  {
    MinuitRes->GetParameter(i,parRes[i],parErrRes[i]);
    //cout << i << " " << parRes[i] << endl;
  }
  Int_t npar=numpar;
  Double_t grad=0;
  Double_t fval=0;
  Int_t flag=1;
  MinuitRes->Eval(npar, &grad, fval, parRes, flag);
  cout << "fval = " << fval << endl;

  MinuitRes->mnstat(aminRes,edmRes,errdefRes,nvparRes,nparxRes,icstatRes);
  MinuitRes->mnprin(numpar,aminRes);

  int nf=MinuitRes->GetNumFreePars();
  if(FREE_ENERGY_FIT) nf-=NEp;

  cout.precision(15);
  cout<<"Minuit Mass= "<<PDGMASS+parRes[2]*2.<<endl;
  cout<<"PDG Mass= "<<PDGMASS<<endl;
  cout.precision(4);
  cout<<"M-Mpdg="<<parRes[2]*2.<< " +- " << parErrRes[2]*2. <<  " MeV." << endl;
  cout<< "chi2/ndf = " <<MinChi2 << "/(" << NpPP<<"-"<<nf<<") = "  << MinChi2/(NpPP-nf) << ", P(chi2)=" << TMath::Prob(MinChi2,NpPP-nf) << endl;
  cout.precision(15);
  cout << "Contribution to chi square:" << endl;
  cout.precision(4);
  cout << "chi2 signal: " << CHI2_SIGNAL <<  " or " << CHI2_SIGNAL/CHI2_TOTAL*100 << "%" << endl;
  cout << "chi2 energy: " << CHI2_ENERGY <<  " or " << CHI2_ENERGY/CHI2_TOTAL*100 << "%" <<  endl;
  cout << "chi2 lum: "    << CHI2_LUM    <<  " or " << CHI2_LUM/CHI2_TOTAL*100 << "%" << endl;
  cout << "Total chi2: " << CHI2_TOTAL << endl;
  Double_t* parPsiPF    = new Double_t [idRNP];
  Double_t* parPsiPF2    = new Double_t [idRNP];
  parPsiPF[idRbg]=parRes[0];
  parPsiPF[idRM]=parRes[2];//parPsiP[Iscan][ippeff];   
  parPsiPF[idRSw]=parRes[3];   
  parPsiPF[idReff]=parRes[1];
  parPsiPF[idRFreeGee]=0;
  parPsiPF[idRTauEff]=0;


  /* 
   *  Fit both scan togegher
   *  This code must be refubrished to use with boost program options.
  if(arguments.scan==3)
  {
    int nep = FREE_ENERGY_FIT ? NEp : 0;
    if(BOTH_FIT==1)
    {
      parPsiPF2[idRbg]=parRes[0+4+nep];
      parPsiPF2[idReff]=parRes[1+4+nep];
      parPsiPF2[idRM]=parRes[2];
      parPsiPF2[idRSw]=parRes[3];   
      parPsiPF2[idRFreeGee]=0;
      parPsiPF2[idRTauEff]=0;
      for(unsigned i=0; i<idRNP;i++)
      {
        cout << i << " " << parPsiPF2[i] << endl;
      }
    }
    if(BOTH_FIT==2)
    {
      parPsiPF2[idRbg]=parRes[0+4+nep];
      parPsiPF2[idReff]=parRes[1+4+nep];
      parPsiPF2[idRM]=parRes[2];
      parPsiPF2[idRSw]=parRes[2+4+nep];   
      parPsiPF2[idRFreeGee]=0;
      parPsiPF2[idRTauEff]=0;
      for(unsigned i=0; i<idRNP;i++)
      {
        cout << i << " " << parPsiPF2[i] << endl;
      }
    }
  }
  */


  double EnergyChi2=0;
	TGraphErrors * dEgr = new TGraphErrors;//energy deviation graph
	TGraphErrors * dNgr = new TGraphErrors;//signal discrepancy graph
	gStyle->SetOptFit(kTRUE);
  for(int is=0;is<NEp;is++)
  {       
    WErrInScan[is]=WErrInScan[is];
    //if(arguments.FreeEnergy==1) WInScan[is]+=2.*+parRes[is+4];
    EnergyChi2+= sq(parRes[is+4]/EErrInScan[is]);
		dEgr->SetPoint(is, WInScan[is], parRes[is+4]*2);
		dEgr->SetPointError(is, 0,  WErrInScan[is]);
		dNgr->SetPoint(is, WInScan[is], SignalDiscrepancy[is]);
		dNgr->SetPointError(is, WErrInScan[is], SignalDiscrepancyError[is]);
  }
	TCanvas * dEc = new TCanvas("dEc", "Energy deviation");
	dEgr->SetMarkerStyle(21);
	dEgr->SetLineWidth(2);
	dEgr->SetMarkerSize(1);
	dEgr->Draw("ap");
	dEgr->GetXaxis()->SetTitle("W, MeV");
	dEgr->GetYaxis()->SetTitle("W_{exp}-W_{CBS}, MeV");
	dEgr->Fit("pol0", "Q");
	TCanvas * dNc = new TCanvas("dNc", "Number of visible events of signal deviation");
	dNgr->SetMarkerStyle(21);
	dNgr->SetLineWidth(2);
	dNgr->SetMarkerSize(1);
	dNgr->Draw("ap");
	dNgr->GetXaxis()->SetTitle("W, MeV");
	dNgr->GetYaxis()->SetTitle("N_{vis} - N_{exp}");
	dNgr->Fit("pol0", "Q");

  GrRes=new TGraphErrors(NEp,WInScan,CrossSInScan,WErrInScan,CrossSErrInScan);
  TF1 * fitfun1=0;
  TF1 * fitfun2=0;
  switch(RESONANCE)
  {
    case JPSIRES:
      fitfun1  = new TF1("FitJPsi",FCrSJpsiAzimov,1540*ScaleEGr,1560*ScaleEGr,idRNP);
      break;
    case PSI2SRES:
      fitfun1  = new TF1("FitPsiP",FCrSPPrimeAzimov,1836.*ScaleEGr,1855*ScaleEGr,idRNP);;
      break;
  }
  TF1* FitPsiP = fitfun1;
  TF1* FitPsiP2=new TF1("FitPsiP2",FCrSPPrimeAzimov,1836.*ScaleEGr,1855*ScaleEGr,idRNP);  
  FitPsiP->SetParameters( parPsiPF);
  FitPsiP2->SetParameters( parPsiPF2);
  TCanvas* TestCanv=new TCanvas("BES_PSIP_SCAN","BES PsiP Scan",900,700); 
  TestCanv->SetGridx();
  TestCanv->SetGridy();
  TestCanv->SetFillColor(0);
  TestCanv->SetBorderMode(0);
  //TestCanv->cd();
  GrRes->Draw("AP");
  gPad->SetBorderMode(0);
  GrRes->GetXaxis()->SetTitle("W, MeV");
  GrRes->GetYaxis()->SetTitle("#sigma, nb");
  double xx=0;
  switch(RESONANCE)
  {
    case JPSIRES:  
      GrRes->SetTitle("J/#psi scan");
      xx=1544*ScaleEGr;
      break;
    case PSI2SRES: 
      GrRes->SetTitle("#psi(2S) scan");
      xx=1838*ScaleEGr;
      break;
  };

  FitPsiP->Draw("SAME");
  if(BOTH_FIT==1 || BOTH_FIT==2)
  {
    FitPsiP2->SetLineColor(kBlue);
    FitPsiP2->Draw("SAME");
    TLegend * l= new TLegend(0.8,0.9,1.0,1.0);
    l->AddEntry(FitPsiP,"First scan", "l");
    l->AddEntry(FitPsiP2,"Second scan", "l");
    l->Draw();
  }
  char Info1[100];
  TLatex*  latexM1=new TLatex();
  latexM1->SetTextSize(0.038);
  latexM1->SetTextColor(2);       
  sprintf(Info1,"#chi^{2} = %3.3f/ (%d -%d) =%3.3f",MinChi2,NpPP,nf,MinChi2/(NpPP-nf)); 
  //double yy = parRes[0]*5;
  double yy = FitPsiP->GetMaximum();
  latexM1->DrawLatex(xx,yy,Info1);
  sprintf(Info1,"#Delta M = %3.3f#pm%3.3f MeV",parRes[2]*2.,parErrRes[2]*2.);
  latexM1->DrawLatex(xx,yy*0.8,Info1);
  TLatex * latexSw = new TLatex();
  latexSw->SetTextSize(0.038);
  latexSw->SetTextColor(2);
  sprintf(Info1,"#sigma_{W} = %1.3f #pm %1.3f MeV",parRes[3],parErrRes[3]); 
  latexSw->DrawLatex(xx,yy*0.6,Info1);

  TLatex * latexProb = new TLatex();
  latexProb->SetTextSize(0.038);
  latexProb->SetTextColor(2);
  sprintf(Info1,"P(#chi^{2}) = %1.4f",TMath::Prob(MinChi2,NpPP-nf)); 
  latexProb->DrawLatex(xx,yy*0.4,Info1);

  TestCanv->Update();
  delete [] En_;   
  delete [] Eerr_;
  delete [] Nmh_;
  delete [] Nbb_;                   
  delete [] Ngg_;                   
  delete [] En;   
  delete [] Eerr;
  delete [] Nmh;
  delete [] Nbb;               
  delete [] Ngg;               
  delete FitRes;
  theApp->Run();
  return 0;
}

static int FCNcall=0;

int FCallCheck=0;

Double_t parprev[10]={0,0,0, 0,0,0, 0,0,0, 0};
void fcnResMult(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //calculate chisquare
  Double_t chisq =0;
  Double_t chisqbb = 0;
  Double_t chisqgg = 0;
  Double_t chisqmh =0 ;
  Double_t sigmaFull;
  Double_t sigmaMH;
  Double_t sigmaBB;
  Double_t sigmaGG;
  Double_t nFull;
  Double_t Nexp; //expected number of luminosity events
  Double_t lumFull;
  Double_t Energy;
  Double_t parmh[idRNP];
  Double_t EnergyChi2=0;
  Double_t SigmaWChi2=0;
  Double_t LumChi2=0;
  Double_t SignalChi2=0;
  parmh[idRbg]=par[0];
  parmh[idReff]=par[1];
  parmh[idRM]=par[2];
  parmh[idRSw]=par[3];   
  parmh[idRFreeGee]=0;       
  parmh[idRTauEff]=0;       
  if(iflag==1)
  {
    cout << setw(4) << "Parnumber" <<  setw(20) << "value" << endl;
    for(int i=0;i<npar;i++)
    {
      cout << setw(4)<< i << setw(20) << par[i] << endl;
    }
    cout << setw(4)<< "point#" << setw(15)<< "energy, MeV" <<  setw(15) << "energy chi2" << setw(15) <<  "signal chi2" << endl;
  }
  if( FCNcall ==0 )
  {
    std::cout << std::setw(3)  << "#pn";
    std::cout << std::setw(10) << "E, MeV";
    std::cout << std::setw(10)  << "Nmh";
    std::cout << std::setw(10)  << "Nlum";
    std::cout << std::setw(10) << "sigBB,nb";
    std::cout << std::setw(10) << "sigMH,nb";
    std::cout << std::setw(10) << "Lmh,1/nb";
    std::cout << std::setw(10) << "Lbb,1/nb";
    std::cout << std::setw(10) << "L,1/nb";
    for(unsigned i=0;i<4;++i) std::cout << setw(7) << "p"+std::string(boost::lexical_cast<std::string>(i));
    std::cout << std::endl;

  }
  for (Int_t i=0;i<NumEpoints;i++)
  {
    Energy=EInScan[i];
    double echi2=0;
    if(FREE_ENERGY_FIT)
    { 
      Energy+=par[4+i];
      echi2=(par[4+i]*par[4+i]/(EErrInScan[i]*EErrInScan[i]));
      EnergyChi2+=echi2;
    }
    if(USE_CBS_SIGMAW)
    {
      //parmh[idRSw]=SigmaWInScan[i];
      SigmaWChi2+= sq((parmh[idRSw] - SigmaWInScan[i])/dSigmaWInScan[i]);
    }
    //Correction to efficiency
    //parmh[idReff]=par[1]*MhadrCor(2*Energy);
    //cout << "Correction to efficiency: " << (MhadrCor(2*Energy)-1)*1000. << " ppm" << endl;
    switch(RESONANCE)
    {
      case JPSIRES:
        sigmaMH=CrSOniumR(_MethodAzimov,_IdJPsi,Energy,parmh);
        break;
      case PSI2SRES:
        sigmaMH=CrSOniumR(_MethodAzimov,_IdPsiPrime,Energy,parmh);
        break;
    };
    sigmaBB=CrossSBhabhaPP(Energy,&CrossBhabha);                        
    sigmaGG=CrossSBhabhaPP(Energy,&CrossGG);                        
    switch(LUMINOSITY)
    {
      case NMHEELUM:
        nFull = (NbbInScan[i]+NmhInScan[i]);
        sigmaFull=sigmaMH+sigmaBB;
        lumFull = nFull/sigmaFull;
        Nexp = sigmaBB*lumFull;
        nFull=1e100;
        break;
      case NMHGGLUM:
        nFull = (NggInScan[i]+NmhInScan[i]);
        sigmaFull=sigmaMH+sigmaGG;
        lumFull = nFull/sigmaFull;
        Nexp = sigmaGG*lumFull;
        nFull=1e100;
        break;
      case NMHEEGGLUM:
        nFull = (NbbInScan[i]+NggInScan[i]+NmhInScan[i]);
        sigmaFull=sigmaMH+sigmaBB+sigmaGG;
        lumFull = nFull/sigmaFull;
        Nexp = (sigmaGG+sigmaBB)*lumFull;
        nFull=1e100;
        break;
      default:
        nFull = NLum[i];
        //lumFull=LumInScan[i]*(LUM_COR ? BBLumCorr[i] : 1);
        //lumFull=LumInScan[i]*(LUM_COR ? BBIntCor(Energy*2, par[3]) : 1);
        lumFull=LumInScan[i];
        //cout << "W=" << Energy*2 << " cor =" << BBIntCor(Energy*2) << endl;
        Nexp=NLum[i];
        break;
    }
    if( FCNcall ==0 )
    {
      std::cout << setw(3) << i+1;
      std::cout << boost::format("%10.3f") % Energy;
      std::cout << setw(10) << NmhInScan[i];
      std::cout << setw(10) << NbbInScan[i];
      std::cout << boost::format("%10.3f") % sigmaBB;
      std::cout << boost::format("%10.3f") % sigmaMH;
      std::cout << boost::format("%10.3f") % (NmhInScan[i]/sigmaMH);
      std::cout << boost::format("%10.3f") % (NbbInScan[i]/sigmaBB);
      std::cout << boost::format("%10.3f") % lumFull;
      for(unsigned i=0;i<4;++i) std::cout << boost::format("%7.2f") % parmh[i];
      std::cout << std::endl;
      NpPP++;
    }
    chisqmh=0;
    //calculate signal contribution
    if(NmhInScan[i]>0)
    {
      double N=sigmaMH*lumFull; //expected number of events
      if(USE_CHI2)
      {
        chisqmh = sq(N - NmhInScan[i])/(NmhInScan[i]*(1. + NmhInScan[i]/nFull)); //just chi square
      } 
      else chisqmh = 2*(NmhInScan[i]*log(NmhInScan[i]/N) +N - NmhInScan[i]); //likelihood for low statistics
			SignalDiscrepancy[i] = NmhInScan[i]-N;
			SignalDiscrepancyError[i] = sqrt(NmhInScan[i]*(1. + NmhInScan[i]/nFull));
    }
    else if(NmhInScan[i]==0) 
    {
        cout << "ERROR = " << endl;
        chisqmh =  sigmaMH*lumFull;
    }
    SignalChi2+=chisqmh;
    chisq+=chisqmh; 

    chisqbb=0;
    if(NLum[i]!=Nexp)
      if(NLum[i]>0)
      {
        double Nobs = NLum[i];
        //cout << "LUMIN="<<LUMINOSITY << endl;
        if(USE_CHI2) chisqbb = sq(Nexp - Nobs)/Nobs; //just chi square
        else chisqbb = 2*(Nobs*log(Nobs/Nexp)+ Nexp - Nobs); //likelihood for low statistics
      }
      else
      {
        cout << "ERROR = " << endl;
        chisqbb =  sigmaBB*lumFull;
      }
    chisq+=chisqbb;
    LumChi2 += chisqbb;

    if(iflag==1)
    {
      cout << setw(4)<< i << setw(15)<< EInScan[i] << setw(15) << echi2 << setw(15) <<  chisqmh << endl;
    }
  }
  f = chisq;
  if(FREE_ENERGY_FIT) f+=EnergyChi2;//;    
  if(USE_CBS_SIGMAW) f+=SigmaWChi2;    
  if(MinChi2>f)
  {
    MinChi2=f;
    CHI2_ENERGY = EnergyChi2;
    CHI2_TOTAL=f;
    CHI2_LUM = LumChi2;
    CHI2_SIGNAL = SignalChi2;
  }
  FCNcall=1;
}

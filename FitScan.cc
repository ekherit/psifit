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
#include <regex>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string_regex.hpp>
#include <boost/algorithm/string/classification.hpp>

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
#include <TRandom.h>

R__EXTERN TSystem *gSystem;
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("ROOTPROGRAM","ROOTPROGRAM", initfuncs);
using namespace std;

#include "FitOniumRCompactLib.hh"

#include "utils.h"

#define  NumMaxP 120
Int_t    NpPP=0;
Double_t CrossBhabha=0;
Double_t CrossGG=0;
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
void print_result(std::ostream & os, TMinuit * minuit, std::string sharp="#");

bool USE_CHI2 = true;
bool USE_CBS_SIGMAW = false; 
bool LUM_COR = true;
bool FREE_ENERGY_FIT=true;
unsigned BOTH_FIT=0;
double DEDIFF=0.02;
double EMS_SCALE=1;

int RANDOM_SEED=0;
double ENERGY_VARIATION=0;
std::string INPUT_FILE = "scan.txt";
std::string OUTPUT_FILE = "fitscan.txt";
std::string RESULT_FILE = "result.txt";
std::vector <double> PAR_INI; //initial parameter value
std::string PAR_INI_STRING;

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
  AUTORES=0,
  JPSIRES=1,
  PSI2SRES=2,
};

unsigned RESONANCE;


double CHI2_TOTAL;
double CHI2_LUM;
double CHI2_ENERGY;
double CHI2_SIGNAL;

double E_CROSS_NORM; //Beam Energy used for normalization of luminosity calculation

void draw_signal_and_energy_deviation(int NEp, double * parRes);

int main(int argc, char **argv)
{
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  opt_desc.add_options()
    ("help,h","Print this help")
    ("verbose,v","Verbose output")
    ("input",po::value<std::string>(&INPUT_FILE)->default_value("scan.txt"),"data file")
    ("lum,l", po::value<std::string>()->default_value("gg"),"luminosity used: bes, ee, gg (default), eegg, mhee, mhgg, mheegg")
    ("resonance,R", po::value<unsigned>(&RESONANCE)->default_value(AUTORES),"resonance: 0 - AUTO, 1 - JPSI, 2 - PSI2S")
    ("free-energy,E","allow energy points to be variated (default)")
    ("nofree-energy","disable energy points to be variated")
    ("both", po::value<unsigned>(&BOTH_FIT)->default_value(0),"Fit mode 0,1,2")
    ("lum-cor","Use MC luminosity correction")
    ("nochi2", "Not use chi2")
    ("cross-section-ee",po::value<double>(&CrossBhabha), "Bhabha cross section, nb")
    ("cross-section-gg", po::value<double>(&CrossGG), "Gamma-gamma cross section, nb")
    ("dE",po::value<double>(&DEDIFF)->default_value(0.02),"Combine run with energy in dE interval, MeV")
    ("skip",po::value< std::vector<unsigned> >(),"list of skipped points")
    ("ems-error",po::value< double >(),"ems energy measurement error for each point")
    ("ems-scale", po::value <double>(&EMS_SCALE)->default_value(1), "Scale ems energy")
    ("variate-energy",po::value<double>(&ENERGY_VARIATION), "Variate point energy")
    ("seed",po::value<int>(&RANDOM_SEED), "Random seed")
    ("output,-o", po::value<std::string>(&OUTPUT_FILE), "Output file with the result")
    ("result,-r", po::value<std::string>(&RESULT_FILE), "Result of the fit accumulated in this file")
    ("exit", "exit after fitting")
    ("par", po::value<std::string> (&PAR_INI_STRING), "Initial parameter values")
    ;
  po::positional_options_description pos;
  pos.add("input",-1);
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


  if(opt.count("nochi2")) USE_CHI2=false;
  std::cout << "Chi2 fit: " << boolalpha << USE_CHI2 << std::endl;

  LUM_COR = opt.count("lum-cor") && LUMINOSITY==NEELUM;
  std::cout << "Luminosity corrections: " << LUM_COR << std::endl;

  if(opt.count("free-energy")) FREE_ENERGY_FIT=true;
  if(opt.count("nofree-energy")) FREE_ENERGY_FIT=false;
  std::cout << "Use free energy fit: " << FREE_ENERGY_FIT << std::endl;

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
  std::cout << "Reading data from file " << data_file_name << "... ";
  npMHFile=GetNumRows(data_file_name.c_str(),dimMHFile);       
  FillArrayFromFile(data_file_name.c_str(),AllMH,dimMHFile,npMHFile);  
  std::cout << "read " << npMHFile << " points." << std::endl;
  if(opt.count("variate-energy"))
  {
    cout << "Variate energy in each point on: " << ENERGY_VARIATION << " MeV" << endl;
    if(RANDOM_SEED==0)
    {
      RANDOM_SEED = time(0);
      cout << "Default random seed (0)  then use current time" << endl;
    }
    cout << "Random seed is " << RANDOM_SEED << endl;
    TRandom r(RANDOM_SEED);
    cout << setw(5) << "point" << setw(20) << "variation, MeV"  << setw(20) << "new energy, MeV" << endl;
    for(int i=0;i<npMHFile;i++)
    {
      double EV = ENERGY_VARIATION == 0 ? AllMH[i][MHEnergyErr] : ENERGY_VARIATION;
      double dW = EV*r.Gaus();
      double W  = AllMH[i][MHEnergy]; 
      W = W + dW;
      AllMH[i][MHEnergy] = W;
      char buf[65535];
      sprintf(buf, "%5d%20.3f%20.3f", i+1, dW, W);
      cout << buf << endl;;
      //cout << setw(5) << i+1 << setw(15) << setprecision(4)<< dW << setw(15) << setprecision(7) <<  W << endl;
    }
  }

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
    AP[Aind][AEnergyErr]=AllMH[i][MHEnergyErr]*0.5*EMS_SCALE;       
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
  if(RESONANCE==AUTORES)
  {
    if(AP[0][AEnergy] < 1800) RESONANCE = JPSIRES;
    else RESONANCE=PSI2SRES;
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
  Double_t LG=0; //Integrated luminosity BES online
  Double_t Lee=0; //Integrated Bhabha luminosity
  Double_t Lgg=0; //integrated Gamma Gamma lum
  Double_t Leegg=0; //integrated BB+GG luminosity
  //fill global arrayes which needed to access from the FCN function
  double ems_error=0;
  if(opt.count("ems-error")) 
  {
    ems_error = opt["ems-error"].as<double>();
    cout << "Set all ems errors to : " << ems_error << endl;
  }

  E_CROSS_NORM = PDGMASS/2.;

  if(!opt.count("cross-section-gg")) CrossGG = estimate_cross_section2(E_CROSS_NORM, NEp, En, Le, Ngg);
  if(!opt.count("cross-section-ee")) CrossBhabha = estimate_cross_section2(E_CROSS_NORM, NEp, En, Le, Nbb);
  std::vector<double> Neegg(NEp); //sum of Nbb and Ngg
  for(int i=0;i<NEp;i++) Neegg[i]=Nbb[i]+Ngg[i];
  double CrossEEGG = estimate_cross_section2(E_CROSS_NORM, NEp, En, Le, &Neegg[0]);
  cout << "Estimate cross section for luminosity processes: " << endl;
  cout << setw(25) << "Gamma-gamma cross section"<< setw(20)  <<  CrossGG << " nb. " << endl;
  cout << setw(25) << "Bhabha cross section" << setw(20) << CrossBhabha << " nb" << endl;
	cout << "Luminosity used: ";
  std::string lumstr="gg";
  lumstr=opt["lum"].as<string>();
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
  for(int is=0;is<NEp;is++)
  {
    EInScan[is]=En[is];
    WInScan[is]=2.*EInScan[is];
    if(opt.count("ems-error")) EErrInScan[is]=ems_error;
    else EErrInScan[is]=Eerr[is];
    WErrInScan[is]=EErrInScan[is]*2;
    SigmaWInScan[is]=SW[is];
    dSigmaWInScan[is]=dSW[is];
    NmhInScan[is]=Nmh[is];   
    NbbInScan[is]=Nbb[is];       
    NggInScan[is]=Ngg[is];       
    //ECorrBB=1./CrossSBhabhaPP(En[is],&CrossBhabha);
    //ECorrGG=1./CrossSBhabhaPP(En[is],&CrossGG);
    //LumInScanEE[is]=NbbInScan[is]*ECorrBB;
    //LumInScanGG[is]=NggInScan[is]*ECorrGG;
    //LumInScanEEGG[is]=(NggInScan[is] + NbbInScan[is])/(1./ECorrBB + 1./ECorrGG);
    LumInScanEE[is]=NbbInScan[is]*sq(En[is]/E_CROSS_NORM)/CrossBhabha;
    LumInScanGG[is]=NggInScan[is]*sq(En[is]/E_CROSS_NORM)/CrossGG;
    LumInScanEEGG[is]=(NggInScan[is] + NbbInScan[is])*sq(En[is]/E_CROSS_NORM)/CrossEEGG;
    LG+=TMath::Max(Le[is],Lp[is]);
    Lee+= LumInScanEE[is];
    Lgg+= LumInScanGG[is];
    Leegg+=LumInScanEEGG[is];
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
  std::cout << "Total luminosity: Lbes=" << LG << "/nb, Lee=" << Lee << "/nb, Lgg="<<Lgg<<"/nb, Leegg="<< Leegg << "/nb" <<std::endl;
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

  //Double_t vstartRes[5]= {0,0.5,0,1.439,LUM_CROSS_SECTION};   
  Double_t vstartRes[5]= {0,0.5,0,1.0,LUM_CROSS_SECTION};   

  if(opt.count("par"))
  {
    cout << "Set initial parameter values: ";
    istringstream is(PAR_INI_STRING);
    PAR_INI.resize(0);
    while(!is.eof())
    {
      double tmp;
      is >> tmp;
      PAR_INI.push_back(tmp);
    }
    
    for(int i=0;i<5 && i<PAR_INI.size();i++)
    {
      cout << PAR_INI[i] << ", ";
      vstartRes[i]=PAR_INI[i];
    }
    cout << endl;
  }

  Double_t stepRes[5] =  {1, 0.1,0.1,0.1,0.0};



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
      MinuitRes->DefineParameter(j+4,NameP,0,0.1,-2.0,+2.0);        
    }
  }

  int nep = FREE_ENERGY_FIT==false ? 0 : NEp;

  MinuitRes->mnexcm("MIGRAD", arglistRes,numpar,ierflgRes);
  MinuitRes->mnimpr();
  MinuitRes->mnexcm("HESSE", arglistRes,0,ierflgRes);
  MinuitRes->mnexcm("MINOs 10000000 3 3", arglistRes,0,ierflgRes);
  // Print results
  Double_t aminRes,edmRes,errdefRes;
  Int_t nvparRes,nparxRes,icstatRes;
  Double_t * parRes= new Double_t [numpar] ;
  Double_t * parErrRes= new Double_t [numpar] ;    
  for(Int_t i=0;i<numpar;i++)
  {
    MinuitRes->GetParameter(i,parRes[i],parErrRes[i]);
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

  //cout.precision(15);
  //cout<<"Minuit Mass= "<<PDGMASS+parRes[2]*2.<<endl;
  //cout<<"PDG Mass= "<<PDGMASS<<endl;
  //cout.precision(4);
  //cout<<"M-Mpdg="<<parRes[2]*2.<< " +- " << parErrRes[2]*2. <<  " MeV." << endl;
  //cout<< "chi2/ndf = " <<MinChi2 << "/(" << NpPP<<"-"<<nf<<") = "  << MinChi2/(NpPP-nf) << ", P(chi2)=" << TMath::Prob(MinChi2,NpPP-nf) << endl;
  //cout.precision(15);
  //cout << "Contribution to chi square:" << endl;
  //cout.precision(4);
  //cout << "chi2 signal: " << CHI2_SIGNAL <<  " or " << CHI2_SIGNAL/CHI2_TOTAL*100 << "%" << endl;
  //cout << "chi2 energy: " << CHI2_ENERGY <<  " or " << CHI2_ENERGY/CHI2_TOTAL*100 << "%" <<  endl;
  //cout << "chi2 lum: "    << CHI2_LUM    <<  " or " << CHI2_LUM/CHI2_TOTAL*100 << "%" << endl;
  //cout << "Total chi2: " << CHI2_TOTAL << endl;
  

  print_result(std::cout, MinuitRes,"#");
  //copy input file
  std::ifstream input_file(INPUT_FILE, std::ios::binary);
  std::ofstream output_file(OUTPUT_FILE, std::ios::binary);
  output_file << input_file.rdbuf();
  print_result(output_file, MinuitRes,"#");

  //print result to input file
  Double_t* parPsiPF    = new Double_t [idRNP];
  Double_t* parPsiPF2    = new Double_t [idRNP];
  parPsiPF[idRbg]=parRes[0];
  parPsiPF[idRM]=parRes[2];//parPsiP[Iscan][ippeff];   
  parPsiPF[idRSw]=parRes[3];   
  parPsiPF[idReff]=parRes[1];
  parPsiPF[idRFreeGee]=0;
  parPsiPF[idRTauEff]=0;

  ofstream result_file(RESULT_FILE.c_str(), fstream::app);
  if(!result_file)
  {
    cerr << "Unable to open file " <<RESULT_FILE << endl;
  }
  else 
  {
    cout << "Print fit result to file: " << RESULT_FILE << endl;
    TMinuit * minuit = MinuitRes;
    char result_string[65535];
    double chi2= MinChi2/(NpPP-nf);
    double prob = TMath::Prob(MinChi2,NpPP-nf);
    sprintf(result_string,"%8.3f %5.3f   %5.3f %5.3f   %5.2f %5.2f  %4.2f  %4.2f  %4.2f  %10.5f",
        parRes[2]*2,     parErrRes[2]*2,
        parRes[3],       parErrRes[3],
        parRes[1]*100,   parErrRes[1]*100,
        fabs(parRes[0]),       parErrRes[0],
        chi2,
        prob
        );
    result_file << result_string << endl;
    result_file.close();
  }

  draw_signal_and_energy_deviation(NEp,parRes);


  GrRes=new TGraphErrors(NEp,WInScan,CrossSInScan,WErrInScan,CrossSErrInScan);
  GrRes->SetMarkerStyle(20);
  GrRes->SetMarkerSize(1);
  GrRes->SetMarkerColor(kBlack);
  GrRes->SetLineColor(kBlack);
  GrRes->SetLineWidth(2);
  TF1 * fitfun1=0;
  TF1 * fitfun2=0;
  double Emin = TMath::MinElement(GrRes->GetN(), GrRes->GetX())-1;
  double Emax = TMath::MaxElement(GrRes->GetN(), GrRes->GetX())+1;
  double Erng = TMath::Max(fabs(Emax-PDGMASS), fabs(Emin-PDGMASS));
  Emin = PDGMASS-Erng;
  Emax = PDGMASS+Erng;
  TCanvas* TestCanv=new TCanvas("BES_PSIP_SCAN","BES PsiP Scan",900,700); 
  TestCanv->SetGridx();
  TestCanv->SetGridy();
  TestCanv->SetFillColor(0);
  TestCanv->SetBorderMode(0);
  switch(RESONANCE)
  {
    case JPSIRES:
      fitfun1  = new TF1("FitJPsi",FCrSJpsiAzimov,Emin/2.*ScaleEGr,Emax/2.*ScaleEGr,idRNP);
      fitfun1->SetTitle("J/#psi");
      GrRes->SetTitle("J/#psi scan");
      TestCanv->SetTitle("J/psi");
      break;
    case PSI2SRES:
      fitfun1  = new TF1("FitPsiP",FCrSPPrimeAzimov,Emin/2.*ScaleEGr,Emax/2.*ScaleEGr,idRNP);
      fitfun1->SetTitle("#psi(2S)");
      GrRes->SetTitle("#psi(2S) scan");
      TestCanv->SetTitle("psi(2S)");
      break;
  }
  fitfun1->SetLineColor(kRed);
  TF1* FitPsiP = fitfun1;
  FitPsiP->SetParameters( parPsiPF);
  FitPsiP->Draw();
  FitPsiP->GetXaxis()->SetTitle("W, MeV");
  FitPsiP->GetYaxis()->SetTitle("#sigma_{obs}, nb");
  GrRes->Draw("p");
  gPad->SetBorderMode(0);
  GrRes->GetXaxis()->SetTitle("W, MeV");
  GrRes->GetYaxis()->SetTitle("#sigma, nb");

  //calculate positions
  double xx=(Emin+2)/2.*ScaleEGr;
  double yy = FitPsiP->GetMaximum()*0.9;

  char Info1[100];
  //Show chi2 / ndf
  TLatex*  latexM1=new TLatex();
  latexM1->SetTextSize(0.038);
  latexM1->SetTextColor(kBlack);       
  sprintf(Info1,"#chi^{2}/ndf = %3.2f / (%d -%d) =%4.2f",MinChi2,NpPP,nf,MinChi2/(NpPP-nf)); 
  latexM1->DrawLatex(xx,yy,Info1);

  //Show mass shift
  sprintf(Info1,"M - M_{PDG} = %3.3f #pm %3.3f MeV",parRes[2]*2.,parErrRes[2]*2.);
  latexM1->DrawLatex(xx,yy*0.8,Info1);

  //Show energy spread
  TLatex * latexSw = new TLatex();
  latexSw->SetTextSize(0.038);
  latexSw->SetTextColor(kBlack);
  sprintf(Info1,"#sigma_{W} = %1.3f #pm %1.3f MeV",parRes[3],parErrRes[3]); 
  latexSw->DrawLatex(xx,yy*0.6,Info1);

  //Show probability
  TLatex * latexProb = new TLatex();
  latexProb->SetTextSize(0.038);
  latexProb->SetTextColor(kBlack);
  sprintf(Info1,"P(#chi^{2}) = %3.1f%%",TMath::Prob(MinChi2,NpPP-nf)*100); 
  latexProb->DrawLatex(xx,yy*0.4,Info1);

  TestCanv->Update();
  char pdf_file[1024];
  char root_file[1024];
  sprintf(pdf_file, "%s.pdf", OUTPUT_FILE.c_str());
  sprintf(root_file, "%s.root", OUTPUT_FILE.c_str());
  gPad->SaveAs(pdf_file);
  gPad->SaveAs(root_file);
  
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
  if(!opt.count("exit")) theApp->Run();
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
    std::cout << std::setw(10) << "L,1/nb";
    std::cout << std::setw(10) << "Nmh";
    std::cout << std::setw(10) << "sigMH,nb";
    std::cout << std::setw(10) << "Lmh,1/nb";
    std::cout << std::setw(10) << "Nlum";
    std::cout << std::setw(10) << "sigLum,nb";
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
    //sigmaBB=CrossSBhabhaPP(Energy,&CrossBhabha);                        
    //sigmaGG=CrossSBhabhaPP(Energy,&CrossGG);                        
    sigmaBB = CrossBhabha * sq(E_CROSS_NORM/Energy);                     
    sigmaGG = CrossGG * sq(E_CROSS_NORM/Energy);                     
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
      std::cout << boost::format("%10.3f") % lumFull;
      std::cout << setw(10) << NmhInScan[i];
      std::cout << boost::format("%10.3f") % sigmaMH;
      std::cout << boost::format("%10.3f") % (NmhInScan[i]/sigmaMH);
      std::cout << setw(10) << NLum[i];
      std::cout << boost::format("%10.3f") % (NLum[i]/LumInScan[i]);
      std::cout << boost::format("%10.3f") % LumInScan[i];
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


void draw_signal_and_energy_deviation(int NEp, double * parRes)
{
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

  cout << "Draw energy deviation" << endl;
	TCanvas * dEc = new TCanvas("dEc", "Energy deviation", 640,480);
	dEgr->SetMarkerStyle(21);
	dEgr->SetLineWidth(2);
	dEgr->SetMarkerSize(1);
	dEgr->Draw("ap");
	dEgr->GetXaxis()->SetTitle("W, MeV");
	dEgr->GetYaxis()->SetTitle("W_{exp}-W_{CBS}, MeV");
	dEgr->Fit("pol0", "Q");
  cout << "Draw signal deviation" << endl;
	TCanvas * dNc = new TCanvas("dNc", "Number of visible events of signal deviation",640,480);
	dNgr->SetMarkerStyle(21);
	dNgr->SetLineWidth(2);
	dNgr->SetMarkerSize(1);
	dNgr->Draw("ap");
	dNgr->GetXaxis()->SetTitle("W, MeV");
	dNgr->GetYaxis()->SetTitle("N_{vis} - N_{exp}");
	dNgr->Fit("pol0", "Q");
}



void print_result(std::ostream & os, TMinuit * minuit, string sharp)
{
  size_t N=minuit->GetNumPars();
  std::vector<double> P(N); //paramers
  std::vector<double> dP(N); //paramers errors
  for(int i=0;i<N;i++)
  {
    minuit->GetParameter(i, P[i], dP[i]);
  }
  //double grad=0;
  //double fval=0;
  //int flag=1;
  //minuit->Eval(N, &grad, fval, &P[0], flag);
  int nf = minuit->GetNumFreePars();
  if(FREE_ENERGY_FIT) nf-=NumEpoints;
  int ndf = (NpPP-nf);
  double chi2= MinChi2/ndf;
  double prob = TMath::Prob(MinChi2,NpPP-nf);
  typedef boost::format bf;
  os << sharp << " Contribution to chi square:" << endl;
  os << sharp << "   chi2 signal: " << CHI2_SIGNAL <<  " or " << CHI2_SIGNAL/CHI2_TOTAL*100 << "%" << endl;
  os << sharp << "   chi2 energy: " << CHI2_ENERGY <<  " or " << CHI2_ENERGY/CHI2_TOTAL*100 << "%" <<  endl;
  os << sharp << "      chi2 lum: " << CHI2_LUM    <<  " or " << CHI2_LUM/CHI2_TOTAL*100 << "%" << endl;
  os << sharp << "    Total chi2: " << CHI2_TOTAL   << endl;
  os << sharp << " ndf    = " << ndf  << endl;
  os << sharp << " chi2/ndf = " <<MinChi2 << "/(" << NpPP<<"-"<<nf<<") = "  << chi2 <<  endl;
  os << sharp << " P(chi2) =" << TMath::Prob(MinChi2,NpPP-nf) << endl;
  os << sharp << " Minuit Mass = "<< bf("%8.3f") % (PDGMASS+P[2]*2.)<<endl;
  os << sharp << " PDG Mass = "   << bf("%8.3f") % PDGMASS <<  endl;
  os << sharp << " " << bf("M-Mpdg = %8.3f +- %5.3f MeV") % (P[2]*2.) %  (dP[2]*2.) << endl;
  os << sharp << " " << bf("    SW = %5.3f +- %5.3f MeV") % P[3] %  dP[3] << endl;
  os << sharp << " " << bf("   eps = %5.2f +- %4.2f %%") % (P[1]*100) %  (dP[1]*100) << endl;
  os << sharp << " " << bf("    bg = %5.2f +- %5.2f nb") % P[0] %  dP[0] << endl;
}

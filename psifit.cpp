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
#include <TFile.h>

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

std::vector<double> CrossSInScan;
std::vector<double> CrossSErrInScan;
std::vector<double> CrossSBBInScan;
std::vector<double> CrossSBBErrInScan;
std::vector<double> CrossSGGInScan;
std::vector<double> CrossSGGErrInScan;
std::vector<double> LumInScan, LumInScanError;
std::vector<double> LumInScanGG;
std::vector<double> LumInScanEE;
std::vector<double> LumInScanEEGG;
std::vector<double> NmhInScan;
std::vector<double> LumLgammaInScan;
std::vector<double> NbbInScan;
std::vector<double> NggInScan;
std::vector<double> NLum;
std::vector<double> EInScan;
std::vector<double> SigmaWInScan;
std::vector<double> dSigmaWInScan;
std::vector<double> WInScan;
std::vector<double> EErrInScan;
std::vector<double> WErrInScan;
std::vector<double> SignalDiscrepancy;
std::vector<double> SignalDiscrepancyError;
std::vector<double> BBLumCorr;

struct ScanPoint_t 
{
  double E,dE; //beam energy, MeV
  double W,dW; //cm energy
  double Sw,dSw; //energy spread
  double L;//online lum
  double Lgg; //gamma-gamma luminosity
  double Lee; //Bhabha luminosity
  double Leegg; //luminosity based on Lee and Lgg
  double XSmh, dXSmh; //cross section of the signal
  double XSee, dXSee; //cross section of bhabha
  double XSgg, dXSgg; //cross section of gamma-gamma
  Long64_t Nmh; //number of multihadronic (signal)
  Long64_t Nee; //bhabha
  Long64_t Ngg; //gamma-gamma
  Long64_t Nlum;//number of events of luminosity process (could be = Ngg or Nee, or Ngg+Nee)
  double Lcor; //поправки к светимости
  double DNmh, dDNmh; //разница между фитом и данными в сисле событий
  double DE, dDE; //разница между фитом и данными для энергии
};

void set_number_of_ponints(int NEp)
{
  CrossSInScan.resize(NEp);
  CrossSErrInScan.resize(NEp);
  CrossSBBInScan.resize(NEp);
  CrossSBBErrInScan.resize(NEp);
  CrossSGGInScan.resize(NEp);
  CrossSGGErrInScan.resize(NEp);
  LumInScan.resize(NEp);
  LumInScanError.resize(NEp);
  LumInScanGG.resize(NEp);
  LumInScanEE.resize(NEp);
  LumInScanEEGG.resize(NEp);
  NmhInScan.resize(NEp);
  LumLgammaInScan.resize(NEp);
  NbbInScan.resize(NEp);
  NggInScan.resize(NEp);
  NLum.resize(NEp);
  EInScan.resize(NEp);
  SigmaWInScan.resize(NEp);
  dSigmaWInScan.resize(NEp);
  WInScan.resize(NEp);
  EErrInScan.resize(NEp);
  WErrInScan.resize(NEp);
  SignalDiscrepancy.resize(NEp);
  SignalDiscrepancyError.resize(NEp);
  BBLumCorr.resize(NEp);
}

Double_t MinChi2=1e+100;
void fcnResMult    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void print_result(std::ostream & os, TMinuit * minuit, std::string sharp="#");

bool USE_CHI2 = true;
bool USE_CBS_SIGMAW = false; 
bool LUM_COR = true;
bool FREE_ENERGY_FIT=true;
unsigned BOTH_FIT=0;
double ENERGY_BIN=0.02;
double EMS_SCALE=1;


int RANDOM_SEED=0;
double ENERGY_VARIATION=0;
std::string INPUT_FILE = "scan.txt";
std::string OUTPUT_FILE = "fitscan.txt";
std::string RESULT_FILE = "result.txt";
std::vector <double> PAR_INI; //initial parameter value
std::string PAR_INI_STRING;
std::map<std::string,int> PAR_INDEX;
int PAR_INDEX_DSIGMA=0;
double PAR_DM;

double AVERAGE_SIGMAW_CBS;
double AVERAGE_SIGMAW_CBS_ERROR;


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
double CHI2_CBS_SIGMAW;

double E_CROSS_NORM; //Beam Energy used for normalization of luminosity calculation

TF1 * get_result_function(const std::vector<double> & parRes, double Emin, double Emax);
TCanvas * draw_signal_and_energy_deviation(const std::vector<double> & parRes);


void print_cross_section(void)
{
  Double_t parmh[idRNP];
  parmh[idRbg]=0;
  parmh[idReff]=1;
  parmh[idRM]=0;
  parmh[idRSw]=1.5;   
  parmh[idRFreeGee]=0;       
  parmh[idRTauEff]=0;       
  vector<double> E = 
  {
    3675.9,
    3683.7,
    3685.1,
    3686.3,
    3687.6,
    3688.8,
    3693.5
  };
  cout << " Printout multihadronic crossection"<< endl;
  cout << " spread = " << parmh[idRSw] << " MeV" <<  endl;
  cout << " M  = " <<  _MPsiPrime << " MeV" << endl;
  cout << " Gtot  = " <<  _GtotPsiPrime << " MeV" << endl;
  cout << " Gee  = " <<  _GeePsiPrime << " MeV" << endl;
  cout << setw(15) << "Wcm, MeV" << setw(15) << " sigma, nb" << endl;
  for(int i=0;i < E.size() ;i++)
  {
    double sigmaMH=CrSOniumR(_MethodAzimov,_IdPsiPrime,E[i]/2.,parmh);
    cout << setw(15) << E[i] << setw(15) << sigmaMH << endl;
  }
}

std::string CFG_TITLE;

int main(int argc, char **argv)
{
  std::map<std::string, double> fixed_parameters;
  std::string fixed_parameters_string;
  std::string lumstr="gg";
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  opt_desc.add_options()
    ("help,h","Print this help")
    ("verbose,v","Verbose output")
    ("input",po::value<std::string>(&INPUT_FILE)->default_value("scan.txt"),"data file")
    ("output", po::value<std::string>(&OUTPUT_FILE)->default_value("fitscan.txt"), "Output file with the result")
    ("lum,l", po::value<std::string>(&lumstr)->default_value("gg"),"luminosity used: bes, ee, gg (default), eegg, mhee, mhgg, mheegg")
    ("resonance,R", po::value<unsigned>(&RESONANCE)->default_value(AUTORES),"resonance: 0 - AUTO, 1 - JPSI, 2 - PSI2S")
    ("free-energy,E","allow energy points to be variated (default)")
    ("nofree-energy","disable energy points to be variated")
    ("variate-energy",po::value<double>(&ENERGY_VARIATION), "Variate point energy")
    ("both", po::value<unsigned>(&BOTH_FIT)->default_value(0),"Fit mode 0,1,2")
    ("lum-cor","Use MC luminosity correction")
    ("nochi2", "Not use chi2")
    ("cross-section-ee",po::value<double>(&CrossBhabha), "Bhabha cross section, nb")
    ("cross-section-gg", po::value<double>(&CrossGG), "Gamma-gamma cross section, nb")
    ("dE",po::value<double>(&ENERGY_BIN)->default_value(0.0),"Combine run with energy in dE interval, MeV")
    ("skip",po::value< std::vector<unsigned> >(),"list of skipped points")
    ("ems-error",po::value< double >(),"ems energy measurement error for each point")
    ("ems-scale", po::value <double>(&EMS_SCALE)->default_value(1), "Scale ems energy")
    ("seed",po::value<int>(&RANDOM_SEED), "Random seed")
    ("result,-r", po::value<std::string>(&RESULT_FILE), "Result of the fit accumulated in this file")
    ("exit", "exit after fitting")
    ("par", po::value<std::string> (&PAR_INI_STRING), "Initial parameter values")
    ("par-dm", po::value<double>(&PAR_DM)->default_value(0), "initial par value for mass ")
    ("fix",po::value<std::string>(&fixed_parameters_string), "Fix paremeter: name=<value>[,name2=<value>]...") 
    ("print","Print cross section")
    ("cbs_sigmaw_each_point","Use different cbs energy spread for each energy point")
    ("title",po::value<std::string>(&CFG_TITLE)->default_value(""), "Add title to canvas with fit result")
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
  cout << "Before first vari energy" << endl;
  opt.count("variate-energy");
  cout << "after first vari energy" << endl;

  if(opt.count("help"))
  {
    std::clog << opt_desc;
    return 0;
  }

  if(opt.count("print"))
  {
    print_cross_section();
    return 0;
  }


  if(opt.count("nochi2")) USE_CHI2=false;
  std::cout << "Chi2 fit: " << boolalpha << USE_CHI2 << std::endl;

  LUM_COR = opt.count("lum-cor") && LUMINOSITY==NEELUM;
  std::cout << "Luminosity corrections: " << LUM_COR << std::endl;

  if(opt.count("free-energy")) FREE_ENERGY_FIT=true;
  if(opt.count("nofree-energy")) FREE_ENERGY_FIT=false;
  USE_CBS_SIGMAW_EACH_POINT  = opt.count("cbs_sigmaw_each_point");
  std::cout << "Use free energy fit: " << FREE_ENERGY_FIT << std::endl;



  TMinuit*  MinuitRes=0; 
  TF1* FitRes=0;
  TF1* FitResBG=0;
  TGraphErrors* GrRes=0;
  TLine* LineRes=0;
  int dimMHFile=11;
  int MHRun=0;
  int MHLum=1;
  int MHLumError=2;
  int MHEnergy=3;
  int MHEnergyErr=4;
  int MHSigmaW=5;
  int MHdSigmaW=6;
  int MHNmh=7; 
  int MHNee=8;
  int MHNgg=9; 
  int MHLumCor=10; 

  /* AllMH хранит массив читаемых из файла данных
   * первый индекс - номер строчки, второй индекс - номер столбца */
  std::vector< std::vector<double> > AllMH; 
  std::cout << "Reading data from file " << INPUT_FILE << "... ";
  try 
  {
    FillArrayFromFile(INPUT_FILE,dimMHFile, AllMH);  
  }
  catch(std::runtime_error e)
  {
    std::cerr << std::endl << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
  std::cout << "read " << AllMH.size() << " points." << std::endl;
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
    for(int i=0;i<AllMH.size();i++)
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

  int dimAP= 20;
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
  int ALumError = 19;
  std::vector < std::vector<double> > AP(AllMH.size());
  for(int i=0;i<AP.size();i++ )
  {      
    AP[i].resize(dimAP);
    //AP[Aind]=new double [dimAP];    
    AP[i][ARun]=AllMH[i][MHRun];                    
    AP[i][AEnergy]=AllMH[i][MHEnergy]*0.5;        
    AP[i][AEnergyErr]=AllMH[i][MHEnergyErr]*0.5*EMS_SCALE;       
    AP[i][ALe]=AllMH[i][MHLum];//Fill luminocity from bes online lum
    AP[i][ALp]=AllMH[i][MHLum];//This is the same for both electron and positron
    AP[i][AMHEv]=AllMH[i][MHNmh];
    AP[i][ASigmaW] = AllMH[i][MHSigmaW];
    AP[i][AdSigmaW] = AllMH[i][MHdSigmaW];
    AP[i][AEE]=AllMH[i][MHNee];
    AP[i][AGG]=AllMH[i][MHNgg];
    AP[i][ALumCor]=AllMH[i][MHLumCor];
    AP[i][ALumError]=AllMH[i][MHLumError];
  }
  if(RESONANCE==AUTORES)
  {
    if(AP[0][AEnergy] < 1800) RESONANCE = JPSIRES;
    else RESONANCE=PSI2SRES;
  }
  std::map<std::string,double> initial_par_value;
  initial_par_value["dm"]=opt["par-dm"].as<double>(); //half mass shift
  initial_par_value["bg"]=10;//nb //background (continuum)
  initial_par_value["eps"]=0.6; //efficiency
  initial_par_value["sigma"]=1.3; //energy spread
  std::cout << "Fit resonance: ";
  switch(RESONANCE)
  {
    case JPSIRES:
      PDGMASS=_MJPsi;
//      initial_par_value["dm"]=0.112; //half mass shift
      initial_par_value["dm"]=PAR_DM; //half mass shift
      initial_par_value["sigma"]=1.1; //half mass shift
      std::cout << "JPSI";
      break;
    case PSI2SRES:
      PDGMASS=_MPsiPrime;
      initial_par_value["dm"]=0.234; //half mass shift
      initial_par_value["sigma"]=1.44; //half mass shift
      std::cout << "PSI2S";
      break;
  }
  std::cout << std::endl;
  std::vector<double> En(AP.size());
  std::vector<double> Eerr(AP.size());
  std::vector<double> SigmaW(AP.size());
  std::vector<double> dSigmaW(AP.size());
  std::vector<double> Nmh(AP.size());
  std::vector<double> Nbb(AP.size());
  std::vector<double> Ngg(AP.size());
  std::vector<double> Le(AP.size());
  std::vector<double> Lp(AP.size());
  std::vector<double> Lcor(AP.size());
  std::vector<double> Lerror(AP.size());
  for(int i=0;i<En.size();i++)
  {   
    En[i]=AP[i][AEnergy];
    Eerr[i]=AP[i][AEnergyErr];
    Nmh[i]=AP[i][AMHEv];
    Nbb[i]=AP[i][AEE];  //here will be Nee or Ngg
    Ngg[i]=AP[i][AGG];  //here will be Nee or Ngg
    Le[i]=AP[i][ALe];   //here will be BES online lum
    Lp[i]=AP[i][ALp];  
    SigmaW[i]=AP[i][ASigmaW];
    dSigmaW[i]=AP[i][AdSigmaW];
		Lcor[i] = AP[i][ALumCor];
    Lerror[i] = AP[i][ALumError];
  }

  int numpar=4;
  if(FREE_ENERGY_FIT) numpar+=En.size();
  if(USE_CBS_SIGMAW_EACH_POINT) numpar+=En.size();
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

  if(!opt.count("cross-section-gg")) CrossGG = estimate_cross_section2(E_CROSS_NORM,  En, Le, Ngg);
  if(!opt.count("cross-section-ee")) CrossBhabha = estimate_cross_section2(E_CROSS_NORM,  En, Le, Nbb);
  std::vector<double> Neegg(En.size()); //sum of Nbb and Ngg
  for(int i=0;i<Neegg.size();i++) Neegg[i]=Nbb[i]+Ngg[i];
  double CrossEEGG = estimate_cross_section2(E_CROSS_NORM, En, Le, Neegg);
  cout << "Estimate cross section for luminosity processes: " << endl;
  cout << setw(25) << "Gamma-gamma cross section"<< setw(20)  <<  CrossGG << " nb. " << endl;
  cout << setw(25) << "Bhabha cross section" << setw(20) << CrossBhabha << " nb" << endl;
	cout << "Luminosity used: ";
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
  set_number_of_ponints(En.size());
  for(int is=0;is<EInScan.size();is++)
  {
    EInScan[is]=En[is];
    WInScan[is]=2.*EInScan[is];
    if(opt.count("ems-error")) EErrInScan[is]=ems_error;
    else EErrInScan[is]=Eerr[is];
    WErrInScan[is]=EErrInScan[is]*2;
    SigmaWInScan[is]=SigmaW[is];
    dSigmaWInScan[is]=dSigmaW[is];
    NmhInScan[is]=Nmh[is];   
    NbbInScan[is]=Nbb[is];       
    NggInScan[is]=Ngg[is];       
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
        LumInScanError[is]=Lerror[is];
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
    for(int is=0;is<En.size();is++)
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

  PAR_INDEX["bg"]    = 0;
  PAR_INDEX["eps"]   = 1;
  PAR_INDEX["dm"]    = 2;
  PAR_INDEX["sigma"] = 3;

  Double_t vstartRes[5]= {
    initial_par_value["bg"],
    initial_par_value["eps"],
    initial_par_value["dm"]/2.,
    initial_par_value["sigma"],
    LUM_CROSS_SECTION};   

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


  Double_t stepRes[5] =  {1, 0.1,1,0.1,0.0};

  MinuitRes->SetMaxIterations(100000000);                         


  //define parameters

  MinuitRes->DefineParameter(0,"bg",vstartRes[0],stepRes[0],-100,+100);
  MinuitRes->DefineParameter(1,"eff",vstartRes[1],stepRes[1],0.0,1.0);      
  MinuitRes->DefineParameter(2,"dM/2.",vstartRes[2],stepRes[2],-10,10);      
  MinuitRes->DefineParameter(3,"SigmaW",vstartRes[3],stepRes[3],0.5,2.0);
  std::cout << "Initial par value: " << std::endl;
  for(int i=0;i<5; i++)
  {
    cout << "par["<< i << "]=" << vstartRes[i] << "  step=" << stepRes[i] << std::endl;
  }

  //Fix parameter
  if(opt.count("fix"))
  {
    std::list<std::string> fixed_parameter_list;
    boost::split(fixed_parameter_list,fixed_parameters_string,boost::is_any_of(", "));
    for(auto item : fixed_parameter_list)
    {
      std::vector<std::string> name_value;
      boost::split(name_value,item,boost::is_any_of("= "));
      fixed_parameters[name_value[0]] = stod(name_value[1]);
    }

    for(auto fix : fixed_parameters)
    {
      std::string name = fix.first;
      double value = fix.second;
      int index = PAR_INDEX[name];
      MinuitRes->DefineParameter(index,name.c_str(),value,0,0,0);
      MinuitRes->FixParameter(index);
    }
  }




  int index=4;
  //if(USE_CBS_SIGMAW) MinuitRes->FixParameter(3);
  if(FREE_ENERGY_FIT)
  {
    PAR_INDEX["dE"]=index;
    for(int j=0;j<En.size();j++,index++)
    {
      char  NameP[10];
      sprintf(NameP,"dE%d",j);         
      MinuitRes->DefineParameter(index,NameP,0,0.1,-0.5,+0.5);        
//      if(j==5) MinuitRes->DefineParameter(j+4,NameP,0,0.1,0,+2.0);        
    }
  }

  if(USE_CBS_SIGMAW_EACH_POINT)
  {
    PAR_INDEX["dsigma"]=index;
    PAR_INDEX_DSIGMA = index;
    MinuitRes->FixParameter(PAR_INDEX["sigma"]);
    for(int j=0;j<En.size();j++,index++)
    {
      char  NameP[10];
      sprintf(NameP,"sW%d",j);         
      MinuitRes->DefineParameter(index,NameP,0,0.1,-0.5,+0.5);        
    }
  }
  MinuitRes->mnprin(3,0);

  //int nep = FREE_ENERGY_FIT==false ? 0 : EInScan.size();

  MinuitRes->mnexcm("MIGRAD", arglistRes,numpar,ierflgRes);
  MinuitRes->mnimpr();
  MinuitRes->mnexcm("HESSE", arglistRes,0,ierflgRes);
  MinuitRes->mnexcm("MINOs 10000000 3 3", arglistRes,0,ierflgRes);
  // Print results
  Double_t aminRes,edmRes,errdefRes;
  Int_t nvparRes,nparxRes,icstatRes;
  std::vector<double> parRes(numpar);
  std::vector<double> parErrRes(numpar);
  //Double_t * parRes= new Double_t [numpar] ;
  //Double_t * parErrRes= new Double_t [numpar] ;    
  for(Int_t i=0;i<numpar;i++)
  {
    MinuitRes->GetParameter(i,parRes[i],parErrRes[i]);
  }

  MinuitRes->mnstat(aminRes,edmRes,errdefRes,nvparRes,nparxRes,icstatRes);
  MinuitRes->mnprin(numpar,aminRes);

  int nf=MinuitRes->GetNumFreePars();
  if(FREE_ENERGY_FIT) nf-=EInScan.size();
  if(USE_CBS_SIGMAW_EACH_POINT) nf-=EInScan.size();

  print_result(std::cout, MinuitRes,"#");
  //copy input file
  std::ifstream input_file(INPUT_FILE, std::ios::binary);
  std::ofstream output_file(OUTPUT_FILE, std::ios::binary);
  output_file << input_file.rdbuf();
  print_result(output_file, MinuitRes,"#");


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
  TApplication* theApp=new  TApplication("App", &argc, argv);

  TFile f((OUTPUT_FILE+".root").c_str(),"RECREATE");
  TCanvas * diviation_c  = draw_signal_and_energy_deviation(parRes);


  GrRes=new TGraphErrors(EInScan.size(),&WInScan[0],&CrossSInScan[0],&WErrInScan[0],&CrossSErrInScan[0]);
  GrRes->SetMarkerStyle(20);
  GrRes->SetMarkerSize(1);
  GrRes->SetMarkerColor(kBlack);
  GrRes->SetLineColor(kBlack);
  GrRes->SetLineWidth(2);
  /* Вот таким изъебистым  путем я нахожу максимальное расстояние от
   * табличного значения края до границы по энергии для данных.
   * Это называется С++ головного мозга */
  double range = fabs(PDGMASS-*max_element ( WInScan.begin(), WInScan.end(),
      [](double x, double y){return fabs(x-PDGMASS) < fabs(y-PDGMASS);}
      ));
  //Зато картинка будет симметрична относительно табличного значения
  double Emin = PDGMASS-range-1;
  double Emax = PDGMASS+range+1;
  TF1 * FitPsiP = get_result_function(parRes, Emin, Emax);

  //Main canvas with fit result //
  TCanvas* mcanvas=new TCanvas("psi","Resonance",777,480); 
  mcanvas->SetFillColor(0);
  mcanvas->SetBorderMode(0);
  //mcanvas->SetGridx();
  //mcanvas->SetGridy();

  FitPsiP->Draw();
  FitPsiP->GetXaxis()->SetTitle("W, MeV");
  FitPsiP->GetYaxis()->SetTitle("#sigma_{obs}, nb");
  GrRes->Draw("p");
  gPad->SetBorderMode(0);
  GrRes->GetXaxis()->SetTitle("W, MeV");
  GrRes->GetYaxis()->SetTitle("#sigma_{obs}, nb");

  //calculate positions
  double xx=(Emin+1)/2.*ScaleEGr;
  double yy = FitPsiP->GetMaximum()*0.9;
  std::map<std::string, double> y;
  y["chi2"] = yy;
  y["M"]    = yy*0.8;
  y["Sw"]   = yy*0.6;
  y["<SW>"] = yy*0.4;
  y["P"]    = yy*0.2;

  typedef boost::format fmt;
  std::map<std::string, boost::format> format;
  format["chi2"] = fmt("#chi^{2}/ndf = %3.2f / (%d -%d) =%4.2f") % MinChi2  % NpPP % nf % (MinChi2/(NpPP-nf));
  format["M"]    = fmt("M - M_{PDG} = %3.3f #pm %3.3f MeV") % (parRes[2]*2.)  % (parErrRes[2]*2.);
  format["Sw"]   = fmt("#sigma_{W} = %1.3f #pm %1.3f MeV") % parRes[3] % parErrRes[3];
  format["P"]    = fmt("P(#chi^{2}) = %3.1f%%") % (TMath::Prob(MinChi2,NpPP-nf)*100);
  format["<SW>"] = fmt("<#sigma_{W}> = %1.3f #pm %1.3f MeV") %  AVERAGE_SIGMAW_CBS % AVERAGE_SIGMAW_CBS_ERROR;
  std::map<std::string, TLatex> latex;
  for(auto & f : format)
  {
    latex[f.first].SetTextSize(0.038);
    latex[f.first].SetTextColor(kBlack);       
    latex[f.first].DrawLatex(xx,y[f.first],f.second.str().c_str());
  }

  string pdf_file = OUTPUT_FILE+".pdf";
  string root_file = OUTPUT_FILE+".root";
  gPad->SaveAs(pdf_file.c_str());
  //gPad->SaveAs(root_file.c_str());
  //diviation_c->Write();
  mcanvas->Write();
  f.Close();
  if(!opt.count("exit")) theApp->Run();
  return 0;
}

static int FCNcall=0;

void fcnResMult(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //calculate chisquare
  //for(int i=0;i<18;i++)
  //{
  //  std::cout << "par" << i <<  " = " << par[i] << std::endl;
  //}
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
    std::cout << std::setw(10) << "W, MeV";
    std::cout << std::setw(10) << "L,1/nb";
    std::cout << std::setw(20) << "Nmh";
    std::cout << std::setw(10) << "sigMH,nb";
    std::cout << std::setw(10) << "Lmh,1/nb";
    std::cout << std::setw(12) << "Nlum";
    std::cout << std::setw(10) << "sigLum,nb";
    std::cout << std::setw(10) << "L,1/nb";
    for(unsigned i=0;i<4;++i) std::cout << setw(7) << "p"+std::string(boost::lexical_cast<std::string>(i));
    std::cout << std::endl;

  }
  for (Int_t i=0;i<EInScan.size();i++)
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
      SigmaWChi2+= sq((parmh[idRSw] - SigmaWInScan[i])/dSigmaWInScan[i]);
    }
    if(USE_CBS_SIGMAW_EACH_POINT)
    {
      double dS = par[PAR_INDEX_DSIGMA + i];
      SigmaWChi2+= sq(dS/dSigmaWInScan[i]);
      CBS_SIGMA_W_IN_CURRENT_POINT = SigmaWInScan[i] + dS;
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
      case BESLUM:
        lumFull=LumInScan[i];
        Nexp=NLum[i];
        nFull = sq(LumInScan[i]/LumInScanError[i]);
        break;
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
      std::cout << boost::format("%10.3f") % (2*Energy);
      std::cout << boost::format("%10.3f") % lumFull;
      std::cout << boost::format("%20d")   %  (unsigned long long)(NmhInScan[i]);
      std::cout << boost::format("%10.3f") % sigmaMH;
      std::cout << boost::format("%10.3f") % (NmhInScan[i]/sigmaMH);
      std::cout << boost::format("%12d")   %  NLum[i];
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
        //chisqmh = sq(N - NmhInScan[i])/(NmhInScan[i] + N*N/nFull); //just chi square
      } 
      else chisqmh = 2*(NmhInScan[i]*log(NmhInScan[i]/N) +N - NmhInScan[i]); //likelihood for low statistics
			SignalDiscrepancy[i] = NmhInScan[i]-N;
			//SignalDiscrepancy[i] = (NmhInScan[i]-N)/sqrt(NmhInScan[i]*(1. + NmhInScan[i]/nFull));
			SignalDiscrepancyError[i] = sqrt(NmhInScan[i]*(1. + NmhInScan[i]/nFull));
			//SignalDiscrepancyError[i] = sqrt(NmhInScan[i]);
			//SignalDiscrepancyError[i] = 1;
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
  if(USE_CBS_SIGMAW_EACH_POINT) 
  {
    f+=SigmaWChi2;
    //std::cout << "f=" << f << "  sigmaw_chi2 = " << SigmaWChi2 <<std::endl;
    /*  
    for(int i=0;i<npar;i++)
    {
      std::cout << " par" << i << " = " << parmh[i] << std::endl;
    }
    for(int i=0;i<7;i++)
    {

      std::cout << " sW" << i << " = " << SigmaWInScan[i] << "  +- " << dSigmaWInScan[i] << std::endl;
    }
    */
  }
  if(MinChi2>f)
  {
    MinChi2=f;
    CHI2_ENERGY = EnergyChi2;
    CHI2_TOTAL=f;
    CHI2_LUM = LumChi2;
    CHI2_SIGNAL = SignalChi2;
    CHI2_CBS_SIGMAW = SigmaWChi2;
  }
  FCNcall=1;
}


TCanvas * draw_signal_and_energy_deviation(const std::vector<double> & parRes)
{

  double EnergyChi2=0;
	TGraphErrors * dEgr = new TGraphErrors;//energy deviation graph
  dEgr->SetName("dE");
	TGraphErrors * dNgr = new TGraphErrors;//signal discrepancy graph
  dNgr->SetName("dN");
	gStyle->SetOptFit(kTRUE);
  for(int is=0;is<EInScan.size();is++)
  {       
    WErrInScan[is]=WErrInScan[is];
    //if(arguments.FreeEnergy==1) WInScan[is]+=2.*+parRes[is+4];
    EnergyChi2+= sq(parRes[is+4]/EErrInScan[is]);
		dEgr->SetPoint(is, WInScan[is], parRes[is+4]*2);
		dEgr->SetPointError(is, 0,  WErrInScan[is]);
		dNgr->SetPoint(is, WInScan[is], SignalDiscrepancy[is]);
		dNgr->SetPointError(is, WErrInScan[is], SignalDiscrepancyError[is]);
  }

  TCanvas * canvas = new TCanvas("dc","Canvas deviation", 777, 800);

  if(FREE_ENERGY_FIT)
  {
    cout << "Draw energy deviation" << endl;
    //TCanvas * dEc = new TCanvas("dEc", "Energy deviation", 777,480);
    canvas->Divide(1,2);
    canvas->cd(1);
    dEgr->SetMarkerStyle(21);
    dEgr->SetLineWidth(2);
    dEgr->SetMarkerSize(1);
    dEgr->Draw("ap");
    dEgr->GetXaxis()->SetTitle("W, MeV");
    dEgr->GetYaxis()->SetTitle("W_{exp}-W_{CBS}, MeV");
    dEgr->Fit("pol0", "Q");
    canvas->cd(2);
  }
  //cout << "Draw signal deviation" << endl;
	//TCanvas * dNc = new TCanvas("dNc", "Number of visible events of signal deviation",777,480);
	dNgr->SetMarkerStyle(21);
	dNgr->SetLineWidth(2);
	dNgr->SetMarkerSize(1);
	dNgr->Draw("ap");
	dNgr->GetXaxis()->SetTitle("W, MeV");
	dNgr->GetYaxis()->SetTitle("N_{vis} - N_{exp}");
	dNgr->Fit("pol0", "Q");
  dEgr->Write();
  dNgr->Write();
  canvas->Write();
  return canvas;
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
  os << sharp << boost::format(" --par=\"%5.1f %4.2f %6.3f %4.3f\"") % P[0] % P[1] %P[2] % P[3] << std::endl;
  //double grad=0;
  //double fval=0;
  //int flag=1;
  //minuit->Eval(N, &grad, fval, &P[0], flag);

  //now calculate and print the average SigmaW from CBS
  ibn::phys_averager sw_aver;
  //TGraphErrors * g = new TGraphErrors;
  for( int i=0;i<SigmaWInScan.size();i++)
  {
    sw_aver.add(SigmaWInScan[i], dSigmaWInScan[i]);
    //std::cout << "sw = " << SigmaWInScan[i] <<  " +- " << dSigmaWInScan[i] << std::endl;
    //g->SetPoint(i, i, SigmaWInScan[i]);
    //g->SetPointError(i,0, dSigmaWInScan[i]);
  }
  //os << sharp  <<  " " << bf( "<Sw> = %5.3f +- %5.3f MeV") %  sw_aver.average() % sw_aver.sigma_average() << std::endl;
  //g->Fit("pol0");

  int nf = minuit->GetNumFreePars();
  if(FREE_ENERGY_FIT) nf-=EInScan.size();
  int ndf = (NpPP-nf);
  double chi2= MinChi2/ndf;
  double prob = TMath::Prob(MinChi2,NpPP-nf);
  typedef boost::format bf;
  os << sharp << " Contribution to chi square:" << endl;
  os << sharp << "   chi2 signal: " << CHI2_SIGNAL <<  " or " << CHI2_SIGNAL/CHI2_TOTAL*100 << "%" << endl;
  os << sharp << "   chi2 energy: " << CHI2_ENERGY <<  " or " << CHI2_ENERGY/CHI2_TOTAL*100 << "%" <<  endl;
  os << sharp << "      chi2 lum: " << CHI2_LUM    <<  " or " << CHI2_LUM/CHI2_TOTAL*100 << "%" << endl;
  os << sharp << "chi2 cbs sigma w: " << CHI2_CBS_SIGMAW  <<  " or " << CHI2_CBS_SIGMAW/CHI2_TOTAL*100 << "%" << endl;
  os << sharp << "    Total chi2: " << CHI2_TOTAL   << endl;
  os << sharp << " ndf    = " << ndf  << endl;
  os << sharp << " chi2/ndf = " <<MinChi2 << "/(" << NpPP<<"-"<<nf<<") = "  << chi2 <<  endl;
  os << sharp << " P(chi2) =" << TMath::Prob(MinChi2,NpPP-nf) << endl;
  os << sharp << " Minuit Mass = "<< bf("%8.3f") % (PDGMASS+P[2]*2.)<<endl;
  os << sharp << " PDG Mass = "   << bf("%8.3f") % PDGMASS <<  endl;
  os << sharp << " " << bf("M-Mpdg = %8.3f +- %5.3f MeV") % (P[2]*2.) %  (dP[2]*2.) << endl;
  os << sharp << " " << bf("       SW = %5.3f +- %5.3f MeV") % P[3] %  dP[3] << endl;
  os << sharp << " " << bf(" <SW_CBS> = %5.3f +- %5.3f MeV") %   sw_aver.average() % sw_aver.sigma_average() << endl;
  os << sharp << " " << bf("      eps = %5.2f +- %4.2f %%") % (P[1]*100) %  (dP[1]*100) << endl;
  os << sharp << " " << bf("       bg = %5.2f +- %5.2f nb") % P[0] %  dP[0] << endl;

  AVERAGE_SIGMAW_CBS = sw_aver.average();
  AVERAGE_SIGMAW_CBS_ERROR = sw_aver.sigma_average();
}

TF1 * get_result_function(const std::vector<double> & parRes, double Emin, double Emax)
{
  TF1 * FitPsiP=0;
  switch(RESONANCE)
  {
    case JPSIRES:
      FitPsiP  = new TF1("FitJPsi",FCrSJpsiAzimov,Emin/2.*ScaleEGr,Emax/2.*ScaleEGr,idRNP);
      FitPsiP->SetTitle(("J/#psi "+CFG_TITLE).c_str());
      //GrRes->SetTitle("J/#psi scan");
      //mcanvas->SetTitle("J/psi");
      break;
    case PSI2SRES:
      FitPsiP  = new TF1("FitPsiP",FCrSPPrimeAzimov,Emin/2.*ScaleEGr,Emax/2.*ScaleEGr,idRNP);
      FitPsiP->SetTitle(("#psi(2S) "+CFG_TITLE).c_str());
      //GrRes->SetTitle("#psi(2S) scan");
      //mcanvas->SetTitle("psi(2S)");
      break;
  }
  FitPsiP->SetLineColor(kRed);

  std::vector<double> parPsiPF(idRNP);
  parPsiPF[idRbg]=parRes[0];
  parPsiPF[idReff]=parRes[1];
  parPsiPF[idRM]=parRes[2];//parPsiP[Iscan][ippeff];   
  parPsiPF[idRSw]=parRes[3];   
  for(int i=4; i<idRNP;i++) parPsiPF[i]=0;
  //parPsiPF[idRFreeGee]=0;
  //parPsiPF[idRTauEff]=0;
  FitPsiP->SetParameters(&parPsiPF[0]);
  return FitPsiP;
}

/*
 * =====================================================================================
 *
 *       Filename:  sim.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  15.04.2018 08:27:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#include "../psifit/FitOniumRCompactLib.hh"
#include <vector>
#include <iostream>
#include <chrono>
using namespace std;

#include <TF1.h>

enum ResonanceType
{
  AUTORES=0,
  JPSIRES=1,
  PSI2SRES=2,
};

unsigned RESONANCE=PSI2SRES;

TF1 * get_result_function(const std::vector<double> & parRes, double Emin, double Emax)
{
  TF1 * FitPsiP=0;
  switch(RESONANCE)
  {
    case JPSIRES:
      FitPsiP  = new TF1("FitJPsi",FCrSJpsiAzimov,Emin/2.*ScaleEGr,Emax/2.*ScaleEGr,idRNP);
      FitPsiP->SetTitle("J/#psi");
      //GrRes->SetTitle("J/#psi scan");
      //mcanvas->SetTitle("J/psi");
      break;
    case PSI2SRES:
      FitPsiP  = new TF1("FitPsiP",FCrSPPrimeAzimov,Emin/2.*ScaleEGr,Emax/2.*ScaleEGr,idRNP);
      FitPsiP->SetTitle("#psi(2S) ");
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

/* 
#define  idRbg       0
#define  idReff      1
#define  idRM        2
#define  idRSw       3
*/

double cm_energy(double E1, double E2)
{
  return 2*sqrt(E1*E2)*cos(0.011);
}
int main(int argc, char ** argv)
{
  const double MJPSI=_MJPsi;
  const double MPSI = _MPsiPrime;
  double Lpoint = 1000; //nb
  double bg = 0;
  double M = 0;

 

  //vector<double> Epoints = { 1544, 1547.8, 1548.2, 1548.6, 1549, 1549.4, 1552};
  //vector<double> Epoints = { 1838, 1841.9, 1842.5, 1843.1, 1843.8, 1844.5, 1847};
  //vector<double> Epoints = { 1838, 1841.9, 1842.5, 1843.1, 1843.8, 1844.5, 1847};
  vector<double> Epoints = { 1838, 1841.9, 1842.2, 1842.5, 1843.1, 1843.8,1844.15, 1844.5, 1847};
  vector<double> dWpoint = { 0.3,     0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,  0.3};
  vector<double> Lpoints;
  for(auto E : Epoints)
  {
    std::cout <<  E*2 << "  " <<  E*2-MPSI << "  "  << E-MPSI/2. << endl;
  }
  ifstream ifs("scenario.txt");
  int n;
  double E, dE, L;
  Epoints.resize(0);
  dWpoint.resize(0);
  while(ifs >> n >> E >> dE >> L )
  {
    Epoints.push_back(E);
    dWpoint.push_back(dE);
    Lpoints.push_back(L);
    std::cout << n << " " << E << "  " << dE << "  " << L << endl;
  }
  for(auto & E: Epoints) { E+=MPSI/2; }

  std::vector<double> par(idRNP);
  par[idRbg] = 7;
  par[idReff] = 0.60;
  par[idRM] = 0;
  par[idRSw] = 1.324;
  for(int i=4;i<idRNP; i++) par[i] = 0;
  //auto f  = get_result_function(par,1540*2,1560*2);
  double MASS = _MPsiPrime;
  auto f  = get_result_function(par,MASS-5*2,MASS+5*2);
  f->SetNpx(1000);
  TRandom R;
  R.SetSeed(chrono::system_clock::now().time_since_epoch().count());
  double sigma_gg = 16; //nb
  double sigma_ee = 16*10; //nb
  int point = 1;
  TGraphErrors * g = new TGraphErrors;
  for (auto & Eset : Epoints)
  {
    double dE=0;
    if( point == 1) dWpoint[point-1]*0.5;
    else  dE = dWpoint[point-2]*0.5;
    //we go to the expected point with the error from previose point
    //E is the real energy
    double  E = R.Gaus(Eset,dE);
    double sigma = f->Eval(E*2);
    double EXTRA_ERROR=0.01;
    Lpoint = Lpoints[point-1];
    unsigned long  Nmh = R.Poisson(sigma*Lpoint)*(1.0+R.Gaus(0,EXTRA_ERROR));
    unsigned long  Ngg = R.Poisson(sigma_gg * Lpoint * pow(E/MASS,-2.))*(1.+R.Gaus(0,EXTRA_ERROR));
    unsigned long  Nee = R.Poisson(sigma_ee * Lpoint * pow(E/MASS,-2.))*(1.+R.Gaus(0,EXTRA_ERROR));
    double dWcbs = dWpoint[point-1];
    double Wcbs = R.Gaus(2*E,dWcbs);
    Wcbs=2*E;
    g->SetPoint(point-1, Wcbs, Nmh/Lpoint);
    g->SetPointError(point-1, dWcbs, sqrt(Nmh)/Lpoint);
    cout << setw(5) << point 
      << setw(10) << Lpoint
      << setw(5) << 0
      << setw(10) << Wcbs
      << setw(5) << dWcbs
      << setw(10) << 1.01
      << setw(5) << 0.1
      << setw(10) <<  Nmh 
      << setw(10) <<  Nee
      << setw(10) << Ngg 
      << endl;
    ++point;
  }
  if(argc > 1)
  {
    TApplication* theApp=new  TApplication("App", &argc, argv);
    g->Draw("a*");
    f->Draw("same");
    theApp->Run();
  }
}

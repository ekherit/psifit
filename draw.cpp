/*
 * =====================================================================================
 *
 *       Filename:  draw.cpp
 *
 *    Description:  Различные вспомогательные функции.
 *
 *        Version:  1.0
 *        Created:  05.12.2011 18:20:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */
#include "draw.h"

#include <iostream>


#include <TROOT.h>
#include <TSystem.h>
#include<TApplication.h>

R__EXTERN TSystem *gSystem;
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("draw","draw", initfuncs);

inline double sq(double x) { return x*x;}
inline double cb(double x) { return x*x*x;}

int main(int argc, char ** argv)
{
  TApplication * theApp=new  TApplication("App", &argc, argv);
  vector <double> EsPsiP;
  vector <double> EsJPsi;
  double sW=1.56; // at PsiPrime
  double dm=0;
  EsPsiP =  {1838.0+dm, 1841.8+dm, 1842.4+dm, 1843.0+dm, 1843.7+dm, 1844.4+dm, 1847.0+dm};
        draw_res(_IdPsiPrime,sW,EsPsiP);
  EsJPsi = {1544.0, 1547.7, 1548.1, 1548.5, 1548.9, 1549.3,1552.0};
  vector<double>  jpsi_add = { 1547.5, 1547.9, 1547.9, 1548.3, 1648.7};
  //EsJPsi.add(jpsi_add);
  vector<double> EsTau = { 1771, 1776.5, 1776.9, 1780.3,1792};
  draw_res(_IdJPsi,sW*sq(_MJPsi/_MPsiPrime),EsJPsi);
  for(double E: EsPsiP) cout << E << " " << E/cos(0.022/2)<< "\t"<<E/cos(0.022/2)-E << endl;
  for(double E: EsJPsi) cout << E << " " << E/cos(0.022/2)<< "\t"<<E/cos(0.022/2)-E << endl;
  for(double E: EsTau) cout << E << " " <<  E/cos(0.022/2)<< "\t"<<E/cos(0.022/2)-E << endl;
  theApp->Run();
}

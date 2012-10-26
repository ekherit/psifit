#ifndef IBN_PSIFIT_UTILS_H
#define IBN_PSIFIT_UTILS_H
// This file is my utils
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <ibn/averager.h>
#include <ibn/math.h>

using ibn::sq;
using ibn::cb;

/* Estimate cross section using BES online lum */
inline double estimate_cross_section 
  (
    double E0 /* normalization energy */,
    size_t size /* array size */, 
    double * E /*Energy */, 
    double * L /* online lum */, 
    double * N /* nomber of events*/
  )
{
  ibn::averager<double> XS; //cross section averager
  for(int i=0; i<size;i++)
  {
    if(L[i]==0) continue;
    XS.add(N[i]*ibn::sq(E[i]/E0)/L[i]);
  }
  return XS.average();
}

/* Keep total luminosity constant */
inline double estimate_cross_section2 
  (
    double E0 /* normalization energy */,
    size_t size /* array size */, 
    double * E /*Energy */, 
    double * L /* online lum */, 
    double * N /* nomber of events*/
  )
{
  double IL=0; //integrated luminosity
  double IN=0; //integrated energy normalized number of events 
  for(int i=0;i<size; i++)
  {
    IL+=L[i];
    IN+=double(N[i])*sq(E[i]/E0);
  }
  if(IL==0) return 0;
  return IN/IL;
}

/* Keep total luminosity constant */
inline double estimate_cross_section2 
  (
    double E0 /* normalization energy */,
    const std::vector<double> &E /*Energy */, 
    const std::vector<double> &L /* online lum */, 
    const std::vector<double> &N /* nomber of events*/
  )
{
  double IL=0; //integrated luminosity
  double IN=0; //integrated energy normalized number of events 
  for(int i=0;i<E.size(); i++)
  {
    IL+=L[i];
    IN+=double(N[i])*sq(E[i]/E0);
  }
  if(IL==0) return 0;
  return IN/IL;
}

//find average of the array
inline double average(size_t size, double *X)
{
  ibn::averager<double> A;
  for(int i=0; i<size;i++)
  {
    A.add(X[i]);
  }
  return A.average();
}

class format
{
  struct Item
  {
    boost::format head_format;
    boost::format data_format;
    Item(string name, unsigned size, string fstr)
    {
      //create format size%s
      head_format = boost::format(string("%")+boost::lexical_cast<string>(size)+"s");

    }
  };
  public:
    format(void) 
    {
    }
};

#endif

//  file: gauss.cpp
//
//  This program determines points and weights for Gaussian quadrature                 
//                                                                     
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      04-Jan-2004  original version, for 780.20 Computational Physics
//
//  Notes:  
//   * compile with:  "g++ -Wall -c gauss.cpp"
//   * adapted from: "Projects in Computational Physics" by Landau and Paez  
//             copyrighted by John Wiley and Sons, New York               
//             code copyrighted by RH Landau  
//   * Needs more careful treatment of integers <--> doubles                         
// 
//************************************************************************

// include files
#include <cmath>
#include "Arsenal.h"

//************************************************************************
void
gauss (int npts, int job, double a, double b, double xpts[], double weights[])
{
  //     npts     number of points                                       
  //     job = 0  rescaling uniformly between (a,b)                      
  //           1  for integral (0,b) with 50% points inside (0, ab/(a+b))
  //           2  for integral (a,inf) with 50% inside (a,b+2a)          
  //     xpts, weights     output grid points and weights.                         

  const double pi = M_PI;
  const double eps = 3.e-10;	// limit for accuracy 

  double t = 0., t1 = 0., p1 = 0., p2 = 0., p3 = 0., pp = 0.;

  int m = (npts + 1) / 2;
  for (int i = 1; i <= m; i++)
  {
    t = cos (pi * (i - 0.25) / (npts + 0.5));
    t1 = 1;
    while ((fabs (t - t1)) >= eps)
    {
      p1 = 1.0;
      p2 = 0.0;
      for (int j = 1; j <= npts; j++)
      {
	p3 = p2;
	p2 = p1;
	p1 = ((2 * j - 1) * t * p2 - (j - 1) * p3) / j;
      }
      pp = npts * (t * p1 - p2) / (t * t - 1);
      t1 = t;
      t = t1 - p1 / pp;
    }
    xpts[i - 1] = -t;
    xpts[npts - i] = t;
    weights[i - 1] = 2.0 / ((1 - t * t) * pp * pp);
    weights[npts - i] = weights[i - 1];
  }

  if (job == 0)		// rescaling uniformly between (a,b) 
  {
    for (int i = 0; i < npts; i++)
    {
      xpts[i] = xpts[i] * (b - a) / 2.0 + (b + a) / 2.0;
      weights[i] = weights[i] * (b - a) / 2.0;
    }
  }
  
  if (job == 1)		// integral (0,b) with 50% points inside (0, ab/(a+b))
  {
    for (int i = 0; i < npts; i++)
    {
      t = (b + a) - (b - a) * xpts[i];
      xpts[i] = a * b * (1 + xpts[i]) / t;
      weights[i] = weights[i] * 2.0 * a * b * b / (t * t);
    }
  }
  
  if (job == 2)		// integral (a,inf) with 50% inside (a,b+2a)
  {
    for (int i = 0; i < npts; i++)
    {
      t = 1.0 - xpts[i];
      xpts[i] = (b * xpts[i] + b + a + a) / t;
      weights[i] = weights[i] * 2.0 * (a + b) / (t * t);
    }
  }
}

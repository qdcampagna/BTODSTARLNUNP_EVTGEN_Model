
/***********************************************************************
* Copyright 1998-2020 CERN for the benefit of the EvtGen authors       *
*                                                                      *
* This file is part of EvtGen.                                         *
*                                                                      *
* EvtGen is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* EvtGen is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     *
***********************************************************************/

#include "EvtGenBase/EvtSpinDensity.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cstring>
using std::endl;
using std::ostream;

EvtSpinDensity::EvtSpinDensity( const EvtSpinDensity& density )
{
    dim = 0;
    rho = 0;
    setDim( density.dim );
    //    memmove(rho,density.rho,dim*dim*sizeof(EvtComplex));
    for(int i=0;i<dim*dim;i++)rho[i] = density.rho[i];
}

EvtSpinDensity& EvtSpinDensity::operator=( const EvtSpinDensity& density )
{
    setDim( density.dim );
    //    memmove(rho,density.rho,dim*dim*sizeof(EvtComplex));
    for(int i=0;i<dim*dim;i++)rho[i] = density.rho[i];
    return *this;
}

EvtSpinDensity::~EvtSpinDensity()
{
  if ( dim != 0 ) delete[] rho;
}

// EvtSpinDensity::EvtSpinDensity()
// {
//     dim = 0;
//     rho = 0;
// }

void EvtSpinDensity::setDim( int n )
{
  if ( dim == n ) return;
  if ( dim != 0 ) {
    delete [] rho;
    rho = 0;
    dim = 0;
  }
  if ( n == 0 ) return;
  dim = n;
  rho = new EvtComplex[n*n];
}

void EvtSpinDensity::setDiag( int n )
{
  setDim( n );
  for (int i = 0; i < n*n; i++ ) rho[i] = 0;
  for (int i = 0; i < n; i++ ) rho[i*dim + i] = 1.0;
}

double EvtSpinDensity::normalizedProb( const EvtSpinDensity& d ) const
{
    if ( dim != d.dim ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Not matching dimensions in NormalizedProb" << endl;
        ::abort();
    }

    double norm = 0.0;
    for (int i = 0; i < dim; i++ ) norm += real( rho[i*dim+i] );

    EvtComplex prob;
    for (int i = 0, imax = dim*dim; i < imax; i++ ) prob += rho[i] * d.rho[i];

    if ( imag( prob ) > 0.00000001 * real( prob ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Imaginary probability:" << prob << " " << norm << endl;
    }
    if ( real( prob ) < 0.0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Negative probability:" << prob << " " << norm << endl;
    }

    return real( prob ) / norm;
}

int EvtSpinDensity::check()
{
    if ( dim < 1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "dim=" << dim << "in SpinDensity::Check" << endl;
    }

    int i, j;

    double trace( 0.0 );

    for ( i = 0; i < dim; i++ ) {
        trace += abs( rho[i*dim+i] );
    }

    for ( i = 0; i < dim; i++ ) {
        if ( real( rho[i*dim+i] ) < 0.0 )
            return 0;
        if ( imag( rho[i*dim+i] ) * 1000000.0 > trace ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << *this << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << trace << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Failing 1" << endl;
            return 0;
        }
    }

    for ( i = 0; i < dim; i++ ) {
        for ( j = i + 1; j < dim; j++ ) {
            if ( fabs( real( rho[i*dim+j] - rho[j*dim+i] ) ) >
                 0.00000001 * ( abs( rho[i*dim+i] ) + abs( rho[j*dim+j] ) ) ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Failing 2" << endl;
                return 0;
            }
            if ( fabs( imag( rho[i*dim+j] + rho[j*dim+i] ) ) >
                 0.00000001 * ( abs( rho[i*dim+i] ) + abs( rho[j*dim+j] ) ) ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Failing 3" << endl;
                return 0;
            }
        }
    }

    return 1;
}

ostream& operator<<( ostream& s, const EvtSpinDensity& d )
{
    s << endl;
    s << "Dimension:" << d.dim << endl;

    for (int i = 0; i < d.dim; i++ ) {
        for (int j = 0; j < d.dim; j++ ) {
	  s << d.get(i,j) << " ";
        }
        s << endl;
    }

    return s;
}

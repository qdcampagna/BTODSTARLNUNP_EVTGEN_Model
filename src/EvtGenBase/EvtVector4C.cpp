
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

#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using std::ostream;


EvtVector4C rotateEuler( const EvtVector4C& rs, double alpha, double beta,
                         double gamma )
{
    EvtVector4C tmp( rs );
    tmp.applyRotateEuler( alpha, beta, gamma );
    return tmp;
}

EvtVector4C boostTo( const EvtVector4C& rs, const EvtVector4R p4 )
{
    EvtVector4C tmp( rs );
    tmp.applyBoostTo( p4 );
    return tmp;
}

EvtVector4C boostTo( const EvtVector4C& rs, const EvtVector3R boost )
{
    EvtVector4C tmp( rs );
    tmp.applyBoostTo( boost );
    return tmp;
}

void EvtVector4C::applyBoostTo( const EvtVector4R& p4 )
{
    double e = p4.get( 0 );

    EvtVector3R boost( p4.get( 1 ) / e, p4.get( 2 ) / e, p4.get( 3 ) / e );

    applyBoostTo( boost );

    return;
}

void EvtVector4C::applyBoostTo( const EvtVector3R& boost )
{
    double bx = boost.get( 0 );
    double by = boost.get( 1 );
    double bz = boost.get( 2 );
    double b2 = bx * bx + by * by + bz * bz;

    assert( b2 < 1.0 );

    double gamma2 = 1/(1 - b2);
    double gamma = sqrt(gamma2);
    double F = gamma2/(gamma + 1);

    EvtComplex bp = bx * v[1] + by * v[2] + bz * v[3];
    EvtComplex c = F * bp + gamma * v[0];

    v[1] += c * bx;
    v[2] += c * by;
    v[3] += c * bz;
    v[0] = gamma * (v[0] + bp);
}

void EvtVector4C::applyRotateEuler( double phi, double theta, double ksi )
{
    double sp = sin( phi );
    double st = sin( theta );
    double sk = sin( ksi );
    double cp = cos( phi );
    double ct = cos( theta );
    double ck = cos( ksi );

    EvtComplex x = ( ck * ct * cp - sk * sp ) * v[1] +
                   ( -sk * ct * cp - ck * sp ) * v[2] + st * cp * v[3];
    EvtComplex y = ( ck * ct * sp + sk * cp ) * v[1] +
                   ( -sk * ct * sp + ck * cp ) * v[2] + st * sp * v[3];
    EvtComplex z = -ck * st * v[1] + sk * st * v[2] + ct * v[3];

    v[1] = x;
    v[2] = y;
    v[3] = z;
}

ostream& operator<<( ostream& s, const EvtVector4C& v )
{
    s << "(" << v.v[0] << "," << v.v[1] << "," << v.v[2] << "," << v.v[3] << ")";

    return s;
}

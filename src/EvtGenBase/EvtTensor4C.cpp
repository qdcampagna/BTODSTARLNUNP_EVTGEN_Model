
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

#include "EvtGenBase/EvtTensor4C.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using std::endl;
using std::ostream;

EvtTensor4C::EvtTensor4C( const EvtTensor4C& t1 )
{
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] = t1.t[i][j];
        }
    }
}

const EvtTensor4C& EvtTensor4C::g()
{
    static EvtTensor4C g_metric( 1.0, -1.0, -1.0, -1.0 );

    return g_metric;
}

EvtTensor4C& EvtTensor4C::operator=( const EvtTensor4C& t1 )
{
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] = t1.t[i][j];
        }
    }
    return *this;
}

EvtTensor4C EvtTensor4C::conj() const
{
    EvtTensor4C temp;

    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            temp.set( j, i, ::conj( t[i][j] ) );
        }
    }
    return temp;
}

EvtTensor4C rotateEuler( const EvtTensor4C& rs, double alpha, double beta,
                         double gamma )
{
    EvtTensor4C tmp( rs );
    tmp.applyRotateEuler( alpha, beta, gamma );
    return tmp;
}

EvtTensor4C boostTo( const EvtTensor4C& rs, const EvtVector4R p4 )
{
    EvtTensor4C tmp( rs );
    tmp.applyBoostTo( p4 );
    return tmp;
}

EvtTensor4C boostTo( const EvtTensor4C& rs, const EvtVector3R boost )
{
    EvtTensor4C tmp( rs );
    tmp.applyBoostTo( boost );
    return tmp;
}

void EvtTensor4C::applyBoostTo( const EvtVector4R& p4 )
{
    double e = p4.get( 0 );

    EvtVector3R boost( p4.get( 1 ) / e, p4.get( 2 ) / e, p4.get( 3 ) / e );

    applyBoostTo( boost );

    return;
}

void EvtTensor4C::applyBoostTo( const EvtVector3R& boost )
{
    double bx, by, bz, gamma, b2;
    double lambda[4][4];
    EvtComplex tt[4][4];

    bx = boost.get( 0 );
    by = boost.get( 1 );
    bz = boost.get( 2 );

    double bxx = bx * bx;
    double byy = by * by;
    double bzz = bz * bz;

    b2 = bxx + byy + bzz;

    if ( b2 == 0.0 ) {
        return;
    }

    assert( b2 < 1.0 );

    gamma = 1.0 / sqrt( 1 - b2 );

    int i, j, k;

    if ( b2 == 0.0 ) {
        return;
    }

    lambda[0][0] = gamma;
    lambda[0][1] = gamma * bx;
    lambda[1][0] = gamma * bx;
    lambda[0][2] = gamma * by;
    lambda[2][0] = gamma * by;
    lambda[0][3] = gamma * bz;
    lambda[3][0] = gamma * bz;

    lambda[1][1] = 1.0 + ( gamma - 1.0 ) * bx * bx / b2;
    lambda[2][2] = 1.0 + ( gamma - 1.0 ) * by * by / b2;
    lambda[3][3] = 1.0 + ( gamma - 1.0 ) * bz * bz / b2;

    lambda[1][2] = ( gamma - 1.0 ) * bx * by / b2;
    lambda[2][1] = ( gamma - 1.0 ) * bx * by / b2;

    lambda[1][3] = ( gamma - 1.0 ) * bx * bz / b2;
    lambda[3][1] = ( gamma - 1.0 ) * bx * bz / b2;

    lambda[3][2] = ( gamma - 1.0 ) * bz * by / b2;
    lambda[2][3] = ( gamma - 1.0 ) * bz * by / b2;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            tt[i][j] = EvtComplex( 0.0 );
            for ( k = 0; k < 4; k++ ) {
                tt[i][j] = tt[i][j] + lambda[j][k] * t[i][k];
            }
        }
    }

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] = EvtComplex( 0.0 );
            for ( k = 0; k < 4; k++ ) {
                t[i][j] = t[i][j] + lambda[i][k] * tt[k][j];
            }
        }
    }
}

void EvtTensor4C::zero()
{
    int i, j;
    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] = EvtComplex( 0.0, 0.0 );
        }
    }
}

ostream& operator<<( ostream& s, const EvtTensor4C& t )
{
    int i, j;
    s << endl;
    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            s << t.t[i][j];
        }
        s << endl;
    }
    return s;
}

void EvtTensor4C::setdiag( double g00, double g11, double g22, double g33 )
{
    t[0][0] = EvtComplex( g00 );
    t[1][1] = EvtComplex( g11 );
    t[2][2] = EvtComplex( g22 );
    t[3][3] = EvtComplex( g33 );
    t[0][1] = EvtComplex( 0.0 );
    t[0][2] = EvtComplex( 0.0 );
    t[0][3] = EvtComplex( 0.0 );
    t[1][0] = EvtComplex( 0.0 );
    t[1][2] = EvtComplex( 0.0 );
    t[1][3] = EvtComplex( 0.0 );
    t[2][0] = EvtComplex( 0.0 );
    t[2][1] = EvtComplex( 0.0 );
    t[2][3] = EvtComplex( 0.0 );
    t[3][0] = EvtComplex( 0.0 );
    t[3][1] = EvtComplex( 0.0 );
    t[3][2] = EvtComplex( 0.0 );
}

EvtTensor4C& EvtTensor4C::operator+=( const EvtTensor4C& t2 )
{
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] += t2.get( i, j );
        }
    }
    return *this;
}

EvtTensor4C& EvtTensor4C::operator-=( const EvtTensor4C& t2 )
{
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] -= t2.get( i, j );
        }
    }
    return *this;
}

EvtTensor4C& EvtTensor4C::operator*=( const EvtComplex& c )
{
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] *= c;
        }
    }
    return *this;
}

EvtTensor4C operator*( const EvtTensor4C& t1, const EvtComplex& c )
{
    return EvtTensor4C( t1 ) *= c;
}

EvtTensor4C operator*( const EvtComplex& c, const EvtTensor4C& t1 )
{
    return EvtTensor4C( t1 ) *= c;
}

EvtTensor4C& EvtTensor4C::operator*=( double d )
{
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] *= d;
        }
    }
    return *this;
}

EvtTensor4C operator*( const EvtTensor4C& t1, double d )
{
    return EvtTensor4C( t1 ) *= d;
}

EvtTensor4C operator*( double d, const EvtTensor4C& t1 )
{
    return EvtTensor4C( t1 ) *= d;
}

EvtComplex cont( const EvtTensor4C& t1, const EvtTensor4C& t2 )
{
    EvtComplex sum( 0.0, 0.0 );
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            if ( ( i == 0 && j != 0 ) || ( j == 0 && i != 0 ) ) {
                sum -= t1.t[i][j] * t2.t[i][j];
            } else {
                sum += t1.t[i][j] * t2.t[i][j];
            }
        }
    }

    return sum;
}

EvtTensor4C EvtGenFunctions::directProd( const EvtVector4C& c1,
                                         const EvtVector4C& c2 )
{
    EvtTensor4C temp;
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            temp.set( i, j, c1.get( i ) * c2.get( j ) );
        }
    }
    return temp;
}

EvtTensor4C EvtGenFunctions::directProd( const EvtVector4C& c1,
                                         const EvtVector4R& c2 )
{
    EvtTensor4C temp;
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            temp.set( i, j, c1.get( i ) * c2.get( j ) );
        }
    }
    return temp;
}

EvtTensor4C EvtGenFunctions::directProd( const EvtVector4R& c1,
                                         const EvtVector4R& c2 )
{
    EvtTensor4C temp;
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            temp.t[i][j] = EvtComplex( c1.get( i ) * c2.get( j ), 0.0 );
        }
    }
    return temp;
}

EvtTensor4C& EvtTensor4C::addDirProd( const EvtVector4R& p1,
                                      const EvtVector4R& p2 )
{
    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] += p1.get( i ) * p2.get( j );
        }
    }
    return *this;
}

EvtTensor4C dual( const EvtTensor4C& t2 )
{
    EvtTensor4C temp;

    temp.set( 0, 0, EvtComplex( 0.0, 0.0 ) );
    temp.set( 1, 1, EvtComplex( 0.0, 0.0 ) );
    temp.set( 2, 2, EvtComplex( 0.0, 0.0 ) );
    temp.set( 3, 3, EvtComplex( 0.0, 0.0 ) );

    temp.set( 0, 1, t2.get( 3, 2 ) - t2.get( 2, 3 ) );
    temp.set( 0, 2, -t2.get( 3, 1 ) + t2.get( 1, 3 ) );
    temp.set( 0, 3, t2.get( 2, 1 ) - t2.get( 1, 2 ) );

    temp.set( 1, 2, -t2.get( 3, 0 ) + t2.get( 0, 3 ) );
    temp.set( 1, 3, t2.get( 2, 0 ) - t2.get( 0, 2 ) );

    temp.set( 2, 3, -t2.get( 1, 0 ) + t2.get( 0, 1 ) );

    temp.set( 1, 0, -temp.get( 0, 1 ) );
    temp.set( 2, 0, -temp.get( 0, 2 ) );
    temp.set( 3, 0, -temp.get( 0, 3 ) );

    temp.set( 2, 1, -temp.get( 1, 2 ) );
    temp.set( 3, 1, -temp.get( 1, 3 ) );

    temp.set( 3, 2, -temp.get( 2, 3 ) );

    return temp;
}

EvtTensor4C conj( const EvtTensor4C& t2 )
{
    EvtTensor4C temp;

    int i, j;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            temp.set( i, j, ::conj( ( t2.get( i, j ) ) ) );
        }
    }

    return temp;
}

EvtTensor4C cont22( const EvtTensor4C& t1, const EvtTensor4C& t2 )
{
    EvtTensor4C temp;

    int i, j;
    EvtComplex c;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            c = t1.get( i, 0 ) * t2.get( j, 0 ) -
                t1.get( i, 1 ) * t2.get( j, 1 ) -
                t1.get( i, 2 ) * t2.get( j, 2 ) - t1.get( i, 3 ) * t2.get( j, 3 );
            temp.set( i, j, c );
        }
    }

    return temp;
}

EvtTensor4C cont11( const EvtTensor4C& t1, const EvtTensor4C& t2 )
{
    EvtTensor4C temp;

    int i, j;
    EvtComplex c;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            c = t1.get( 0, i ) * t2.get( 0, j ) -
                t1.get( 1, i ) * t2.get( 1, j ) -
                t1.get( 2, i ) * t2.get( 2, j ) - t1.get( 3, i ) * t2.get( 3, j );
            temp.set( i, j, c );
        }
    }

    return temp;
}

EvtVector4C EvtTensor4C::cont1( const EvtVector4C& v4 ) const
{
    EvtVector4C temp;

    int i;

    for ( i = 0; i < 4; i++ ) {
        temp.set( i, t[0][i] * v4.get( 0 ) - t[1][i] * v4.get( 1 ) -
                         t[2][i] * v4.get( 2 ) - t[3][i] * v4.get( 3 ) );
    }

    return temp;
}

EvtVector4C EvtTensor4C::cont2( const EvtVector4C& v4 ) const
{
    EvtVector4C temp;

    int i;

    for ( i = 0; i < 4; i++ ) {
        temp.set( i, t[i][0] * v4.get( 0 ) - t[i][1] * v4.get( 1 ) -
                         t[i][2] * v4.get( 2 ) - t[i][3] * v4.get( 3 ) );
    }

    return temp;
}

EvtVector4C EvtTensor4C::cont1( const EvtVector4R& v4 ) const
{
    EvtVector4C temp;

    int i;

    for ( i = 0; i < 4; i++ ) {
        temp.set( i, t[0][i] * v4.get( 0 ) - t[1][i] * v4.get( 1 ) -
                         t[2][i] * v4.get( 2 ) - t[3][i] * v4.get( 3 ) );
    }

    return temp;
}

EvtVector4C EvtTensor4C::cont2( const EvtVector4R& v4 ) const
{
    EvtVector4C temp;

    int i;

    for ( i = 0; i < 4; i++ ) {
        temp.set( i, t[i][0] * v4.get( 0 ) - t[i][1] * v4.get( 1 ) -
                         t[i][2] * v4.get( 2 ) - t[i][3] * v4.get( 3 ) );
    }

    return temp;
}

void EvtTensor4C::applyRotateEuler( double phi, double theta, double ksi )
{
    EvtComplex tt[4][4];
    double sp, st, sk, cp, ct, ck;
    double lambda[4][4];

    sp = sin( phi );
    st = sin( theta );
    sk = sin( ksi );
    cp = cos( phi );
    ct = cos( theta );
    ck = cos( ksi );

    lambda[0][0] = 1.0;
    lambda[0][1] = 0.0;
    lambda[1][0] = 0.0;
    lambda[0][2] = 0.0;
    lambda[2][0] = 0.0;
    lambda[0][3] = 0.0;
    lambda[3][0] = 0.0;

    lambda[1][1] = ck * ct * cp - sk * sp;
    lambda[1][2] = -sk * ct * cp - ck * sp;
    lambda[1][3] = st * cp;

    lambda[2][1] = ck * ct * sp + sk * cp;
    lambda[2][2] = -sk * ct * sp + ck * cp;
    lambda[2][3] = st * sp;

    lambda[3][1] = -ck * st;
    lambda[3][2] = sk * st;
    lambda[3][3] = ct;

    int i, j, k;

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            tt[i][j] = EvtComplex( 0.0 );
            for ( k = 0; k < 4; k++ ) {
                tt[i][j] += lambda[j][k] * t[i][k];
            }
        }
    }

    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            t[i][j] = EvtComplex( 0.0 );
            for ( k = 0; k < 4; k++ ) {
                t[i][j] += lambda[i][k] * tt[k][j];
            }
        }
    }
}

EvtTensor4C& EvtTensor4C::addScaled(EvtComplex s, const EvtTensor4C& a){
  EvtComplex *u0 = reinterpret_cast<EvtComplex*>(t);
  const EvtComplex *u1 = reinterpret_cast<const EvtComplex*>(a.t);
  for(int i=0;i<16;i++) u0[i] += s*u1[i];
  return *this;
}

EvtTensor4C EvtGenFunctions::asymProd( const EvtVector4R& a, const EvtVector4R& b){
  // t_{ij} = eps_{ijkl} a^{k} b^{l}
  // eps_{ijkl} is the Levi-Civita symbol -- antisymmetric tensor

  EvtTensor4C res;
  EvtComplex (&t)[4][4] = res.t;
  double t01 = a.get(2)*b.get(3) - a.get(3)*b.get(2); t[0][1] = t01; t[1][0] = -t01;
  double t02 = a.get(3)*b.get(1) - a.get(1)*b.get(3); t[0][2] = t02; t[2][0] = -t02;
  double t03 = a.get(1)*b.get(2) - a.get(2)*b.get(1); t[0][3] = t03; t[3][0] = -t03;
  double t12 = a.get(0)*b.get(3) - a.get(3)*b.get(0); t[1][2] = t12; t[2][1] = -t12;
  double t13 = a.get(2)*b.get(0) - a.get(0)*b.get(2); t[1][3] = t13; t[3][1] = -t13;
  double t23 = a.get(0)*b.get(1) - a.get(1)*b.get(0); t[2][3] = t23; t[3][2] = -t23;
  return res;
}

EvtVector4R EvtGenFunctions::asymProd(const EvtVector4R& _a, const EvtVector4R& _b, const EvtVector4R& _c){
  // t_{ij} = eps_{ijkl} a^{j} b^{k} c^{l}
  // eps_{ijkl} is the Levi-Civita symbol -- antisymmetric tensor

  double a[4] = {_a.get(0), _a.get(1), _a.get(2), _a.get(3)};
  /*
  double b[4] = {_b.get(0), _b.get(1), _b.get(2), _b.get(3)};
  double c[4] = {_c.get(0), _c.get(1), _c.get(2), _c.get(3)};
  double t[4] = {0};

  t[0] += a[1]*b[2]*c[3];
  t[0] -= a[1]*b[3]*c[2];
  t[0] -= a[2]*b[1]*c[3];
  t[0] += a[2]*b[3]*c[1];
  t[0] += a[3]*b[1]*c[2];
  t[0] -= a[3]*b[2]*c[1];

  t[1] -= a[0]*b[2]*c[3];
  t[1] += a[0]*b[3]*c[2];
  t[1] += a[2]*b[0]*c[3];
  t[1] -= a[2]*b[3]*c[0];
  t[1] -= a[3]*b[0]*c[2];
  t[1] += a[3]*b[2]*c[0];

  t[2] += a[0]*b[1]*c[3];
  t[2] -= a[0]*b[3]*c[1];
  t[2] -= a[1]*b[0]*c[3];
  t[2] += a[1]*b[3]*c[0];
  t[2] += a[3]*b[0]*c[1];
  t[2] -= a[3]*b[1]*c[0];

  t[3] -= a[0]*b[1]*c[2];
  t[3] += a[0]*b[2]*c[1];
  t[3] += a[1]*b[0]*c[2];
  t[3] -= a[1]*b[2]*c[0];
  t[3] -= a[2]*b[0]*c[1];
  t[3] += a[2]*b[1]*c[0];

  EvtVector4R res(t[0],t[1],t[2],t[3]);
  //  std::cout<<res<<std::endl;
  return res;
  */
  double t01 = _b.get(2)*_c.get(3) - _b.get(3)*_c.get(2);
  double t02 = _b.get(3)*_c.get(1) - _b.get(1)*_c.get(3);
  double t03 = _b.get(1)*_c.get(2) - _b.get(2)*_c.get(1);
  double t12 = _b.get(0)*_c.get(3) - _b.get(3)*_c.get(0);
  double t13 = _b.get(2)*_c.get(0) - _b.get(0)*_c.get(2);
  double t23 = _b.get(0)*_c.get(1) - _b.get(1)*_c.get(0);
  double u[4] = {0};
  u[0] =   t01*a[1] +  t02*a[2] + t03*a[3];
  u[1] =   t12*a[2] +  t13*a[3] - t01*a[0];
  u[2] =   t23*a[3] - (t02*a[0] + t12*a[1]);
  u[3] = -(t03*a[0] +  t13*a[1] + t23*a[2]);
  //    std::cout<<u[0]<<" "<<u[1]<<" "<<u[2]<<" "<<u[3]<<"\n";
  return EvtVector4R(u[0], u[1], u[2], u[3]);
}

EvtVector4C EvtGenFunctions::asymProd(const EvtVector4C& _a, const EvtVector4R& _b, const EvtVector4R& _c){
  // t_{ij} = eps_{ijkl} a^{j} b^{k} c^{l}
  // eps_{ijkl} is the Levi-Civita symbol -- antisymmetric tensor

  EvtComplex a[4] = {_a.get(0), _a.get(1), _a.get(2), _a.get(3)};
  double t01 = _b.get(2)*_c.get(3) - _b.get(3)*_c.get(2);
  double t02 = _b.get(3)*_c.get(1) - _b.get(1)*_c.get(3);
  double t03 = _b.get(1)*_c.get(2) - _b.get(2)*_c.get(1);
  double t12 = _b.get(0)*_c.get(3) - _b.get(3)*_c.get(0);
  double t13 = _b.get(2)*_c.get(0) - _b.get(0)*_c.get(2);
  double t23 = _b.get(0)*_c.get(1) - _b.get(1)*_c.get(0);
  return EvtVector4C(t01*a[1] +  t02*a[2] + t03*a[3],
		     t12*a[2] +  t13*a[3] - t01*a[0],
		     t23*a[3] - (t02*a[0] + t12*a[1]),
		     -(t03*a[0] +  t13*a[1] + t23*a[2]));
}

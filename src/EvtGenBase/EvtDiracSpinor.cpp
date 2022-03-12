
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

#include "EvtGenBase/EvtDiracSpinor.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <assert.h>
#include <math.h>
using std::ostream;

ostream& operator<<( ostream& s, const EvtDiracSpinor& sp )
{
    s << "[" << sp.spinor[0] << "," << sp.spinor[1] << "," << sp.spinor[2]
      << "," << sp.spinor[3] << "]";
    return s;
}

EvtDiracSpinor rotateEuler( const EvtDiracSpinor& sp, double alpha, double beta,
                            double gamma )
{
    EvtDiracSpinor tmp( sp );
    tmp.applyRotateEuler( alpha, beta, gamma );
    return tmp;
}

EvtDiracSpinor boostTo( const EvtDiracSpinor& sp, const EvtVector4R p4 )
{
    EvtDiracSpinor tmp( sp );
    tmp.applyBoostTo( p4 );
    return tmp;
}

EvtDiracSpinor boostTo( const EvtDiracSpinor& sp, const EvtVector3R boost )
{
    EvtDiracSpinor tmp( sp );
    tmp.applyBoostTo( boost );
    return tmp;
}

void EvtDiracSpinor::applyBoostTo( const EvtVector4R& p4 )
{
#if 0
    double e = p4.get( 0 );

    EvtVector3R boost( p4.get( 1 ) / e, p4.get( 2 ) / e, p4.get( 3 ) / e );

    applyBoostTo( boost );

    return;
#else
    double px = p4.get(1), py = p4.get(2), pz = p4.get(3), e = p4.get(0);
    double p2 = px*px + py*py + pz*pz, m2 = e*e - p2, m = sqrt(m2), em = e + m;
    double f2 = sqrt(0.5/(em*m)), f1 = f2*em;
    double bx = f2*px, by = f2*py, bz = f2*pz;
    EvtComplex t = spinor[0], x = spinor[1], y = spinor[2], z = spinor[3], bxy(bx, by);

    spinor[0] = f1 * t + (bxy & z) + bz * y;
    spinor[1] = f1 * x + (bxy * y) - bz * z;
    spinor[2] = f1 * y + (bxy & x) + bz * t;
    spinor[3] = f1 * z + (bxy * t) - bz * x;
#endif
}

void EvtDiracSpinor::applyBoostTo( const EvtVector3R& boost )
{
#if 0
    double bx, by, bz, gamma, b2, f1, f2;
    EvtComplex spinorp[4];

    bx = boost.get( 0 );
    by = boost.get( 1 );
    bz = boost.get( 2 );
    b2 = bx * bx + by * by + bz * bz;

    if ( b2 == 0.0 ) {
        return;
    }

    //assert(b2<1.0);

    gamma = 1.0;
    if ( b2 < 1.0 ) {
        gamma = 1.0 / sqrt( 1.0 - b2 );
    }

    f1 = sqrt( ( gamma + 1.0 ) / 2.0 );
    f2 = f1 * gamma / ( gamma + 1.0 );

    spinorp[0] = f1 * spinor[0] + f2 * bz * spinor[2] +
                 f2 * EvtComplex( bx, -by ) * spinor[3];
    spinorp[1] = f1 * spinor[1] + f2 * EvtComplex( bx, by ) * spinor[2] -
                 f2 * bz * spinor[3];
    spinorp[2] = f2 * bz * spinor[0] + f2 * EvtComplex( bx, -by ) * spinor[1] +
                 f1 * spinor[2];
    spinorp[3] = f2 * EvtComplex( bx, by ) * spinor[0] - f2 * bz * spinor[1] +
                 f1 * spinor[3];

    spinor[0] = spinorp[0];
    spinor[1] = spinorp[1];
    spinor[2] = spinorp[2];
    spinor[3] = spinorp[3];

    return;
#else
    double bx = boost.get( 0 ), by = boost.get( 1 ), bz = boost.get( 2 );
    double b2 = bx * bx + by * by + bz * bz, c2 = 1 - b2, c = sqrt( c2 );
    double f2 = sqrt( 0.5 / (c2 + c) );
    bx *= f2;
    by *= f2;
    bz *= f2;
    double f1 = f2 + f2*c;
    EvtComplex t = spinor[0], x = spinor[1], y = spinor[2], z = spinor[3], bxy(bx, by);

    spinor[0] = f1 * t + (bxy & z) + bz * y;
    spinor[1] = f1 * x + (bxy * y) - bz * z;
    spinor[2] = f1 * y + (bxy & x) + bz * t;
    spinor[3] = f1 * z + (bxy * t) - bz * x;
#endif
}

void EvtDiracSpinor::applyRotateEuler( double alpha, double beta, double gamma )
{
    EvtComplex retVal[4];

    double cb2 = cos( 0.5 * beta );
    double sb2 = sin( 0.5 * beta );
    double capg2 = cos( 0.5 * ( alpha + gamma ) );
    double camg2 = cos( 0.5 * ( alpha - gamma ) );
    double sapg2 = sin( 0.5 * ( alpha + gamma ) );
    double samg2 = sin( 0.5 * ( alpha - gamma ) );

    EvtComplex m11( cb2 * capg2, -cb2 * sapg2 );
    EvtComplex m12( -sb2 * camg2, sb2 * samg2 );
    EvtComplex m21( sb2 * camg2, sb2 * samg2 );
    EvtComplex m22( cb2 * capg2, cb2 * sapg2 );

    retVal[0] = m11 * spinor[0] + m12 * spinor[1];
    retVal[1] = m21 * spinor[0] + m22 * spinor[1];
    retVal[2] = m11 * spinor[2] + m12 * spinor[3];
    retVal[3] = m21 * spinor[2] + m22 * spinor[3];

    spinor[0] = retVal[0];
    spinor[1] = retVal[1];
    spinor[2] = retVal[2];
    spinor[3] = retVal[3];

    return;
}

EvtDiracSpinor EvtDiracSpinor::conj() const
{
    EvtDiracSpinor sp;

    for ( int i = 0; i < 4; i++ )
        sp.set_spinor( i, ::conj( spinor[i] ) );

    return sp;
}

EvtVector4C EvtLeptonVACurrent( const EvtDiracSpinor& d, const EvtDiracSpinor& dp )
{
  // Old code; below is a new specialized code that does it more efficiently.
  // EvtVector4C tmp;
  // tmp.set(0,d*(EvtGammaMatrix::va0()*dp));
  // tmp.set(1,d*(EvtGammaMatrix::va1()*dp));
  // tmp.set(2,d*(EvtGammaMatrix::va2()*dp));
  // tmp.set(3,d*(EvtGammaMatrix::va3()*dp));
  // return tmp;

  EvtComplex u02( real(d[0]) - real(d[2]), imag(d[2]) - imag(d[0]) );
  EvtComplex u13( real(d[1]) - real(d[3]), imag(d[3]) - imag(d[1]) );

  EvtComplex v02 = dp[2] - dp[0];
  EvtComplex v13 = dp[1] - dp[3];

  EvtComplex a = u02 * v02;
  EvtComplex b = u13 * v13;

  EvtComplex c = u02 * v13;
  EvtComplex e = u13 * v02;
  EvtComplex ce = c + e;
  EvtVector4C res( b - a, e - c, EvtComplex( -imag(ce), real(ce) ), b + a );
  //    std::cout<< res - tmp <<std::endl; // check that everything is correct
  return res;
}

EvtVector4C EvtLeptonVCurrent( const EvtDiracSpinor& d, const EvtDiracSpinor& dp )
{
    EvtVector4C temp;

    // no conjugate here; done in the multiplication
    // yes this is stupid and fooled me to for a long time (ryd)

    temp.set( 0, d * ( EvtGammaMatrix::v0() * dp ) );
    temp.set( 1, d * ( EvtGammaMatrix::v1() * dp ) );
    temp.set( 2, d * ( EvtGammaMatrix::v2() * dp ) );
    temp.set( 3, d * ( EvtGammaMatrix::v3() * dp ) );

    return temp;
/*
    double a0r = real(d[0]), a0i = imag(d[0]);
    double a1r = real(d[1]), a1i = imag(d[1]);
    double a2r = real(d[2]), a2i = imag(d[2]);
    double a3r = real(d[3]), a3i = imag(d[3]);

    double b0r = real(dp[0]), b0i = imag(dp[0]);
    double b1r = real(dp[1]), b1i = imag(dp[1]);
    double b2r = real(dp[2]), b2i = imag(dp[2]);
    double b3r = real(dp[3]), b3i = imag(dp[3]);

    double w0r = a0r*b0r + a0i*b0i, w0i = a0r*b0i - a0i*b0r;
    double w1r = a1r*b1r + a1i*b1i, w1i = a1r*b1i - a1i*b1r;
    double w2r = a2r*b2r + a2i*b2i, w2i = a2r*b2i - a2i*b2r;
    double w3r = a3r*b3r + a3i*b3i, w3i = a3r*b3i - a3i*b3r;
    double t0r = a0r*b3r + a0i*b3i, t0i = a0r*b3i - a0i*b3r;
    double t1r = a1r*b2r + a1i*b2i, t1i = a1r*b2i - a1i*b2r;
    double t2r = a2r*b1r + a2i*b1i, t2i = a2r*b1i - a2i*b1r;
    double t3r = a3r*b0r + a3i*b0i, t3i = a3r*b0i - a3i*b0r;
    double q0r = a0r*b2r + a0i*b2i, q0i = a0r*b2i - a0i*b2r;
    double q1r = a1r*b3r + a1i*b3i, q1i = a1r*b3i - a1i*b3r;
    double q2r = a2r*b0r + a2i*b0i, q2i = a2r*b0i - a2i*b0r;
    double q3r = a3r*b1r + a3i*b1i, q3i = a3r*b1i - a3i*b1r;

    EvtComplex c0(w0r + w1r + (w2r + w3r), w0i + w1i + (w2i + w3i));
    EvtComplex c1(t0r + t1r + (t2r + t3r), t0i + t1i + (t2i + t3i));
    EvtComplex c2(t0i + t2i - (t1i + t3i), t1r + t3r - (t0r + t2r));
    EvtComplex c3(q0r + q2r - (q1r + q3r), q0i + q2i - (q1i + q3i));

    return EvtVector4C(c0,c1,c2,c3);

    // EvtComplex t0 = temp.get(0), t1 = temp.get(1), t2 = temp.get(2), t3 = temp.get(3);
    // std::cout<<c0-t0<<" "<<c1-t1<<" "<<c2-t2<<" "<<c3-t3<<std::endl;
 */
}

EvtVector4C EvtLeptonACurrent( const EvtDiracSpinor& d, const EvtDiracSpinor& dp )
{
    EvtVector4C temp;

    EvtGammaMatrix mat;

    // no conjugate here; done in the multiplication
    // yes this is stupid and fooled me to for a long time (ryd)

    mat = EvtGammaMatrix::v0() - EvtGammaMatrix::va0();
    temp.set( 0, d * ( mat * dp ) );

    mat = EvtGammaMatrix::v1() - EvtGammaMatrix::va1();
    temp.set( 1, d * ( mat * dp ) );

    mat = EvtGammaMatrix::v2() - EvtGammaMatrix::va2();
    temp.set( 2, d * ( mat * dp ) );

    mat = EvtGammaMatrix::v3() - EvtGammaMatrix::va3();
    temp.set( 3, d * ( mat * dp ) );
    
    return temp;
/*
    double a0r = real(d[0]), a0i = imag(d[0]);
    double a1r = real(d[1]), a1i = imag(d[1]);
    double a2r = real(d[2]), a2i = imag(d[2]);
    double a3r = real(d[3]), a3i = imag(d[3]);

    double b0r = real(dp[0]), b0i = imag(dp[0]);
    double b1r = real(dp[1]), b1i = imag(dp[1]);
    double b2r = real(dp[2]), b2i = imag(dp[2]);
    double b3r = real(dp[3]), b3i = imag(dp[3]);

    double w0r = a0r*b0r + a0i*b0i, w0i = a0r*b0i - a0i*b0r;
    double w1r = a1r*b1r + a1i*b1i, w1i = a1r*b1i - a1i*b1r;
    double w2r = a2r*b2r + a2i*b2i, w2i = a2r*b2i - a2i*b2r;
    double w3r = a3r*b3r + a3i*b3i, w3i = a3r*b3i - a3i*b3r;

    double q0r = a0r*b2r + a0i*b2i, q0i = a0r*b2i - a0i*b2r;
    double q1r = a1r*b3r + a1i*b3i, q1i = a1r*b3i - a1i*b3r;
    double q2r = a2r*b0r + a2i*b0i, q2i = a2r*b0i - a2i*b0r;
    double q3r = a3r*b1r + a3i*b1i, q3i = a3r*b1i - a3i*b1r;

    double t0r = a0r*b1r + a0i*b1i, t0i = a0r*b1i - a0i*b1r;
    double t1r = a1r*b0r + a1i*b0i, t1i = a1r*b0i - a1i*b0r;
    double t2r = a2r*b3r + a2i*b3i, t2i = a2r*b3i - a2i*b3r;
    double t3r = a3r*b2r + a3i*b2i, t3i = a3r*b2i - a3i*b2r;

    EvtComplex c0(q0r + q2r + (q1r + q3r), q0i + q2i + (q1i + q3i));
    EvtComplex c1(t0r + t1r + (t2r + t3r), t0i + t1i + (t2i + t3i));
    EvtComplex c2(t0i + t2i - (t1i + t3i), t1r + t3r - (t0r + t2r));
    EvtComplex c3(w0r + w2r - (w1r + w3r), w0i + w2i - (w1i + w3i));

    return EvtVector4C(c0,c1,c2,c3);

    //    std::cout<<EvtVector4C(c0,c1,c2,c3)-temp<<std::endl;

    */
}

void EvtLeptonVandACurrents(EvtVector4C &v, EvtVector4C &a, const EvtDiracSpinor& x, const EvtDiracSpinor& y )
{
  EvtComplex w0 = x[0]&y[0], w1 = x[1]&y[1], w2 = x[2]&y[2], w3 = x[3]&y[3];
  EvtComplex w02 = w0 + w2, w13 = w1 + w3;
  EvtComplex W1 = w02 + w13, W2 = w02 - w13;
  EvtComplex q0 = x[0]&y[2], q1 = x[1]&y[3], q2 = x[2]&y[0], q3 = x[3]&y[1];
  EvtComplex q02 = q0 + q2, q13 = q1 + q3;
  EvtComplex Q1 = q02 + q13, Q2 = q02 - q13;
  EvtComplex e0 = x[0]&y[3], e1 = x[1]&y[2], e2 = x[2]&y[1], e3 = x[3]&y[0];
  EvtComplex e20 = e0 + e2, e13 = e1 + e3;
  EvtComplex E1 = e13 + e20, E2 = e13 - e20;
  v.set(W1, E1, EvtComplex(-imag(E2), real(E2)), Q2);
  EvtComplex t0 = x[0]&y[1], t1 = x[1]&y[0], t2 = x[2]&y[3], t3 = x[3]&y[2];
  EvtComplex t20 = t0 + t2, t13 = t1 + t3;
  EvtComplex T1 = t13 + t20, T2 = t13 - t20;
  a.set(Q1, T1, EvtComplex(-imag(T2), real(T2)), W2);
}

EvtComplex EvtLeptonSCurrent( const EvtDiracSpinor& d, const EvtDiracSpinor& dp )
{
    EvtComplex temp;

    // no conjugate here; done in the multiplication
    // yes this is stupid and fooled me to for a long time (ryd)

    temp = d * ( EvtGammaMatrix::g0() * dp );

    return temp;
/*
    double a0r = real(d[0]), a0i = imag(d[0]);
    double a1r = real(d[1]), a1i = imag(d[1]);
    double a2r = real(d[2]), a2i = imag(d[2]);
    double a3r = real(d[3]), a3i = imag(d[3]);

    double b0r = real(dp[0]), b0i = imag(dp[0]);
    double b1r = real(dp[1]), b1i = imag(dp[1]);
    double b2r = real(dp[2]), b2i = imag(dp[2]);
    double b3r = real(dp[3]), b3i = imag(dp[3]);

    EvtComplex res((a0r*b0r+a0i*b0i) + (a1r*b1r+a1i*b1i) - (a2r*b2r+a2i*b2i) + (a3r*b3r+a3i*b3i),
		   (a0r*b0i-a0i*b0r) + (a1r*b1i-a1i*b1r) - (a2r*b2i-a2i*b2r) + (a3r*b3i-a3i*b3r));
    //    std::cout<<res-temp<<std::endl;
    return res;
    */
}

EvtComplex EvtLeptonPCurrent( const EvtDiracSpinor& d, const EvtDiracSpinor& dp )
{
    EvtComplex temp;

    // no conjugate here; done in the multiplication
    // yes this is stupid and fooled me to for a long time (ryd)
    static EvtGammaMatrix m = EvtGammaMatrix::g0() * EvtGammaMatrix::g5();
    temp = d * ( m * dp );

    return temp;
/*
    double a0r = real(d[0]), a0i = imag(d[0]);
    double a1r = real(d[1]), a1i = imag(d[1]);
    double a2r = real(d[2]), a2i = imag(d[2]);
    double a3r = real(d[3]), a3i = imag(d[3]);

    double b0r = real(dp[0]), b0i = imag(dp[0]);
    double b1r = real(dp[1]), b1i = imag(dp[1]);
    double b2r = real(dp[2]), b2i = imag(dp[2]);
    double b3r = real(dp[3]), b3i = imag(dp[3]);

    EvtComplex res((a0r*b2r+a0i*b2i) + (a1r*b3r+a1i*b3i) - (a2r*b0r+a2i*b0i) - (a3r*b1r+a3i*b1i),
		   (a0r*b2i-a0i*b2r) + (a1r*b3i-a1i*b3r) - (a2r*b0i-a2i*b0r) - (a3r*b1i-a3i*b1r));
    //    std::cout<<res-temp<<std::endl;
    return res;
    */
}

EvtTensor4C EvtLeptonTCurrent( const EvtDiracSpinor& d, const EvtDiracSpinor& dp )
{
    // EvtTensor4C temp;
    // temp.zero();
    // EvtComplex i2( 0, 0.5 );

    // static EvtGammaMatrix mat01 =
    //     EvtGammaMatrix::g0() * ( EvtGammaMatrix::g0() * EvtGammaMatrix::g1() -
    //                              EvtGammaMatrix::g1() * EvtGammaMatrix::g0() );
    // static EvtGammaMatrix mat02 =
    //     EvtGammaMatrix::g0() * ( EvtGammaMatrix::g0() * EvtGammaMatrix::g2() -
    //                              EvtGammaMatrix::g2() * EvtGammaMatrix::g0() );
    // static EvtGammaMatrix mat03 =
    //     EvtGammaMatrix::g0() * ( EvtGammaMatrix::g0() * EvtGammaMatrix::g3() -
    //                              EvtGammaMatrix::g3() * EvtGammaMatrix::g0() );
    // static EvtGammaMatrix mat12 =
    //     EvtGammaMatrix::g0() * ( EvtGammaMatrix::g1() * EvtGammaMatrix::g2() -
    //                              EvtGammaMatrix::g2() * EvtGammaMatrix::g1() );
    // static EvtGammaMatrix mat13 =
    //     EvtGammaMatrix::g0() * ( EvtGammaMatrix::g1() * EvtGammaMatrix::g3() -
    //                              EvtGammaMatrix::g3() * EvtGammaMatrix::g1() );
    // static EvtGammaMatrix mat23 =
    //     EvtGammaMatrix::g0() * ( EvtGammaMatrix::g2() * EvtGammaMatrix::g3() -
    //                              EvtGammaMatrix::g3() * EvtGammaMatrix::g2() );

    // temp.set( 0, 1, i2 * ( d * ( mat01 * dp ) ) );
    // temp.set( 1, 0, -temp.get( 0, 1 ) );

    // temp.set( 0, 2, i2 * ( d * ( mat02 * dp ) ) );
    // temp.set( 2, 0, -temp.get( 0, 2 ) );

    // temp.set( 0, 3, i2 * ( d * ( mat03 * dp ) ) );
    // temp.set( 3, 0, -temp.get( 0, 3 ) );

    // temp.set( 1, 2, i2 * ( d * ( mat12 * dp ) ) );
    // temp.set( 2, 1, -temp.get( 1, 2 ) );

    // temp.set( 1, 3, i2 * ( d * ( mat13 * dp ) ) );
    // temp.set( 3, 1, -temp.get( 1, 3 ) );

    // temp.set( 2, 3, i2 * ( d * ( mat23 * dp ) ) );
    // temp.set( 3, 2, -temp.get( 2, 3 ) );

    //    return temp;

    double a0r = real(d[0]), a0i = imag(d[0]);
    double a1r = real(d[1]), a1i = imag(d[1]);
    double a2r = real(d[2]), a2i = imag(d[2]);
    double a3r = real(d[3]), a3i = imag(d[3]);

    double b0r = real(dp[0]), b0i = imag(dp[0]);
    double b1r = real(dp[1]), b1i = imag(dp[1]);
    double b2r = real(dp[2]), b2i = imag(dp[2]);
    double b3r = real(dp[3]), b3i = imag(dp[3]);

    double t0r = a0r*b3r + a0i*b3i, t0i = a0r*b3i - a0i*b3r;
    double t1r = a1r*b2r + a1i*b2i, t1i = a1r*b2i - a1i*b2r;
    double t2r = a2r*b1r + a2i*b1i, t2i = a2r*b1i - a2i*b1r;
    double t3r = a3r*b0r + a3i*b0i, t3i = a3r*b0i - a3i*b0r;

    double t01r = t2i + t3i - (t0i + t1i);
    double t01i = t0r + t1r - (t2r + t3r);

    double t02r = t0r + t3r - (t1r + t2r);
    double t02i = t0i + t3i - (t1i + t2i);

    double t03r = (a1r*b3i-a1i*b3r) + (a2r*b0i-a2i*b0r) - (a0r*b2i-a0i*b2r) - (a3r*b1i-a3i*b1r);
    double t03i = (a0r*b2r+a0i*b2i) - (a1r*b3r+a1i*b3i) - (a2r*b0r+a2i*b0i) + (a3r*b1r+a3i*b1i);

    double t12r = (a0r*b0r+a0i*b0i) - (a1r*b1r+a1i*b1i) - (a2r*b2r+a2i*b2i) + (a3r*b3r+a3i*b3i);
    double t12i = (a0r*b0i-a0i*b0r) - (a1r*b1i-a1i*b1r) - (a2r*b2i-a2i*b2r) + (a3r*b3i-a3i*b3r);

    double q0r = a0r*b1r + a0i*b1i, q0i = a0r*b1i - a0i*b1r;
    double q1r = a1r*b0r + a1i*b0i, q1i = a1r*b0i - a1i*b0r;
    double q2r = a2r*b3r + a2i*b3i, q2i = a2r*b3i - a2i*b3r;
    double q3r = a3r*b2r + a3i*b2i, q3i = a3r*b2i - a3i*b2r;

    double t13r = q1i + q2i - (q0i + q3i);
    double t13i = q0r + q3r - (q1r + q2r);

    double t23r = q0r + q1r - (q2r + q3r);
    double t23i = q0i + q1i - (q2i + q3i);

    EvtTensor4C res;
    res.set( 0, 1, EvtComplex( t01r, t01i));
    res.set( 1, 0, EvtComplex(-t01r,-t01i));

    res.set( 0, 2, EvtComplex( t02r, t02i));
    res.set( 2, 0, EvtComplex(-t02r,-t02i));

    res.set( 0, 3, EvtComplex( t03r, t03i));
    res.set( 3, 0, EvtComplex(-t03r,-t03i));

    res.set( 1, 2, EvtComplex( t12r, t12i));
    res.set( 2, 1, EvtComplex(-t12r,-t12i));

    res.set( 1, 3, EvtComplex( t13r, t13i));
    res.set( 3, 1, EvtComplex(-t13r,-t13i));

    res.set( 2, 3, EvtComplex( t23r, t23i));
    res.set( 3, 2, EvtComplex(-t23r,-t23i));

    //    std::cout<<res - temp<<std::endl;

    return res;
}

EvtTensor4C EvtLeptonTACurrent( const EvtDiracSpinor& a, const EvtDiracSpinor& b ){
  // <a|sigma_mu_nu*g5|b>
  double a0r = real(a[0]), a0i = imag(a[0]);
  double a1r = real(a[1]), a1i = imag(a[1]);
  double a2r = real(a[2]), a2i = imag(a[2]);
  double a3r = real(a[3]), a3i = imag(a[3]);
  
  double b0r = real(b[0]), b0i = imag(b[0]);
  double b1r = real(b[1]), b1i = imag(b[1]);
  double b2r = real(b[2]), b2i = imag(b[2]);
  double b3r = real(b[3]), b3i = imag(b[3]);

  double p01r = a0r*b1r + a0i*b1i, p01i = a0r*b1i - a0i*b1r;
  double p03r = a0r*b3r + a0i*b3i, p03i = a0r*b3i - a0i*b3r;
  double p10r = a1r*b0r + a1i*b0i, p10i = a1r*b0i - a1i*b0r;
  double p12r = a1r*b2r + a1i*b2i, p12i = a1r*b2i - a1i*b2r;
  double p23r = a2r*b3r + a2i*b3i, p23i = a2r*b3i - a2i*b3r;
  double p21r = a2r*b1r + a2i*b1i, p21i = a2r*b1i - a2i*b1r;
  double p32r = a3r*b2r + a3i*b2i, p32i = a3r*b2i - a3i*b2r;
  double p30r = a3r*b0r + a3i*b0i, p30i = a3r*b0i - a3i*b0r;
  
  double t01r = (p23i) + (p32i) - (p01i) - (p10i);
  double t01i = (p01r) + (p10r) - (p23r) - (p32r);
  
  double t02r = (p01r) + (p32r) - (p10r) - (p23r);
  double t02i = (p01i) + (p32i) - (p10i) - (p23i);

  double t03r = (a1r*b1i-a1i*b1r) + (a2r*b2i-a2i*b2r) - (a0r*b0i-a0i*b0r) - (a3r*b3i-a3i*b3r);
  double t03i = (a0r*b0r+a0i*b0i) + (a3r*b3r+a3i*b3i) - (a1r*b1r+a1i*b1i) - (a2r*b2r+a2i*b2i);

  double t12r = (a0r*b2r+a0i*b2i) + (a3r*b1r+a3i*b1i) - (a1r*b3r+a1i*b3i) - (a2r*b0r+a2i*b0i);
  double t12i = (a0r*b2i-a0i*b2r) + (a3r*b1i-a3i*b1r) - (a1r*b3i-a1i*b3r) - (a2r*b0i-a2i*b0r);

  double t13r = (p12i) + (p21i) - (p03i) - (p30i);
  double t13i = (p03r) + (p30r) - (p12r) - (p21r);  

  double t23r = (p03r) + (p12r) - (p21r) - (p30r);
  double t23i = (p03i) + (p12i) - (p21i) - (p30i);

  EvtTensor4C res;

  res.set( 0, 1, EvtComplex( t01r, t01i));
  res.set( 1, 0, EvtComplex(-t01r,-t01i));

  res.set( 0, 2, EvtComplex( t02r, t02i));
  res.set( 2, 0, EvtComplex(-t02r,-t02i));

  res.set( 0, 3, EvtComplex( t03r, t03i));
  res.set( 3, 0, EvtComplex(-t03r,-t03i));

  res.set( 1, 2, EvtComplex( t12r, t12i));
  res.set( 2, 1, EvtComplex(-t12r,-t12i));

  res.set( 1, 3, EvtComplex( t13r, t13i));
  res.set( 3, 1, EvtComplex(-t13r,-t13i));

  res.set( 2, 3, EvtComplex( t23r, t23i));
  res.set( 3, 2, EvtComplex(-t23r,-t23i));

  return res;
}

EvtDiracSpinor operator*( const EvtComplex& c, const EvtDiracSpinor& d )
{
    // EvtDiracSpinor result;
    // result.spinor[0] = c * d.spinor[0];
    // result.spinor[1] = c * d.spinor[1];
    // result.spinor[2] = c * d.spinor[2];
    // result.spinor[3] = c * d.spinor[3];
    // return result;
  return EvtDiracSpinor(c*d[0], c*d[1], c*d[2], c*d[3]);
}

EvtDiracSpinor EvtDiracSpinor::adjoint() const
{
    // EvtDiracSpinor d = this->conj();    // first conjugate, then multiply with gamma0
    // EvtGammaMatrix g0 = EvtGammaMatrix::g0();
    // EvtDiracSpinor result;    // automatically initialized to 0

    // for ( int i = 0; i < 4; ++i )
    //     for ( int j = 0; j < 4; ++j )
    //         result.spinor[i] += d.spinor[j] * g0._gamma[i][j];

    // return result;
  const EvtDiracSpinor &d = *this;
  double a0r = real(d[0]), a0i = imag(d[0]);
  double a1r = real(d[1]), a1i = imag(d[1]);
  double a2r = real(d[2]), a2i = imag(d[2]);
  double a3r = real(d[3]), a3i = imag(d[3]);
  return EvtDiracSpinor(EvtComplex(a0r,-a0i),EvtComplex(a1r,-a1i),EvtComplex(-a2r,a2i),EvtComplex(-a3r,a3i));
}

EvtComplex operator*( const EvtDiracSpinor& d, const EvtDiracSpinor& dp )
{
  // int i;
  // EvtComplex temp;

  // temp = EvtComplex( 0.0, 0.0 );

  // for ( i = 0; i < 4; i++ ) {
  //     temp += conj( d.get_spinor( i ) ) * dp.get_spinor( i );
  // }
  // return temp;

  double a0r = real(d[0]), a0i = imag(d[0]);
  double a1r = real(d[1]), a1i = imag(d[1]);
  double a2r = real(d[2]), a2i = imag(d[2]);
  double a3r = real(d[3]), a3i = imag(d[3]);

  double b0r = real(dp[0]), b0i = imag(dp[0]);
  double b1r = real(dp[1]), b1i = imag(dp[1]);
  double b2r = real(dp[2]), b2i = imag(dp[2]);
  double b3r = real(dp[3]), b3i = imag(dp[3]);

  double w0r = a0r*b0r + a0i*b0i, w0i = a0r*b0i - a0i*b0r;
  double w1r = a1r*b1r + a1i*b1i, w1i = a1r*b1i - a1i*b1r;
  double w2r = a2r*b2r + a2i*b2i, w2i = a2r*b2i - a2i*b2r;
  double w3r = a3r*b3r + a3i*b3i, w3i = a3r*b3i - a3i*b3r;

  return EvtComplex(w0r+w1r+(w2r+w3r), w0i+w1i+(w2i+w3i));
}

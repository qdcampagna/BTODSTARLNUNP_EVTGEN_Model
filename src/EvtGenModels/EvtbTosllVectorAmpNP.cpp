
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

#include "EvtGenModels/EvtbTosllVectorAmpNP.hh"

#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenModels/EvtbTosllAmp.hh"
#include "EvtGenModels/EvtbTosllFF.hh"

#if 0
#include "EvtGenBase/EvtRandom.hh"
void TestCalcAmp(EvtParticle*, EvtAmp&, EvtbTosllFF*);
void Rusa       (EvtParticle*, EvtAmp&, EvtbTosllFF*);
void Rusa2      (EvtParticle*, EvtAmp&, EvtbTosllFF*);
#endif

static EvtId IdB0,  IdaB0, IdBp,  IdBm, IdBs,  IdaBs, IdRhop, IdRhom, IdRho0, IdOmega, IdKst0, IdaKst0, IdKstp, IdKstm;
static bool cafirst = true;
void EvtbTosllVectorAmpNP::CalcAmp( EvtParticle* parent, EvtAmp& amp, EvtbTosllFF* formFactors ){
  if(cafirst){
    cafirst = false;
    IdB0 = EvtPDL::getId("B0");  IdaB0 = EvtPDL::getId("anti-B0");
    IdBp = EvtPDL::getId("B+");  IdBm = EvtPDL::getId("B-");
    IdBs = EvtPDL::getId("B_s0");  IdaBs = EvtPDL::getId("anti-B_s0");
    IdRhop = EvtPDL::getId("rho+"); IdRhom =  EvtPDL::getId("rho-"); IdRho0 =  EvtPDL::getId("rho0");
    IdOmega =  EvtPDL::getId("omega"); IdKst0 =  EvtPDL::getId("K*0" ); IdaKst0 =  EvtPDL::getId("anti-K*0" );
    IdKstp = EvtPDL::getId( "K*+" ); IdKstm = EvtPDL::getId( "K*-" );
  }

  // Add the lepton and neutrino 4 momenta to find q2
  EvtId pId = parent->getId();
  EvtId dId = parent->getDaug( 0 )->getId();
  EvtVector4R q = parent->getDaug( 1 )->getP4() + parent->getDaug( 2 )->getP4();
  double q2 = q.mass2();
  double mesonmass = parent->getDaug( 0 )->mass();

  double a1, a2, a0, v, t1, t2, t3; // form factors
  formFactors->getVectorFF( pId, dId, q2, mesonmass, a1, a2, a0, v, t1, t2, t3 );

  bool btod = false;
  bool nnlo = true;
  if ( ( pId == IdB0 || pId == IdaB0 || pId == IdBp || pId == IdBm ) &&
       ( dId == IdRhop || dId == IdRhom || dId == IdRho0 || dId ==  IdOmega ) ) {
    btod = true;
  }
  if ( ( pId == IdBs || pId == IdaBs ) &&
       ( dId == IdKst0 || dId == IdaKst0 || dId == IdKstp || dId == IdKstm ) ) {
    btod = true;
  }

  const EvtVector4R &p4b = parent->getP4();//( pmass, 0.0, 0.0, 0.0 );
  const EvtVector4R &p4meson = parent->getDaug( 0 )->getP4();

  double pmass = parent->mass(), ipmass = 1/pmass;
  EvtVector4R pbhat = p4b * ipmass;
  EvtVector4R qhat = q * ipmass;
  EvtVector4R pkstarhat = p4meson * ipmass;
  EvtVector4R phat = pbhat + pkstarhat;

  EvtComplex c7  = EvtbTosllAmp::GetC7Eff ( q2, nnlo ); // c7eff
  EvtComplex c9  = EvtbTosllAmp::GetC9Eff ( q2, nnlo, btod ); // c9eff
  EvtComplex c10 = EvtbTosllAmp::GetC10Eff( q2, nnlo ); // c10eff
  EvtComplex I( 0.0, 1.0 );

  double mb = 4.18/*GeV/c^2*/*ipmass, ms = 0.093/*GeV/c^2*/*ipmass;
  double mH = mesonmass * ipmass, oamH = 1 + mH, osmH = 1 - mH, osmH2 = oamH*osmH, iosmH2 = 1/osmH2; // mhatkstar
  double is = pmass * pmass / q2; // 1/shat
  a1 *= oamH;
  a2 *= osmH;
  a0 *= 2*mH;
  double cs0 = ( a1 - a2 - a0 ) * is;
  a2 *= iosmH2;
  v *= 2/oamH;

  EvtComplex a = (c9 + _c9p) * v   + (c7 + _c7p) * (4 * mb * is * t1);
  EvtComplex b = (c9 - _c9p) * a1  + (c7 - _c7p) * (2 * mb * is * osmH2 * t2);
  EvtComplex c = (c9 - _c9p) * a2  + (c7 - _c7p) * (2 * mb * ( t3 * iosmH2 + t2 * is ) );
  EvtComplex d = (c9 - _c9p) * cs0 - (c7 - _c7p) * (2 * mb * is * t3);
  EvtComplex e = (c10 + _c10p) * v;
  EvtComplex f = (c10 - _c10p) * a1;
  EvtComplex g = (c10 - _c10p) * a2;
  EvtComplex h = (c10 - _c10p) * cs0;
  double sscale = a0/(mb + ms);
  EvtComplex hs = _cS*sscale, hp = _cP*sscale;

  int charge1 = EvtPDL::chg3( parent->getDaug( 1 )->getId() );

  EvtParticle* lepPos = ( charge1 > 0 ) ? parent->getDaug( 1 ) : parent->getDaug( 2 );
  EvtParticle* lepNeg = ( charge1 < 0 ) ? parent->getDaug( 1 ) : parent->getDaug( 2 );

  EvtDiracSpinor lp0(lepPos->spParent( 0 )), lp1(lepPos->spParent( 1 ));
  EvtDiracSpinor lm0(lepNeg->spParent( 0 )), lm1(lepNeg->spParent( 1 ));

  EvtVector4C l11, l12, l21, l22, a11, a12, a21, a22;
  EvtComplex s11, s12, s21, s22, p11, p12, p21, p22;
  EvtTensor4C tt0(EvtGenFunctions::asymProd( pbhat, pkstarhat ) );

  EvtTensor4C T1(tt0), T2(tt0);
  const EvtTensor4C &gmn = EvtTensor4C::g();
  EvtTensor4C tt1(EvtGenFunctions::directProd( pbhat, phat ));
  EvtTensor4C tt2(EvtGenFunctions::directProd( pbhat, qhat ));

  b *= I; c *= I; d *= I; f *= I; g *= I; h *= I;
  if ( pId == IdBm || pId == IdaB0 || pId == IdaBs ) {
    T1 *= a; T1.addScaled(-b, gmn); T1.addScaled(c, tt1); T1.addScaled(d, tt2);
    T2 *= e; T2.addScaled(-f, gmn); T2.addScaled(g, tt1); T2.addScaled(h, tt2);

    EvtLeptonVandACurrents(l11, a11, lp0, lm0);
    EvtLeptonVandACurrents(l21, a21, lp1, lm0);
    EvtLeptonVandACurrents(l12, a12, lp0, lm1);
    EvtLeptonVandACurrents(l22, a22, lp1, lm1);

    s11 = EvtLeptonSCurrent(lp0, lm0); p11 = EvtLeptonPCurrent(lp0, lm0);
    s21 = EvtLeptonSCurrent(lp1, lm0); p21 = EvtLeptonPCurrent(lp1, lm0);
    s12 = EvtLeptonSCurrent(lp0, lm1); p12 = EvtLeptonPCurrent(lp0, lm1);
    s22 = EvtLeptonSCurrent(lp1, lm1); p22 = EvtLeptonPCurrent(lp1, lm1);
  } else if ( pId == IdBp || pId == IdB0 || pId == IdBs ) {
    T1 *= -a; T1.addScaled(-b, gmn); T1.addScaled(c, tt1); T1.addScaled(d, tt2);
    T2 *= -e; T2.addScaled(-f, gmn); T2.addScaled(g, tt1); T2.addScaled(h, tt2);

    EvtLeptonVandACurrents(l11, a11, lp1, lm1);
    EvtLeptonVandACurrents(l21, a21, lp0, lm1);
    EvtLeptonVandACurrents(l12, a12, lp1, lm0);
    EvtLeptonVandACurrents(l22, a22, lp0, lm0);

    s11 = EvtLeptonSCurrent(lp1, lm1); p11 = EvtLeptonPCurrent(lp1, lm1);
    s21 = EvtLeptonSCurrent(lp0, lm1); p21 = EvtLeptonPCurrent(lp0, lm1);
    s12 = EvtLeptonSCurrent(lp1, lm0); p12 = EvtLeptonPCurrent(lp1, lm0);
    s22 = EvtLeptonSCurrent(lp0, lm0); p22 = EvtLeptonPCurrent(lp0, lm0);
  } else {
    EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Wrong lepton number\n";
    T1.zero();
    T2.zero();    // Set all tensor terms to zero.
  }

  for (int i = 0; i < 3; i++ ) {
    EvtVector4C eps = parent->getDaug( 0 )->epsParent( i ).conj();

    EvtVector4C E1 = T1.cont1( eps );
    EvtVector4C E2 = T2.cont1( eps );

    EvtComplex epsq = I*(eps*q), epsqs =  epsq * hs, epsqp = epsq * hp;

    amp.vertex( i, 0, 0, l11 * E1 + a11 * E2 - s11 * epsqs - p11 * epsqp );
    amp.vertex( i, 0, 1, l12 * E1 + a12 * E2 - s12 * epsqs - p12 * epsqp );
    amp.vertex( i, 1, 0, l21 * E1 + a21 * E2 - s21 * epsqs - p21 * epsqp );
    amp.vertex( i, 1, 1, l22 * E1 + a22 * E2 - s22 * epsqs - p22 * epsqp );
  }
}

#if 0
static EvtComplex randbuf[10];
static bool firstrun = true;
void TestCalcAmp( EvtParticle* parent, EvtAmp& amp, EvtbTosllFF* formFactors ){
  if(firstrun){
    for(int i=0;i<10;i++) randbuf[i] = EvtComplex(EvtRandom::random(), EvtRandom::random());
  }
  EvtId pId = parent->getId();
  EvtId dId = parent->getDaug( 0 )->getId();
  EvtVector4R q = parent->getDaug( 1 )->getP4() + parent->getDaug( 2 )->getP4();
  double q2 = q.mass2();
  double mesonmass = parent->getDaug( 0 )->mass();

  double a1, a2, a0, v, t1, t2, t3; // form factors
  formFactors->getVectorFF( pId, dId, q2, mesonmass, a1, a2, a0, v, t1, t2, t3 );

  const EvtVector4R &p4b = parent->getP4();//( pmass, 0.0, 0.0, 0.0 );
  const EvtVector4R &p4meson = parent->getDaug( 0 )->getP4();
  std::cout<<"Test "<<p4meson<<std::endl;

  double pmass = parent->mass(), ipmass = 1/pmass;
  EvtVector4R pbhat = p4b * ipmass;
  EvtVector4R qhat = q * ipmass;
  EvtVector4R pkstarhat = p4meson * ipmass;
  EvtVector4R phat = pbhat + pkstarhat;

  EvtComplex c7 = randbuf[0];
  EvtComplex c9 = randbuf[1];
  EvtComplex c10 = randbuf[2];
  EvtComplex I( 0.0, 1.0 );

  double mb = 4.18 * ipmass; // mhatb
  double mH = mesonmass * ipmass, oamH = 1 + mH, osmH = 1 - mH, osmH2 = oamH*osmH, iosmH2 = 1/osmH2; // mhatkstar
  double is = pmass * pmass / q2; // 1/shat
  a1 *= oamH;
  a2 *= osmH;
  double cs0 = ( a1 - a2 - 2 * mH * a0 ) * is;
  a2 *= iosmH2;
  v *= 2/oamH;

  double ms = 0.093 * ipmass;
  EvtComplex c7p = randbuf[3];
  EvtComplex c9p = randbuf[4];
  EvtComplex c10p = randbuf[5];

  EvtComplex a = (c9 + c9p) * v   + (c7 + c7p) * (4 * mb * is * t1);
  EvtComplex b = (c9 - c9p) * a1  + (c7 - c7p) * (2 * mb * is * osmH2 * t2);
  EvtComplex c = (c9 - c9p) * a2  + (c7 - c7p) * (2 * mb * ( t3 * iosmH2 + t2 * is ) );
  EvtComplex d = (c9 - c9p) * cs0 - (c7 - c7p) * (2 * mb * is * t3);
  EvtComplex e = (c10 + c10p) * v;
  EvtComplex f = (c10 - c10p) * a1;
  EvtComplex g = (c10 - c10p) * a2;
  EvtComplex h = (c10 - c10p) * cs0;

  EvtComplex cS(randbuf[6]), cSp(randbuf[7]), cP(randbuf[8]), cPp(randbuf[9]);

  EvtTensor4C tt0(EvtGenFunctions::asymProd( pbhat, pkstarhat ) );
  EvtTensor4C T1(tt0), T2(tt0);
  const EvtTensor4C &gmn = EvtTensor4C::g();
  EvtTensor4C tt1(EvtGenFunctions::directProd( pbhat, phat ));
  EvtTensor4C tt2(EvtGenFunctions::directProd( pbhat, qhat ));

  b *= I; c *= I; d *= I; f *= I; g *= I; h *= I;
  T1 *= a; T1.addScaled(-b, gmn); T1.addScaled(c, tt1); T1.addScaled(d, tt2);
  T2 *= e; T2.addScaled(-f, gmn); T2.addScaled(g, tt1); T2.addScaled(h, tt2);

  //  std::cout<<cS<<" "<<cSp<<" "<<mH<<" "<<a0<<" "<<mb<<" "<<ms<<qhat<<std::endl;
  for (int i = 0; i < 3; i++ ) {
    EvtVector4C eps = parent->getDaug( 0 )->epsParent( i ).conj();

    EvtVector4C E1 = T1.cont1( eps );
    EvtVector4C E2 = T2.cont1( eps );
    EvtComplex epsq = eps*qhat;
    //    std::cout<<i<<" "<<epsq<<std::endl;
    EvtComplex S = -I*(cS-cSp)*epsq*2*mH*a0/(mb + ms), P = -I*(cP-cPp)*epsq*2*mH*a0/(mb + ms);
    std::cout<<E1<<" "<<E2<<" "<<S<<" "<<P<<std::endl;
  }
}


using EvtGenFunctions::asymProd;
void Rusa(EvtParticle* parent, EvtAmp& amp, EvtbTosllFF* formFactors){
  EvtId pId = parent->getId();
  const EvtVector4R &p = parent->getP4();
  const double mB = p.mass();

  const EvtVector4R q = parent->getDaug(1)->getP4() + parent->getDaug(2)->getP4();
  const double q2 = q.mass2();

  EvtId dId = parent->getDaug(0)->getId();
  const EvtVector4R &k = parent->getDaug(0)->getP4();
  std::cout<<"Rusa "<<k<<std::endl;
  const double mK = k.mass();

  const double mb = 4.18, ms = 0.093;

  double A1, A2, A0, V, T1, T2, T3; // form factors
  formFactors->getVectorFF( pId, dId, q2, mK, A1, A2, A0, V, T1, T2, T3 );
  double A3 = (mB + mK)/(2*mK)*A1 - (mB - mK)/(2*mK)*A2;

  EvtComplex C7(randbuf[0]), C9(randbuf[1]), C10(randbuf[2]);
  EvtComplex C7p(randbuf[3]), C9p(randbuf[4]), C10p(randbuf[5]);
  EvtComplex CS(randbuf[6]), CSp(randbuf[7]), CP(randbuf[8]), CPp(randbuf[9]);
  for (int j = 0; j < 3; j++ ) {
    EvtComplex i(0.0, 1.0);
    EvtVector4C eps = parent->getDaug(0)->epsParent(j).conj();

    EvtVector4C H1Lmu = -i*eps*((mB+mK)*A1) + i*(2*p-q)*((eps*q)*A2/(mB+mK))
      + i*q*((eps*q)*2*mK/q2*(A3-A0)) + asymProd(eps, p, k)*(2*V/(mB+mK));
    EvtVector4C H1Rmu =  i*eps*((mB+mK)*A1) - i*(2*p-q)*((eps*q)*A2/(mB+mK))
      - i*q*((eps*q)*2*mK/q2*(A3-A0)) + asymProd(eps, p, k)*(2*V/(mB+mK));

    EvtVector4C H2Lmu = i*asymProd(eps, p, k)*(2*T1)
      - T2*(eps*(mB*mB - mK*mK) - (eps*q)*(2*p-q)) - T3*(eps*q)*(q - q2/(mB*mB-mK*mK)*(2*p-q));
    EvtVector4C H2Rmu = i*asymProd(eps, p, k)*(2*T1)
      + T2*(eps*(mB*mB - mK*mK) - (eps*q)*(2*p-q)) + T3*(eps*q)*(q - q2/(mB*mB-mK*mK)*(2*p-q));

    EvtComplex H3L =  i*(eps*q)*(2*mK/(mb + ms)*A0);
    EvtComplex H3R = -i*(eps*q)*(2*mK/(mb + ms)*A0);

    EvtVector4C Vmu = C9*H1Lmu + C9p*H1Rmu - (2*mb/q2)*(i*C7*H2Rmu + i*C7p*H2Lmu);
    EvtVector4C Amu = C10*H1Lmu + C10p*H1Rmu;
    EvtComplex S = CS*H3R + CSp*H3L, P = CP*H3R + CPp*H3L;

    //    std::cout<<j<<" "<<eps<<" "<<eps*q<<std::endl;
    std::cout<<Vmu*(1/mB)<<" "<<Amu*(1/mB)<<" "<<S*(1/mB)<<" "<<P*(1/mB)<<std::endl;

  }


  // for(int i=0; i<5; i++){
  //   EvtVector4R a, b, c;
  //   a.set(EvtRandom::random(), EvtRandom::random(), EvtRandom::random(), EvtRandom::random());
  //   b.set(EvtRandom::random(), EvtRandom::random(), EvtRandom::random(), EvtRandom::random());
  //   c.set(EvtRandom::random(), EvtRandom::random(), EvtRandom::random(), EvtRandom::random());
  //   std::cout<<asymProd(a, b, c)<<std::endl;
  // }
  // exit(0);
}

void Rusa2(EvtParticle* parent, EvtAmp& amp, EvtbTosllFF* formFactors){
  EvtId pId = parent->getId();
  EvtVector4R p = parent->getP4();
  const double mB = p.mass(), imB = 1/mB;
  p *= imB;

  const EvtVector4R q = (parent->getDaug(1)->getP4() + parent->getDaug(2)->getP4())*imB;
  const double q2 = q.mass2();

  EvtId dId = parent->getDaug(0)->getId();
  const EvtVector4R k = parent->getDaug(0)->getP4() * imB;
  std::cout<<"Rusa2 "<<k<<std::endl;
  const double mK = k.mass();

  const double mb = 4.18*imB, ms = 0.093*imB;

  double A1, A2, A0, V, T1, T2, T3; // form factors
  formFactors->getVectorFF( pId, dId, q2*mB*mB, mK*mB, A1, A2, A0, V, T1, T2, T3 );

  EvtTensor4C tt0(EvtGenFunctions::asymProd( p, k ) );
  const EvtTensor4C &gmn = EvtTensor4C::g();
  EvtTensor4C tt1(EvtGenFunctions::directProd( p, p + k ));
  EvtTensor4C tt2(EvtGenFunctions::directProd( p, q ));

  EvtComplex C7(randbuf[0]), C9(randbuf[1]), C10(randbuf[2]);
  EvtComplex C7p(randbuf[3]), C9p(randbuf[4]), C10p(randbuf[5]);
  EvtComplex CS(randbuf[6]), CSp(randbuf[7]), CP(randbuf[8]), CPp(randbuf[9]);
  //  std::cout<<CS<<" "<<CSp<<" "<<mK<<" "<<A0<<" "<<mb<<" "<<ms<<q<<std::endl;
  for (int j = 0; j < 3; j++ ) {
    EvtComplex i(0.0, 1.0);
    EvtVector4C eps = parent->getDaug(0)->epsParent(j).conj();

    EvtVector4C Vmu =
      -((C9+C9p)*2*V/(1+mK) + (C7+C7p)*4*mb/q2*T1)*tt0.cont2(eps)
      - i*((C9-C9p)*(1+mK)*A1  + (C7-C7p)*2*mb/q2*(1-mK*mK)*T2)*gmn.cont1(eps)
      + i*((C9-C9p)*A2/(1+mK)  + (C7-C7p)*2*mb*(T2/q2 + T3/(1-mK*mK)))*tt1.cont1(eps)
      + i*((C9-C9p)*(((1+mK)*A1 - (1-mK)*A2 - 2*mK*A0)/q2) - (C7-C7p)*2*mb/q2*T3)*tt2.cont1(eps);

    //    std::cout<<tt0.cont1(eps)<<" "<<tt0.cont2(eps)<<std::endl;

    EvtVector4C Amu =
      -(C10+C10p)*2*V/(1+mK)*tt0.cont2(eps)
      - i*(C10-C10p)*(1+mK)*A1*gmn.cont1(eps)
      + i*(C10-C10p)*A2/(1+mK)*tt1.cont1(eps)
      + i*(C10-C10p)*(((1+mK)*A1 - (1-mK)*A2 - 2*mK*A0)/q2)*tt2.cont1(eps);

    //    EvtComplex S = i*(CS-CSp)*2*mK*A0/(mb+ms)*(eps*q), P = i*(CP-CPp)*2*mK*A0/(mb+ms)*(eps*q);
    EvtComplex S = -i*(CS-CSp)*2*mK*A0/(mb + ms)*(eps*q), P = -i*(CP-CPp)*2*mK*A0/(mb + ms)*(eps*q);

    //    std::cout<<j<<" "<<eps<<" "<<eps*q<<std::endl;
    std::cout<<Vmu<<" "<<Amu<<" "<<S<<" "<<P<<std::endl;

  }
}
#endif

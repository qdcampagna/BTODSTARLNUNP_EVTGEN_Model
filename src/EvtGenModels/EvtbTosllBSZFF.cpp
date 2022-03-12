#include "EvtGenModels/EvtbTosllBSZFF.hh"
#include "EvtGenBase/EvtPDL.hh"
#include <cmath>

inline double poly(double x, int n, const double *c){
  double t = c[--n];
  while(n--) t = c[n] + x*t;
  return t;
}

void EvtbTosllBSZFF::getVectorFF(EvtId pId, EvtId, double q2, double mV,
				 double& a1, double& a2, double& a0, double& v, double& t1,  double& t2, double& t3 ){
  const static double alfa[7][4] =
    {
      // coefficients are from  https://arxiv.org/src/1503.05534v3/anc/BKstar_LCSR-Lattice.json
      //  m_res,       c0,        c1,       c2
      {5.415000, 0.376313, -1.165970, 2.424430}, // V
      {5.366000, 0.369196, -1.365840, 0.128191}, // A0
      {5.829000, 0.297250,  0.392378, 1.189160}, // A1
      {5.829000, 0.265375,  0.533638, 0.483166}, // A12
      {5.415000, 0.312055, -1.008930, 1.527200}, // T1
      {5.829000, 0.312055,  0.496846, 1.614310}, // T2
      {5.829000, 0.667412,  1.318120, 3.823340}  // T12
    };

  double mB = EvtPDL::getMeanMass( pId );
  double mBaV = mB + mV, mBsV = mB - mV;
  double tp = mBaV*mBaV; // t_{+} = (m_B + m_V)^2
  double s0 = sqrt(2*mBaV*sqrt(mB*mV)); // sqrt(t_{+} - t_{+}*(1 - sqrt(1 - t_{-}/t_{+})))
  double z0 = (mBaV - s0)/(mBaV + s0); // (sqrt(t_{+}) - s0)/(sqrt(t_{+}) + s0)

  double s = sqrt(tp - q2), z = (s - s0)/(s + s0), dz = z - z0, ff[7];
  for(int j=0; j<7; j++){
    double mR = alfa[j][0], mR2 = mR*mR;
    ff[j] = (mR2/(mR2 - q2)) * poly(dz, 3, alfa[j] + 1);
  }

  // Källén-function
  // arXiv:1503.05534 Eq. D.3
  double lambda = (mBaV*mBaV - q2) * (mBsV*mBsV - q2);

  v = ff[0];
  a0 = ff[1];
  a1 = ff[2];
  // Eq. D.5 arXiv:1503.05534
  // A12 = (mBaV*mBaV*(mBaV*mBsV - q2)*A1 - lambda(q2)*A2)/(16*mB*mV*mV*mBaV);
  double a12 = ff[3];
  a2 = mBaV*((mBaV*(mBaV*mBsV - q2))*a1 - (16*mB*mV*mV)*a12)/lambda;
  t1 = ff[4];
  t2 = ff[5];
  // Eq. D.5 arXiv:1503.05534
  // T23 = mBaV*mBsV*(mB*mB + 3*mV*mV - q2)*T2 - lambda(q2)*T3)/(8*mB*mV*mV*mBsV);
  double t23 = ff[6];
  t3 = mBsV*(mBaV*(mB*mB + 3*mV*mV - q2)*t2 - (8*mB*mV*mV)*t23)/lambda;
}

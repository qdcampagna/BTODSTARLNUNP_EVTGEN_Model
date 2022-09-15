
/**************************************************************************
 * Written by Alexei Sibidanov                                            *
 * Modified & Updated by Lopamudra Mukherjee & Quinn Campagna             *
**************************************************************************/

/**************************************************************************
 * Version date : 17th Feb, 2022                                          *
 * alpha_s/pi, 1/mb, 1/mc and 1/mc^2 corrections included in the HQET FFs *
**************************************************************************/

/**************************************************************************
 * Version date : 11th April, 2022                                        *
 * BGL FF parametrization included *
**************************************************************************/

#include "EvtGenModels/EvtSemiLeptonicVectorAmpNP.hh"

#include <cmath>
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtDiLog.hh"

EvtComplex asymProd(const EvtTensor4C &t1, const EvtTensor4C &t2){
  //returns  eps_ijkl*t1(i,j)*t2(k,l);
  // eps_ijkl is the antisymmetric Levi-Civita tensor
  EvtComplex sum;
  /*sum -= t1.get(0,1)*t2.get(2,3);
  sum += t1.get(0,1)*t2.get(3,2);
  sum += t1.get(0,2)*t2.get(1,3);
  sum -= t1.get(0,2)*t2.get(3,1);
  sum -= t1.get(0,3)*t2.get(1,2);
  sum += t1.get(0,3)*t2.get(2,1);
  sum += t1.get(1,0)*t2.get(2,3);
  sum -= t1.get(1,0)*t2.get(3,2);
  sum -= t1.get(1,2)*t2.get(0,3);
  sum += t1.get(1,2)*t2.get(3,0);
  sum += t1.get(1,3)*t2.get(0,2);
  sum -= t1.get(1,3)*t2.get(2,0);
  sum -= t1.get(2,0)*t2.get(1,3);
  sum += t1.get(2,0)*t2.get(3,1);
  sum += t1.get(2,1)*t2.get(0,3);
  sum -= t1.get(2,1)*t2.get(3,0);
  sum -= t1.get(2,3)*t2.get(0,1);
  sum += t1.get(2,3)*t2.get(1,0);
  sum += t1.get(3,0)*t2.get(1,2);
  sum -= t1.get(3,0)*t2.get(2,1);
  sum -= t1.get(3,1)*t2.get(0,2);
  sum += t1.get(3,1)*t2.get(2,0);
  sum += t1.get(3,2)*t2.get(0,1);
  sum -= t1.get(3,2)*t2.get(1,0);
  return sum;*/
  
  sum += t1.get(0,1)*t2.get(2,3);
  sum -= t1.get(0,1)*t2.get(3,2);
  sum -= t1.get(0,2)*t2.get(1,3);
  sum += t1.get(0,2)*t2.get(3,1);
  sum += t1.get(0,3)*t2.get(1,2);
  sum -= t1.get(0,3)*t2.get(2,1);
  sum -= t1.get(1,0)*t2.get(2,3);
  sum += t1.get(1,0)*t2.get(3,2);
  sum += t1.get(1,2)*t2.get(0,3);
  sum -= t1.get(1,2)*t2.get(3,0);
  sum -= t1.get(1,3)*t2.get(0,2);
  sum += t1.get(1,3)*t2.get(2,0);
  sum += t1.get(2,0)*t2.get(1,3);
  sum -= t1.get(2,0)*t2.get(3,1);
  sum -= t1.get(2,1)*t2.get(0,3);
  sum += t1.get(2,1)*t2.get(3,0);
  sum += t1.get(2,3)*t2.get(0,1);
  sum -= t1.get(2,3)*t2.get(1,0);
  sum -= t1.get(3,0)*t2.get(1,2);
  sum += t1.get(3,0)*t2.get(2,1);
  sum += t1.get(3,1)*t2.get(0,2);
  sum -= t1.get(3,1)*t2.get(2,0);
  sum -= t1.get(3,2)*t2.get(0,1);
  sum += t1.get(3,2)*t2.get(1,0);
  return sum;
}

void getFF(double mB, double mV, double q2,
	   double &V, double &A0, double &A1, double &A2, double &T1, double &T2, double &T3 ){ 
	     
  /* Quark masses from PDG */
  const double m_c = 1.27, // +- 0.02 GeV PDG
    m_b = 4.18; // +- 0.025 GeV PDG
    
  double rmB = 1/mB, rDstar = mV/mB, r2Dstar = rDstar*rDstar;
  double c = 0.5/sqrt(mB*mV), mBaV = mB + mV, mBsV = mB - mV, rmBaV = 1/mBaV, rmBsV = 1/mBsV;
  double w = (mB*mB + mV*mV - q2)/(2*mB*mV), wsq = w*w, sqrtwa1 = sqrt(w+1), ws1 = w-1, ws12 = ws1*ws1;
  double u = 1 + r2Dstar - 2*rDstar*w;
  double z = (sqrtwa1 - sqrt(2))/(sqrtwa1 + sqrt(2));


/************************************************************************************/
/* // CLN Parametrization of the FFs from 1309.0301 begin here :
  
   // CLN FF Parameters based on arXiv:1309.0301
    const double
    rho2Dst = 1.207, // +- 0.026 Eq 44
    R1t = 1.403, // +- 0.033 Eq 44
    R2t = 0.854, // +- 0.020 Eq 44
    hA1t = 0.908; // +- 0.017 Eq 45;
    
  //  double V1w = V1t*(1 - z*(8*rho2D - z*((51*rho2D - 10) - z*(252*rho2D - 84))));
  double hA1w = hA1t*(1 - z*(8*rho2Dst - z*((53*rho2Dst - 15) - z*(231*rho2Dst - 91))));
  double R1w = R1t - 0.12*ws1 + 0.05*ws12;
  double R2w = R2t + 0.11*ws1 - 0.06*ws12;
  double R3w = 1.22 - 0.052*ws1 + 0.026*ws12;
  double hVw = R1w*hA1w;
  double hA2w = (R2w - R3w)/(2*rDstar)*hA1w;
  double hA3w = (R2w + R3w)*0.5*hA1w;
  
 // arXiv:1309.0301 Eq 48b
  double cA1 = (m_b-m_c)*rmBsV*hA1w, cV = (m_b+m_c)*rmBaV*hVw;
  double hT1 = (0.5*((1-rDstar)*(1-rDstar)*(w+1)*cA1
		  -(1+rDstar)*(1+rDstar)*(w-1)*cV))/u;
  double hT2 = (0.5*(1-r2Dstar)*(w+1)*(cA1 - cV))/u;
  double hT3 = (-0.5*(2*rDstar*(w+1)*cA1
				+(m_b-m_c)*rmBsV*u*(hA3w-rDstar*hA2w)
				-(1+rDstar)*(1+rDstar)*cV))/(u*(1+rDstar)); */

// CLN Parametrization ends here.
/************************************************************************************/


/************************************************************************************/   
// HQET Parametrization begins here : 
 
 /* // HQET FF parameters based on Table 2 Column 1 (SM 3/2/1) arxiv:2004.10208
   const double 
   Xi0 = 1,
   Xi1 = -0.93, // +- 0.10
   Xi2 = 1.35, // +- 0.26
   Xi3 = -2.67, // +- 0.75
   chi20 = -0.05, // +- 0.02
   chi21 = 0.01, // +- 0.02
   chi22 = -0.01, // +- 0.02
   chi30 = 0,
   chi31 = -0.05, // +- 0.02
   chi32 = -0.03, // +- 0.03
   eta0 = 0.74, // +- 0.11
   eta1 = 0.05, // +- 0.03
   eta2 = -0.05, // +- 0.05
   lh10 = 0.09, // +- 0.18
   lh11 = 1.20, // +- 2.09
   lh20 = -2.29, // +- 0.33
   lh21 = -3.66, // +- 1.56
   lh30 = -1.90, // +- 12.4
   lh31 = 3.91, // +- 4.35
   lh40 = -2.56, // +- 0.94
   lh41 = 1.78, // +- 0.93
   lh50 = 3.96, // +- 1.17
   lh51 = 2.10, // +- 1.47
   lh60 = 4.96, // +- 5.76
   lh61 = 5.08; // +- 2.97 */
   
 // HQET FF parameters based on Table 2 Column 4 (SM 2/1/0) arxiv:2004.10208
/*   const double 
   Xi0 = 1,
   Xi1 = -1.10, // +- 0.04
   Xi2 = 1.57, // +- 0.10
   Xi3 = 0,
   chi20 = -0.06, // +- 0.02
   chi21 = 0.01,  // +- 0.02
   chi22 = 0,
   chi30 = 0,
   chi31 = -0.03, // +- 0.01
   chi32 = 0,
   eta0 = 0.38, // +- 0.06
   eta1 = 0.08, // +- 0.03
   eta2 = 0,
   lh10 = 0.5, // +- 0.16
   lh11 = 0,
   lh20 = -2.16, // +- 0.29
   lh21 = 0,
   lh30 = -1.14, // +- 2.34
   lh31 = 0,
   lh40 =0.82, // +- 0.47
   lh41 = 0,
   lh50 = 1.39, // +- 0.43
   lh51 = 0,
   lh60 = 0.17, // +- 1.15
   lh61 = 0;
 
   //Expansion coefficients for alpha_s, 1/m and 1/m^2 corrections (2004.10208 , below eq 10)
  const double epsa = 0.0716, epsb = 0.0522, epsc = 0.1807;
  
  //Some factors required for higher order corrections of the FFs
  double zcb = m_c/m_b, zcbsq = zcb*zcb;
  double wcb = 0.5*(zcb + 1/zcb);
  double omep = w + sqrt((w+1)*(w-1));
  double omem = w - sqrt((w+1)*(w-1));
  double rw = log(omep)/sqrt((w+1)*(w-1));
  double epscsq = epsc*epsc;
  double Omew = w*(2*EvtDiLog::DiLog(1-omem*zcb) - 2*EvtDiLog::DiLog(1-omep*zcb) + EvtDiLog::DiLog(1-omep*omep) - EvtDiLog::DiLog(1-omem*omem))/(2*sqrt(wsq-1)) - w*rw*log(zcb) + 1;
  
  //double scalemusq = m_b*m_c;
    const double scalemusq = 4.2*4.2;
  
  // The order alpha_s/pi corrections arxiv:1703.05330 Appendix A1 with scale factors taken from arxiv:2004.10208 Eq. 87-89
  double SFCV1 = -2*(w*rw -1)*log(m_b*m_c/scalemusq)/3;
  double SFCA1 = -2*(w*rw -1)*log(m_b*m_c/scalemusq)/3;
  //double SFCP = -(2*w*rw + 1)*log(m_b*m_c/scalemusq)/3;
  double SFCT1 = -(2*w*rw - 3)*log(m_b*m_c/scalemusq)/3;
  
  
  double CV1w = (2*(w+1)*((3*w-1)*zcb-zcbsq-1)*rw + (12*zcb*(wcb-w) - (zcbsq-1)*log(zcb)) + 4*zcb*(w-wcb)*Omew)/(6*zcb*(w-wcb)) + SFCV1;
  double CV2w = -(((4*wsq + 2*w)*zcbsq - (2*wsq + 5*w -1)*zcb - (w+1)*zcb*zcbsq + 2)*rw + zcb*(2*(zcb-1)*(wcb-w) + (zcbsq - (4*w -2)*zcb + (3-2*w))*log(zcb)))/(6*zcbsq*(w-wcb)*(w-wcb));
  double CV3w = (((2*wsq + 5*w -1)*zcbsq - (4*wsq +2*w)*zcb - 2*zcb*zcbsq + w + 1)*rw + (2*zcb*(zcb-1)*(wcb-w)+((3-2*w)*zcbsq + (2-4*w)*zcb +1)*log(zcb)))/(6*zcb*(w-wcb)*(w-wcb));
  
  double CA1w = (2*(w-1)*((3*w+1)*zcb-zcbsq-1)*rw + (12*zcb*(wcb-w)-(zcbsq-1)*log(zcb)) +4*zcb*(w-wcb)*Omew)/(6*zcb*(w-wcb)) + SFCA1;
  double CA2w = -(((4*wsq-2*w)*zcbsq + (2*wsq-5*w-1)*zcb +(1-w)*zcb*zcbsq +2)*rw + zcb*(2*(zcb+1)*(wcb-w) + (zcbsq - (4*w + 2)*zcb + (3+2*w))*log(zcb)))/(6*zcbsq*(w-wcb)*(w-wcb));
  double CA3w = ((2*zcb*zcbsq + (2*wsq-5*w-1)*zcbsq + (4*wsq-2*w)*zcb -w +1)*rw + (2*zcb*(zcb+1)*(wcb-w)-((3+2*w)*zcbsq -(2+4*w)*zcb +1)*log(zcb)))/(6*zcb*(w-wcb)*(w-wcb));
  
  //double CPw = (2*zcb*(w-wcb)*Omew - (w+1)*(zcb-1)*(zcb-1)*rw +(zcbsq -1)*log(zcb)) + SFCP;
  
  double CT1w = ((w-1)*((4*w+2)*zcb -zcbsq -1)*rw + (6*zcb*(wcb-w)-(zcbsq-1)*log(zcb))+2*zcb*(w-wcb)*Omew)/(3*zcb*(w-wcb)) + SFCT1;
  double CT2w = 2*((1-w*zcb)*rw + zcb*log(zcb))/(3*zcb*(w-wcb));
  double CT3w = 2*((w-zcb)*rw + log(zcb))/(3*(w-wcb));
  
  double dhaV = CV1w;
  double dhaA1 = CA1w;
  double dhaA2 = CA2w;
  double dhaA3 = CA1w + CA3w;
  //double dhaP = CPw;
  double dhaT1 = CT1w + 0.5*(w-1)*(CT2w-CT3w);
  double dhaT2 = 0.5*(w+1)*(CT2w+CT3w);
  double dhaT3 = -0.5*CT2w;
  
 //Definition of Xi-function
  double Xiw = Xi0 + 8*Xi1*z + 16*(Xi1+2*Xi2)*pow(z,2) + 8*(9*Xi1 + 48*Xi2 + 32*Xi3)*pow(z,3);
  
  //1/mb, 1/mc corrections
  double etaw = eta0 + 8*eta1*z + 16*(eta1 + 2*eta2)*pow(z,2);
  double chi2w = chi20 + 8*chi21*z + 16*(chi21 + 2*chi22)*pow(z,2);
  double chi3w = chi30 + 8*chi31*z + 16*(chi31 + 2*chi32)*pow(z,2);
  
  // L-hat functions arxiv:1703.05330 Eq 8 for 1/m corrections.
  double L1hat = -4*(w-1)*chi2w + 12*chi3w;
  double L2hat = -4*chi3w;
  double L3hat = 4*chi2w;
  double L4hat = 2*etaw-1;
  const double L5hat = -1;
  double L6hat = -2*(1+etaw)/(w+1);
  
  // 1/m corrections from arxiv:1703.05330 Eq.15
  double dhbV = L1hat - L4hat;
  double dhcV = L2hat - L5hat;
  
  double dhbA1 = L1hat - L4hat*(w-1)/(w+1);
  double dhcA1 = L2hat - L5hat*(w-1)/(w+1);
  
  const double dhbA2 = 0;
  double dhcA2 = L3hat + L6hat;
  
  double dhbA3 = L1hat - L4hat;
  double dhcA3 = L2hat - L3hat + L6hat - L5hat;
  
  //double dhbP = L1hat - L4hat;
  //double dhcP = L2hat + L3hat*(w-1) + L5hat - L6hat*(w+1);
  
  double dhbT1 = L1hat;
  double dhcT1 = L2hat;
  
  double dhbT2 = -L4hat;
  double dhcT2 = L5hat;
  
  const double dhbT3 = 0;
  double dhcT3 = -0.5*(L6hat - L3hat);
  
  // lhat functions for 1/mc^2 corrections.
  double lh1 = lh10 + 8*lh11*z;
  double lh2 = lh20 + 8*lh21*z;
  double lh3 = lh30 + 8*lh31*z;
  double lh4 = lh40 + 8*lh41*z;
  double lh5 = lh50 + 8*lh51*z;
  double lh6 = lh60 + 8*lh61*z;
  
  // 1/mc^2 corrections from arxiv:2004.10208 Eq.103-114
  // The small lh functions are related to Lhat functions as mention in Eq.11 of arxiv:2004.10208
  double dhccV = lh2 - lh5;
  double dhccA1 = lh2 - lh5*(w-1)/(w+1);
  double dhccA2 = lh3 + lh6;
  double dhccA3 = lh2 - lh3 - lh5 + lh6;
  //double dhccP = lh2 + (w-1)*lh3 + lh5 - lh6*(w+1);
  double dhccT1 = lh2;
  double dhccT2 = lh5;
  double dhccT3 = 0.5*(lh3 - lh6); 
    
  
  // Upto 1/mc^2 corrections included as defined in arxiv:2004.10208 Eq.9 & 10
  double hVhat = 1 + epsa*dhaV + epsb*dhbV + epsc*dhcV + epscsq*dhccV;
  double hA1hat = 1 + epsa*dhaA1 + epsb*dhbA1 + epsc*dhcA1 + epscsq*dhccA1;
  double hA2hat = 0 + epsa*dhaA2 + epsb*dhbA2 + epsc*dhcA2 + epscsq*dhccA2;
  double hA3hat = 1 + epsa*dhaA3 + epsb*dhbA3 + epsc*dhcA3 + epscsq*dhccA3;
  //double hPhat = 1 + epsa*dhaP + epsb*dhbP + epsc*dhcP + epscsq*dhccP;
  double hT1hat = 1 + epsa*dhaT1 + epsb*dhbT1 + epsc*dhcT1 + epscsq*dhccT1;
  double hT2hat = 0 + epsa*dhaT2 + epsb*dhbT2 + epsc*dhcT2 + epscsq*dhccT2;
  double hT3hat = 0 + epsa*dhaT3 + epsb*dhbT3 + epsc*dhcT3 + epscsq*dhccT3;
  
  // Higher-order corrected HQET FFs
  double hVw = Xiw*hVhat;
  double hA1w = Xiw*hA1hat;
  double hA2w = Xiw*hA2hat;
  double hA3w = Xiw*hA3hat;
  double hT1 = Xiw*hT1hat;
  double hT2 = Xiw*hT2hat;
  double hT3 = Xiw*hT3hat; */
  
// HQET Parametrization ends here.
/************************************************************************************/ 
 
/************************************************************************************/
  // BGL parametrization of the FF begins here :
  //FF parameter values taken from 2111.01176
  const double 
   a0f = 0.0123, // +- 0.0001
   a1f = 0.0222, // +- 0.0096
   a2f = -0.522, // +- 0.196
   a0g = 0.0318, // +- 0.0010
   a1g = -0.133, // +- 0.063
   a2g = -0.62,  // +- 1.46
   a1F1 = 0.0021, // +- 0.0015 
   a0F2 = 0.0515, // +- 0.0021
   a1F2 = -0.149, // +- 0.059
   a2F2 = 0.987, // +- 0.932
   nI = 2.6,
   chi1minT0 = 5.131e-4,
   chi1plusT0 = 3.894e-4,
   chi1plusL0 = 1.9421e-2; 
  
   double mP1plus[] = {6.739,6.750,7.145,7.150};
   double mP1min[] = {6.329,6.920,7.020};
   double mP0min[] = {6.275,6.842,7.250};
   
   double zP1plus[4], zP1min[3], zP0min[3];
   
   double P1plus = 1.0;
   double P1min = 1.0;
   double P0min = 1.0;
   
   for(int i=0; i<4; i++){
        zP1plus[i] = (sqrt(mBaV*mBaV - mP1plus[i]*mP1plus[i]) - sqrt(4*mB*mDst))/(sqrt(mBaV*mBaV - mP1plus[i]*mP1plus[i]) + sqrt(4*mB*mDst));
   	P1plus*= (z-zP1plus[i])/(1-z*zP1plus[i])
   }
   	
   for(int i=0; i<3; i++){
        zP1min[i] = (sqrt(mBaV*mBaV - mP1min[i]*mP1min[i]) - sqrt(4*mB*mDst))/(sqrt(mBaV*mBaV - mP1min[i]*mP1min[i]) + sqrt(4*mB*mDst));
        zP0min[i] = (sqrt(mBaV*mBaV - mP0min[i]*mP0min[i]) - sqrt(4*mB*mDst))/(sqrt(mBaV*mBaV - mP0min[i]*mP0min[i]) + sqrt(4*mB*mDst));
   	P1min*= (z-zP1min[i])/(1-z*zP1min[i]);
   	P0min*= (z-zP0min[i])/(1-z*zP0min[i]);
   }
   
   double phig = 16*r2Dstar*sqrt(nI/(3*M_PI*chi1minT0))*((pow(1+z,2))/(sqrt(1-z)*pow((1+rDstar)*(1-z)+2*sqrt(rDstar)*(1+z),4)));
   double phif = 4*rDstar*sqrt(nI/(3*M_PI*chi1plusT0))*(((1+z)*(1-z)*sqrt(1-z))/(pow((1+rDstar)*(1-z)+2*sqrt(rDstar)*(1+z),4)))/(mB*mB);
   double phiF1 = 4*rDstar*sqrt(nI/(6*M_PI*chi1plusT0))*(((1+z)*pow(1-z,2)*sqrt(1-z))/(pow((1+rDstar)*(1-z)+2*sqrt(rDstar)*(1+z),5)))/(pow(mB,3));
   double phiF2 = 8*sqrt(2)*r2Dstar*sqrt(nI/(M_PI*chi1plusL0))*((pow(1+z,2))/(sqrt(1-z)*pow((1+rDstar)*(1-z)+2*sqrt(rDstar)*(1+z),4)));
   
   //Kinetic constraints :
   const double a0F1 = (1-rDstar)/(sqrt(2)*pow(1+sqrt(rDstar),2)); // at w=1
   const double a2F1 = -0.5162*(104.455*a0f+34.7622*a1F1-29.3161*(a0F2+0.0557*a1F2+0.0031*a2F2)); // obtained at w = wmax = 1.5 after substituting all other constants
   
   double gfunc = (a0g + a1g*z + a2g*z*z)/(P1min*phig);
   double ffunc = (a0f + a1f*z + a2f*z*z)/(P1plus*phif);
   double F1func = (a0F1 + a1F1*z + a2F1*z*z)/(P1plus*phiF1);
   double F2func = (a0F2 + a1F2*z + a2F2*z*z)/(P0min*phiF2);
   
   double hVw = mB*sqrt(rDstar)*gfunc;
   double hA1w = ffunc/(mB*sqrt(rDstar)*(1+w));
   double hA2w = hA1w/(1-w) + (F1func*(w-rDstar) - F2func*mB*mB*rDstar*(wsq-1))/(mB*mB*sqrt(rDstar)*u*(wsq-1));
   double hA3w = (F1func*(rDstar*w-1) + mB*mB*sqrt(rDstar)*(1+w)*(F2func*rDstar*sqrt(rDstar)*(w-1) + hA1w*w*u))/(mB*mB*sqrt(rDstar)*u*(wsq-1));
   
  // arXiv:1309.0301 Eq 48b
  double cA1 = (m_b-m_c)*rmBsV*hA1w, cV = (m_b+m_c)*rmBaV*hVw;
  double hT1 = (0.5*((1-rDstar)*(1-rDstar)*(w+1)*cA1
		  -(1+rDstar)*(1+rDstar)*(w-1)*cV))/u;
  double hT2 = (0.5*(1-r2Dstar)*(w+1)*(cA1 - cV))/u;
  double hT3 = (-0.5*(2*rDstar*(w+1)*cA1
				+(m_b-m_c)*rmBsV*u*(hA3w-rDstar*hA2w)
				-(1+rDstar)*(1+rDstar)*cV))/(u*(1+rDstar));
 
// BGL Parametrization ends here.
/************************************************************************************/ 

// Form Factors :  
// arXiv:1309.0301 Eq 39
  T1 = mBaV*hT1 - mBsV*hT2;
  T2 = (mBaV-q2*rmBaV)*hT1 - (mBsV-q2*rmBsV)*hT2;
  T3 = mBsV*hT1 - mBaV*hT2 - 2*mBaV*mBsV*rmB*hT3;
  T1 *= c;
  T2 *= c;
  T3 *= c;

  V = c*mBaV*hVw;
  A1 = c*(mBaV*mBaV-q2)/mBaV*hA1w;
  A2 = c*mBaV*(hA3w+rDstar*hA2w);
  A0 = c*((mBaV*mBaV-q2)/(2*mV)*hA1w - (mBaV*mBsV+q2)/(2*mB)*hA2w - (mBaV*mBsV-q2)/(2*mV)*hA3w);
  
}

static EvtId EM, MUM, TAUM, EP, MUP, TAUP, D0, D0B, DP, DM, DSM, DSP;
static bool cafirst = true;
using EvtGenFunctions::asymProd;
void EvtSemiLeptonicVectorAmpNP::CalcAmp(EvtParticle* parent, EvtAmp& amp, EvtSemiLeptonicFF*){
  //  return  SMCalcAmpVerb(parent, amp, FormFactors);
  if(cafirst){
    cafirst = false;
    EM = EvtPDL::getId( "e-" ); EP = EvtPDL::getId( "e+" );
    MUM = EvtPDL::getId( "mu-" ); MUP = EvtPDL::getId( "mu+" );
    TAUM = EvtPDL::getId( "tau-" ); TAUP = EvtPDL::getId( "tau+" );
    
    D0 = EvtPDL::getId( "D0" ); D0B = EvtPDL::getId( "anti-D0" );
    DP = EvtPDL::getId( "D+" ); DM = EvtPDL::getId( "D-" );
    DSP = EvtPDL::getId( "D_s+" ); DSM = EvtPDL::getId( "D_s-" );
  }
  
  EvtParticle
    *pvm = parent->getDaug( 0 ), // vector meson
    *plp = parent->getDaug( 1 ), // charged lepton
    *pnu = parent->getDaug( 2 ); // neutrino
  EvtId pId = parent->getId(), vId = pvm->getId(), lId = plp->getId();
  
  //Add the lepton and neutrino 4 momenta to find q2
  EvtVector4R q = plp->getP4() + pnu->getP4();
  double q2 = q.mass2(), m_V = pvm->mass(), m_B = parent->mass();
  double m_BaV = m_B + m_V, m_BsV = m_B - m_V;
  
  double a1f, a2f, vf, a0f, a3f;
  /*  Tensor form factors */
  double T1, T2, T3;
  getFF(m_B, m_V, q2, vf, a0f, a1f, a2f, T1, T2, T3);
  a3f = ( m_BaV * a1f - m_BsV * a2f ) / ( 2.0 * m_V );
 
  double c0 = T1;
  double c1 = (m_BaV*m_BsV*(T1-T2))/q2;
  double c2 = (T1 - T2 - T3*q2/(m_BaV*m_BsV))*2/q2;

  if ( pId == D0 || pId == D0B || pId == DP || pId == DM || pId == DSP || pId == DSM)
    vf = -vf;

  EvtVector4R p( m_B, 0., 0., 0. );
  const EvtVector4R &k = pvm->getP4();
  EvtVector4R pak = p + k;
    
  EvtDiracSpinor
    nu = pnu->spParentNeutrino(), lp0 = plp->spParent(0), lp1 = plp->spParent(1);

  /* Various lepton currents */
  EvtComplex s1, s2, p1, p2;
  EvtVector4C l1, l2, a1, a2;
  EvtTensor4C z1, z2, w1, w2;
  if ( lId == EM || lId == MUM || lId == TAUM ) {
    l1 = EvtLeptonVCurrent( nu, lp0 );
    l2 = EvtLeptonVCurrent( nu, lp1 );
    a1 = EvtLeptonACurrent( nu, lp0 );
    a2 = EvtLeptonACurrent( nu, lp1 );
    s1 = EvtLeptonSCurrent( nu, lp0 );
    s2 = EvtLeptonSCurrent( nu, lp1 );
    p1 = EvtLeptonPCurrent( nu, lp0 );
    p2 = EvtLeptonPCurrent( nu, lp1 );
    z1 = EvtLeptonTCurrent( nu, lp0 );
    z2 = EvtLeptonTCurrent( nu, lp1 );
  } else if ( lId == EP || lId == MUP || lId == TAUP ) {
    vf = -vf; /*This takes care of sign flip in vector hadronic current when you conjugate it*/
    l1 = EvtLeptonVCurrent( lp0, nu );
    l2 = EvtLeptonVCurrent( lp1, nu );
    a1 = EvtLeptonACurrent( lp0, nu );
    a2 = EvtLeptonACurrent( lp1, nu );
    s1 = EvtLeptonSCurrent( lp0, nu );
    s2 = EvtLeptonSCurrent( lp1, nu );
    p1 = EvtLeptonPCurrent( lp0, nu );
    p2 = EvtLeptonPCurrent( lp1, nu );
    z1 = EvtLeptonTCurrent( lp0, nu );
    z2 = EvtLeptonTCurrent( lp1, nu );
  }

  /* Quark masses from PDG */
  double m_c = 1.27, m_b = 4.18;
  
  /* New Physics */
  EvtComplex C_V_L = 1, //Standard Model
    C_V_R = _Cvr,
    C_S_L = _Csl, 
    C_S_R = _Csr, 
    C_T = _cT;
  C_V_L += _Cvl;
 
  /*Eq 12 Axial vector hadronic tensor*/
  EvtTensor4C tA = (a1f * m_BaV) * EvtTensor4C::g(); 
  tA.addDirProd(-(a2f/m_BaV) * p, pak);
  tA.addDirProd( -2 * m_V/q2 * (a3f - a0f) * p, q); 


  /*Eq 11 Vector hadronic tensor*/
  EvtTensor4C tV = EvtComplex(0.0, (-2*vf/m_BaV)) * asymProd(p, k);

  EvtTensor4C tVA = 0.5*(C_V_L + C_V_R)*tV + 0.5*(C_V_R - C_V_L)*tA;
  
  EvtComplex tS;
  /*Eq 13 Pseudoscalar hadronic vector*/
  /*Scalar hadronic vector is zero*/
  /*  if ( lId == EM || lId == MUM || lId == TAUM ) {
    tS = -0.5*(C_S_R - C_S_L) * 2*m_V/(m_b + m_c)*a0f;
  } else if ( lId == EP || lId == MUP || lId == TAUP ) {
    tS = 0.5*(C_S_R - C_S_L) * 2*m_V/(m_b + m_c)*a0f;
  }*/
  tS = 0.5*(C_S_R - C_S_L) * 2*m_V/(m_b + m_c)*a0f;
  
  //z1 *= C_T;
  //z2 *= C_T;

  pak *= c0;
  EvtVector4C qc1(q); qc1 *= c1;
  EvtTensor4C tt,tt5;
  EvtComplex sum1, sum2, TT1, TT2;
  for(int i=0;i<3;i++){
    EvtVector4C eps = pvm->epsParent(i).conj(), epsTVA = tVA.cont1( eps );
    EvtComplex epsq = eps*q, epsqc2 = epsq*c2;

    /*Eq 14 Tensor component */
    for(int i2=0;i2<4;i2++){
      for(int i3=0;i3<4;i3++){
      	tt.set(i2, i3, eps.get(i2)*(-pak.get(i3) + qc1.get(i3)) + p.get(i2)*k.get(i3)*epsqc2);	
      	//tt5.set(i2, i3, EvtComplex(0.0,1.0)*(-eps.get(i2)*(pak.get(i3)) + pak.get(i2)*eps.get(i3) + eps.get(i2)*qc1.get(i3) -qc1.get(i2)*eps.get(i3) + epsqc2*(p.get(i2)*k.get(i3)-k.get(i2)*p.get(i3))));
      }
    }
    if ( lId == EM || lId == MUM || lId == TAUM ) {
    	TT1 = EvtComplex(0.0,1.0)*C_T*(asymProd(z1,tt));
    	TT2 = EvtComplex(0.0,1.0)*C_T*(asymProd(z2,tt));
    	sum1 = l1.cont(epsTVA) - a1.cont(epsTVA) - epsq*tS*p1 + TT1; 
    	sum2 = l2.cont(epsTVA) - a2.cont(epsTVA) - epsq*tS*p2 + TT2;
  } else if ( lId == EP || lId == MUP || lId == TAUP ) {
  	TT1 = EvtComplex(0.0,-1.0)*conj(C_T)*(asymProd(z1,tt));
  	TT2 = EvtComplex(0.0,-1.0)*conj(C_T)*(asymProd(z2,tt));
    	sum1 = l1.cont(epsTVA) - a1.cont(epsTVA) - epsq*tS*p1 + TT1; 
   	sum2 = l2.cont(epsTVA) - a2.cont(epsTVA) - epsq*tS*p2 + TT2;
  }
   
    amp.vertex( i, 0, sum1);
    amp.vertex( i, 1, sum2);
  }
}

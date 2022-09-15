#include "EvtGen/EvtGen.hh"

#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtVector3R.hh"

#include "TFile.h"
#include "TNtuple.h"

#include <getopt.h>

static const char *program_name = "";
void usage(int status){
  fprintf(stderr, "Usage: %s -n nevent -b B0|anti-B0|B+|B- [-u user_decay_file] [-o out_root_file]\n", program_name);
  exit(status);
}

bool allowed_arg(const char *arg, int n, const char* args[] ){
  bool match = false;
  for(int i=0;i<n;i++) if(strcmp(arg,args[i])==0) match = true;
  return match;
}

int main( int argc, char* argv[] ){
  program_name = argv[0];
  if(argc<2) usage(EXIT_FAILURE);
  const char *mnames[] = {"B0","anti-B0","B+","B-"};
  const char *lnames[] = {"e","mu","tau"};
  const char *mname = "B0", *lname, *fname=NULL, *ufile = "bb_dstarlnu_np.dec";
  char *endptr;
  long int nevent = 1000;
  int opt, seed = 1;
  while ((opt = getopt(argc, argv, "n:b:o:u:s:")) != -1) {
    switch (opt) {
    case 'n':
      nevent = strtoul(optarg, &endptr, 0);
      if(endptr!=NULL && *endptr!='\0') usage(EXIT_FAILURE);
      break;
    case 'b':
      if(!allowed_arg(optarg,4,mnames)) usage(EXIT_FAILURE);
      mname = optarg;
      break;
    case 'o':
      fname = optarg;
      break;
    case 'u':
      ufile = optarg;
      break;
    case 's':
      seed = atoi(optarg);
      break;
    default: /* '?' */
      usage(EXIT_FAILURE);
    }
  }

  EvtRandomEngine *RandEng = new EvtMTRandomEngine(seed);
  EvtGen Generator( "../myDECAY.DEC", "../evt.pdl", RandEng, NULL, NULL );
  Generator.readUDecay(ufile);
  EvtId Mother = EvtPDL::getId("Upsilon(4S)");
  double pe = 7, pp = 4, thetac = 0.083/2, px = sin(thetac)*(pe+pp), pz = cos(thetac)*(pe-pz), M = EvtPDL::getMass(Mother);
  EvtVector4R pY( sqrt(M*M + px*px + pz*pz), px, 0.0, pz );

  if(fname){
    TFile* file = new TFile( fname, "RECREATE" );
    TNtuple *ntp = new TNtuple("ntp","","mD:q2:ctd:ctl:chi:lct:lp:nuct:nup:Dct:Dp:pict:pip");
    
    long int count = 1;
    do {
      EvtParticle* mY = EvtParticleFactory::particleFactory( Mother, pY );
      mY->setDiagonalSpinDensity();
      Generator.generateDecay( mY );
      
      EvtParticle* mp = mY->getDaug(0);

      EvtVector4R dstar = mp->getDaug( 0 )->getP4Lab();
      EvtVector4R l1 = mp->getDaug( 1 )->getP4Lab();
      EvtVector4R l2 = mp->getDaug( 2 )->getP4Lab();
      EvtVector4R q = l1 + l2;
      EvtVector4R b = mp->getP4Lab();
      EvtVector4R k = mp->getDaug( 0 )->getDaug( 0 )->getP4Lab();
      EvtVector4R pi = mp->getDaug( 0 )->getDaug( 1 )->getP4Lab();
      
      // Kinematic of the decay B -> D*(D pi) l nu is fully described by 4 parameters:
      double q2 = q.mass2(); // q^2 -- hadronic recoil
      double ctk = EvtDecayAngle(b, dstar, k); // cos(theta_D)
      double ctl = EvtDecayAngle(b, q, l1); // cos(theta_ell)
      double chi = EvtDecayAngleChi(b, k, pi, l1, l2); // -pi<chi<pi angle
      
      double mDpi = (k+pi).mass(); // for convenience -- invariant mass of the D and pion system or D*
      double lpx = l1.get(1), lpy = l1.get(2), lpz = l1.get(3), lp = sqrt(lpx*lpx+lpy*lpy+lpz*lpz);
      double pipx = pi.get(1), pipy = pi.get(2), pipz = pi.get(3), pip = sqrt(pipx*pipx+pipy*pipy+pipz*pipz);
      double Dpx = k.get(1), Dpy = k.get(2), Dpz = k.get(3), Dp = sqrt(Dpx*Dpx+Dpy*Dpy+Dpz*Dpz);
      double nupx = l2.get(1), nupy = l2.get(2), nupz = l2.get(3), nup = sqrt(nupx*nupx+nupy*nupy+nupz*nupz);
      ntp->Fill(mDpi, q2, ctk, ctl, chi, lpz/lp, lp, nupz/nup, nup, Dpz/Dp, Dp, pipz/pip, pip);
      mY->deleteTree();
    } while ( count++ < nevent );
    
    file->Write();
    file->Close();
  } else {
    long int count = 1;
    do {
      EvtParticle* mp = EvtParticleFactory::particleFactory( Mother, pY );
      mp->setDiagonalSpinDensity();
      Generator.generateDecay( mp );
      mp->deleteTree();
    } while ( count++ < nevent );
  }
  delete RandEng;
  return 0;
}

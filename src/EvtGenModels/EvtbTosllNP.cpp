
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

#include "EvtGenModels/EvtbTosllNP.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenModels/EvtbTosllAmp.hh"
#include "EvtGenModels/EvtbTosllBSZFF.hh"
#include "EvtGenModels/EvtbTosllVectorAmpNP.hh"

using std::endl;

std::string EvtbTosllNP::getName()
{
  return "BTOSLLNP";
}

EvtDecayBase* EvtbTosllNP::clone()
{
  return new EvtbTosllNP;
}

void EvtbTosllNP::decay( EvtParticle* p )
{
  setWeight( p->initializePhaseSpace( getNDaug(), getDaugs(), false, _poleSize, 1, 2 ) );
  _calcamp->CalcAmp( p, _amp2, _ffmodel.get() );
}

void EvtbTosllNP::initProbMax()
{
  EvtId pId = getParentId(), mId = getDaug( 0 ), l1Id = getDaug( 1 ), l2Id = getDaug( 2 );
  //This routine sets the _poleSize.
  double mymaxprob = _calcamp->CalcMaxProb( pId, mId, l1Id, l2Id, _ffmodel.get(), _poleSize );
  setProbMax( mymaxprob );
}

void EvtbTosllNP::init()
{
  // First choose format of NP coefficients from the .DEC file
  // Cartesian(x,y)(0) or Polar(R,phase)(1)
  int n = getNArg();
  if(!(n==0||(n-1)%3==0)){
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
      << "Error in parameters in the BTOSLLNP decay model." << endl;
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
      << "Will terminate execution!" << endl;
    ::abort();
  }

  checkNDaug( 3 );

  //We expect the parent to be a scalar
  //and the daughters to be K* lepton+ lepton-

  checkSpinParent( EvtSpinType::SCALAR );
  checkSpinDaughter( 1, EvtSpinType::DIRAC );
  checkSpinDaughter( 2, EvtSpinType::DIRAC );

  EvtId mId = getDaug( 0 );
  EvtId IdKst0 = EvtPDL::getId("K*0"), IdaKst0 = EvtPDL::getId("anti-K*0"),
    IdKstp = EvtPDL::getId("K*+"), IdKstm = EvtPDL::getId("K*-");
  if ( mId != IdKst0 && mId != IdaKst0 && mId != IdKstp && mId != IdKstm ) {
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
      << "EvtbTosllNP generator expected a K* 1st daughter, found: "
      << EvtPDL::name( getDaug( 0 ) )  << endl;
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
      << "Will terminate execution!" << endl;
    ::abort();
  }

  _ffmodel = std::make_unique<EvtbTosllBSZFF>();
  _calcamp = std::make_unique<EvtbTosllVectorAmpNP>();

  auto getInteger = [this](int i) -> int {
    double a = getArg(i);
    if(a-static_cast<int>(a) != 0){
      EvtGenReport( EVTGEN_ERROR, "EvtGen" )
	<< "Parameters is not integer in the BTOSLLNP decay model: " << i<<" "<< a << endl;
      EvtGenReport( EVTGEN_ERROR, "EvtGen" )
	<< "Will terminate execution!" << endl;
      ::abort();
    }
    return static_cast<int>(a);
  };
  EvtbTosllVectorAmpNP *amp =  static_cast<EvtbTosllVectorAmpNP*>(_calcamp.get());
  if(n>0) { // parse arguments
    int i = 0, coordsyst = getInteger( i++ );
    while(i<n){
      int id =  getInteger( i++ ); // New Physics cooeficient Id
      double a0 = getArg(i++);
      double a1 = getArg(i++);
      EvtComplex c = (coordsyst)? EvtComplex(a0*cos(a1), a0*sin(a1)):EvtComplex(a0, a1);
      if(id==0) amp->_c7p = c; // C'_7eff -- right hand polarizations
      if(id==1) amp->_c9p = c; // C'_9eff -- right hand polarizations
      if(id==2) amp->_c10p = c; // c'_10eff -- right hand polarizations
      if(id==3) amp->_cS = c; // (C_S - C'_S) -- scalar right and left polarizations
      if(id==4) amp->_cP = c; // (C_P - C'_P) -- pseudo-scalar right and left polarizations
    }
  }

}

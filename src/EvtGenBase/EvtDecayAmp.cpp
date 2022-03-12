
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

#include "EvtGenBase/EvtDecayAmp.hh"

#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;

void EvtDecayAmp::makeDecay( EvtParticle* p, bool recursive ){
  //original default value
  int ntimes = 10000;
  int more;
  EvtSpinDensity rho;
  double prob, prob_max;
  _amp2.init( p->getId(), getNDaug(), getDaugs() );
  do {
    _daugsDecayedByParentModel = false;
    _weight = 1.0;
    decay( p );
    rho = _amp2.getSpinDensity();
    //std::cout<<"phase space loop "<<p->getName()<<p->getSpinDensityForward()<<" "<<rho<<std::endl;
    prob = p->getSpinDensityForward().normalizedProb( rho );
    //std::cout<<prob<<std::endl;
    prob /= _weight;
    //	std::cout<<_weight<<std::endl;
    prob_max = getProbMax( prob );
    p->setDecayProb( prob / prob_max );
    more = prob < EvtRandom::Flat( prob_max );
    ntimes--;
  } while ( ntimes && more );

  // it may be that the parent decay model has already
  // done the decay - this should be rare and the
  // model better know what it is doing..

  if ( !daugsDecayedByParentModel() ) {
    if ( recursive ) {
      EvtSpinDensity rho_list[10];
      rho_list[0] = p->getSpinDensityForward();
      EvtAmp ampcont;
      if ( _amp2._pstates != 1 ) {
	ampcont = _amp2.contract( 0, p->getSpinDensityForward() );
      } else {
	ampcont = _amp2;
      }
      for ( size_t i = 0; i < p->getNDaug(); i++ ) {
	//std::cout<<"loop "<<p->getName()<<" "<<i<<std::endl;
	rho.setDim( _amp2.dstates[i] );
	if ( _amp2.dstates[i] == 1 ) {
	  rho.set( 0, 0, EvtComplex( 1.0, 0.0 ) );
	} else {
	  //std::cout<<"cont "<<p->getName()<<" "<<i<<" "<< _amp2._dnontrivial[i]<<std::endl;
	  rho = ampcont.contract( _amp2._dnontrivial[i], _amp2 );
	}
	// std::cout<<rho;
	p->getDaug( i )->setSpinDensityForward( rho );
	p->getDaug( i )->decay();
	rho_list[i + 1] = p->getDaug( i )->getSpinDensityBackward();
	if ( _amp2.dstates[i] != 1 ) {
	  ampcont = ampcont.contract( _amp2._dnontrivial[i], rho_list[i + 1] );
	}
      }
      p->setSpinDensityBackward( _amp2.getBackwardSpinDensity( rho_list ) );
    }
  }

  if ( getPHOTOS() || EvtRadCorr::alwaysRadCorr() ) {
    int n_daug_orig = p->getNDaug();
    EvtRadCorr::doRadCorr( p );
    int n_daug_new = p->getNDaug();
    for ( int i = n_daug_orig; i < n_daug_new; i++ ) {
      p->getDaug( i )->decay();
    }
  }
}

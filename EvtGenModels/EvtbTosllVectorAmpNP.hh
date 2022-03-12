
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

#ifndef EVTBTOSLLVECTORAMPNP_HH
#define EVTBTOSLLVECTORAMPNP_HH

#include "EvtGenModels/EvtbTosllAmp.hh"
#include "EvtGenBase/EvtComplex.hh"

class EvtAmp;
class EvtParticle;
class EvtbTosllFF;

class EvtbTosllVectorAmpNP : public EvtbTosllAmp {
public:
  EvtbTosllVectorAmpNP():_c7p(0),_c9p(0),_c10p(0),_cS(0),_cP(0){}
  void CalcAmp( EvtParticle* parent, EvtAmp& amp, EvtbTosllFF* formFactors ) override;

  EvtComplex _c7p; // C'_7eff -- right hand polarizations
  EvtComplex _c9p; // C'_9eff -- right hand polarizations
  EvtComplex _c10p; // c'_10eff -- right hand polarizations
  EvtComplex _cS; // (C_S - C'_S) -- scalar right and left polarizations
  EvtComplex _cP; // (C_P - C'_P) -- pseudo-scalar right and left polarizations
};

#endif


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

#ifndef EVTGEN_HH
#define EVTGEN_HH

#include "EvtGenBase/EvtPDL.hh"

#include <list>

class EvtParticle;
class EvtRandomEngine;
class EvtVector4R;
class EvtStdHep;
class EvtSpinDensity;
class EvtAbsRadCorr;
class EvtDecayBase;
class EvtHepMCEvent;

class EvtGen {
  public:
    EvtGen( const char* const decayName, const char* const pdtTableName,
            EvtRandomEngine* randomEngine = 0, EvtAbsRadCorr* isrEngine = 0,
            const std::list<EvtDecayBase*>* extraModels = 0, int mixingType = 1,
            bool useXml = false );

    ~EvtGen();

    void readUDecay( const char* const udecay_name, bool useXml = false );

    EvtHepMCEvent* generateDecay( int PDGid, EvtVector4R refFrameP4,
                                  EvtVector4R translation,
                                  EvtSpinDensity* spinDensity = 0 );

    void generateDecay( EvtParticle* p );

  private:
    EvtPDL _pdl;
    int _mixingType;
};

#endif

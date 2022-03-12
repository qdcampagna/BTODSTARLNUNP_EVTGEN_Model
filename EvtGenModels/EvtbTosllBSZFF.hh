
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

#ifndef EVTBTOSLLBSZFF_HH
#define EVTBTOSLLBSZFF_HH

#include "EvtGenModels/EvtbTosllFF.hh"

class EvtId;

// Description: Form factors for b->sll according to Aoife Bharucha, David M. Straub, Roman Zwicky
// https://arxiv.org/abs/1503.05534

class EvtbTosllBSZFF : public EvtbTosllFF {
  public:
  EvtbTosllBSZFF(){}

  void getVectorFF( EvtId parent, EvtId vmeson, double q2, double vmesonmass,
		    double& a1, double& a2, double& a0, double& v, double& t1,
		    double& t2, double& t3 ) override;
};

#endif

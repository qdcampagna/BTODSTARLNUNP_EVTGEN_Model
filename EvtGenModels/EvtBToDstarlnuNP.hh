/***********************************************************************
 * Written by Alexei Sibidanov                                         *
 ***********************************************************************/

#ifndef EVTBTODSTARLNUNP_HH
#define EVTBTODSTARLNUNP_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenModels/EvtSemiLeptonicVectorAmpNP.hh"

#include <memory>

class EvtParticle;

// Description:Implementation of the B->Dstar l nu with right hand hadronic currents

class EvtBToDstarlnuNP : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void decay( EvtParticle* p ) override;
    void init() override;
    void initProbMax() override;

  private:
    std::unique_ptr<EvtSemiLeptonicVectorAmpNP> _calcamp;
    std::unique_ptr<EvtSemiLeptonicFF> _ffmodel;
};

#endif

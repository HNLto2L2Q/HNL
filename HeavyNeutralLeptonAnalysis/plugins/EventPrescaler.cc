#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class EventPrescaler : public edm::global::EDFilter<> {

  // perhaps we need better structure here (begin run etc)
public:
  explicit EventPrescaler(const edm::ParameterSet &cfg):
    prescale_{cfg.getParameter<int>("prescale")},
    offset_{cfg.getParameter<int>("offset")} {}

  ~EventPrescaler() override {}
  
  bool filter(edm::StreamID, edm::Event& evt, const edm::EventSetup&) const final;

private:
  const int prescale_;
  const int offset_;
};

bool EventPrescaler::filter(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {
  return abs(evt.eventAuxiliary().event() % prescale_) == offset_;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventPrescaler);

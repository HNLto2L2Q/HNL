// -*- C++ -*-
//
// Package:    HNL/HeavyNeutralLeptonAnalysis
// Class:      LeptonFilter
// 
/**\class LeptonFilter LeptonFilter.cc HNL/HeavyNeutralLeptonAnalysis/plugins/LeptonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Angela Taliercio
//         Created:  Mon, 25 Feb 2019 15:43:58 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

//triggers
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//
// class declaration
//

class LeptonFilter : public edm::stream::EDFilter<> {
   public:
      explicit LeptonFilter(const edm::ParameterSet&);
      ~LeptonFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------


  edm::EDGetTokenT            < pat::MuonCollection > muonsMiniAODToken_;
  edm::EDGetTokenT        < pat::ElectronCollection > electronsMiniAODToken_;
  //trigger
  edm::EDGetTokenT            < edm::TriggerResults > triggerResultsToken_;

  unsigned int    minimalNumberOfMuons_;  
  unsigned int    minimalNumberOfElectrons_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
LeptonFilter::LeptonFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  //trackProducerTag_             = iConfig.getParameter<edm::InputTag>("TrackProducerTag");
  //mu_ = consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muon_"));
  //  muonsMiniAODToken_(mayConsume<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  //  electronsMiniAODToken_(mayConsume<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),

  muonsMiniAODToken_ = consumes <pat::MuonCollection> (iConfig.getParameter<edm::InputTag>("muonSrc"));
  electronsMiniAODToken_ = consumes <pat::ElectronCollection> (iConfig.getParameter<edm::InputTag>("electronSrc"));
  triggerResultsToken_ = consumes <edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerResultSrc"));

  minimalNumberOfMuons_  = iConfig.getParameter<unsigned int>("MinimalNumberOfMuons");
  minimalNumberOfElectrons_ = iConfig.getParameter<unsigned int>("MinimalNumberOfElectrons");

}


LeptonFilter::~LeptonFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
LeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  

  edm::Handle < pat::MuonCollection > muonsHandle;
  edm::Handle < pat::ElectronCollection > electronsHandle;
  edm::Handle < edm::TriggerResults > triggers;

  int muons_count = 0, electrons_count = 0;

  iEvent.getByToken(muonsMiniAODToken_, muonsHandle);
  iEvent.getByToken(electronsMiniAODToken_, electronsHandle);
  iEvent.getByToken(triggerResultsToken_, triggers);

  //const edm::TriggerResults triggerResults =  *triggers.product();
  //const edm::TriggerNames&    trigNames  = iEvent.triggerNames(triggerResults);

  bool result = false;
  
  pat::MuonCollection muons;
  pat::ElectronCollection electrons;
  
  //muons section
  if(muonsHandle.isValid()){
    muons = *muonsHandle;
        
    for(auto& mu : muons){
      if (fabs(mu.eta()) < 2.4 && mu.pt() > 5. && mu.isLooseMuon()){
	muons_count++;
	//result = true;
      }
    }
  }

  //electrons selection
  for(auto ele = electronsHandle->begin(); ele != electronsHandle->end(); ++ele){
    if(ele->gsfTrack().isNull() || ele->pt() < 5 || fabs(ele->eta()) > 2.5 ) continue;
    if(ele->full5x5_sigmaIetaIeta() <  0.036 && ele->passConversionVeto() == 1){
      electrons_count++;
    }
  }

  
  //if((muons_count + electrons_count) >= 2){
  if((muons_count >= 2) || (electrons_count >= 2)){ 
   result = true;
  }
  
  
    

  //trigger level cut
  /*
  for (size_t i = 0; i < trigNames.size(); ++i) {
    const std::string &name = trigNames.triggerName(i);
    bool fired = triggerResults.accept(i);
    if(fired && (name.compare("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") == 0 || name.compare("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") == 0))
      result = true;
  }*/

  //number of leptons cut
  /*
    if ( (muons->size() + electrons->size() ) >= 2 && (muons->size() + electrons->size() ) <= 4) {
    result = true;
    }
  */
  return result;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
LeptonFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
LeptonFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
LeptonFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
LeptonFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
LeptonFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
LeptonFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(LeptonFilter);

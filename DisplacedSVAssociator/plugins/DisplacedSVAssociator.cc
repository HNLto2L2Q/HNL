// -*- C++ -*-
//
// Package:    TEST/DisplacedSVAssociator
// Class:      DisplacedSVAssociator
// 
/**\class DisplacedSVAssociator DisplacedSVAssociator.cc TEST/DisplacedSVAssociator/plugins/DisplacedSVAssociator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mohamed Darwish
//         Created:  Sun, 18 Jun 2017 15:14:26 GMT
//
//



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


// kinematics

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "HNL/DisplacedSVAssociator/interface/VertexAssociation.h"


//
// class declaration
//

class DisplacedSVAssociator : public edm::stream::EDProducer<> {
   public:
      explicit DisplacedSVAssociator(const edm::ParameterSet&);
      ~DisplacedSVAssociator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  typedef std::vector<float> VertexScores; 

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::VertexCollection> tag_secondaryVertices_;
  edm::EDGetTokenT<reco::VertexCollection> tag_primaryVertices_;

  std::string algoName_;
  std::string outputLabel_;
  float muonPtCut_;
  int debug_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
DisplacedSVAssociator::DisplacedSVAssociator(const edm::ParameterSet& iConfig):
  muonToken_(consumes<pat::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("muonTag"))),
  tag_secondaryVertices_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices"))),
  tag_primaryVertices_ (consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices")))

{
  algoName_ = iConfig.getUntrackedParameter<std::string>("algoName");
  muonPtCut_ = iConfig.getUntrackedParameter<double>("muonPtCut");
  debug_ = iConfig.getUntrackedParameter<int>("debug");
  outputLabel_ = iConfig.getUntrackedParameter<std::string>("outputLabel");

  produces<std::vector<float> >(outputLabel_);
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


DisplacedSVAssociator::~DisplacedSVAssociator()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DisplacedSVAssociator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //std::auto_ptr<std::vector<VertexScores > > VertexAssociationCollection(new std::vector<VertexScores >);
   std::unique_ptr<std::vector<VertexScores > > VertexAssociationCollection(new std::vector<VertexScores >);
   edm::Handle<pat::MuonCollection> mu;
   edm::Handle<reco::VertexCollection>secondaryVertices;
   edm::Handle<reco::VertexCollection>primaryVertices;

   //get the products from the tags
   iEvent.getByToken(muonToken_, mu);
   iEvent.getByToken(tag_secondaryVertices_, secondaryVertices);
   iEvent.getByToken(tag_primaryVertices_, primaryVertices);

   const pat::MuonCollection&    muon = *(mu.product());
   const reco::VertexCollection &   sv = *(secondaryVertices.product());
   const reco::VertexCollection &   pv = *(primaryVertices.product());

   //use the first primary vertex for now
   const reco::Vertex primaryVertex = pv[0];
   if (debug_ > 1) std::cout << "[DEBUG] [SV Associator] Loop over muons " << std::endl;
   float Muon_pt= 999;
   float Muon_eta = 999;
   float Muon_phi = 999;

   // loop over each jet
   // reco::PFJetCollection::const_iterator jj = jet.begin();
   //for(; jj != jet.end() ; ++jj){
   pat::MuonCollection::const_iterator Mu = muon.begin();
     for(; Mu != muon.end() ; ++Mu){
       reco::TrackRef tunePTrack = Mu->muonBestTrack();
       if(tunePTrack->pt() < Muon_pt){
	 Muon_pt = tunePTrack->pt();
	 Muon_eta = tunePTrack->eta();
	 Muon_phi = tunePTrack->phi();
       }

     // minimum Muon pt
     if (Muon_pt < muonPtCut_) continue;

     //start building the vertex scores
     VertexScores vertexScores; 
     float bestScore = -1;
     reco::Vertex bestVertex;
  
   //loop over each of the vertices to be scored
     if (debug_ > 1) std::cout << "[DEBUG] [SV Associator] Loop over Vertices " << std::endl; 
     reco::VertexCollection::const_iterator  ss = sv.begin();
     size_t iss = 0;
     for(; ss != sv.end(); ++ss, ++iss) {
       float score = 0;             

       //loop over each track in the vertex       
       if (debug_ > 1) std::cout << "[DEBUG] [SV Associator] Loop over Tracks " << std::endl; 
       std::vector<reco::TrackBaseRef>::const_iterator tt = ss->tracks_begin();
       for(; tt != ss->tracks_end(); ++tt) {
	 float track_outerEta = (*tt)->outerEta();
	 float track_outerPhi = (*tt)->outerPhi();

	 float dR = reco::deltaR(Muon_eta, Muon_phi, track_outerEta, track_outerPhi);

	 if (dR > 1.0) continue;
	 else score += 1.0 / dR;    
       }

       if (debug_ > 1) std::cout << "[DEBUG] [SV Associator] Inserting Scores " << std::endl;   

       //parse a vertex reference using the handle and size
       //std::pair<reco::Vertex, float> scorePair(*ss, score);
       vertexScores.push_back(score);

       if (score > bestScore) {
	 bestScore = score;
	 bestVertex = *ss;
       }
     }

       if (debug_ > 1) std::cout << "[DEBUG] [SV Associator] Building Association" << std::endl; 
     // put together the scores with the best vertex 
     //VertexAssociation association(vertexScores, *jj, bestVertex, algoName_);     
           (*VertexAssociationCollection).push_back(vertexScores);         
	   //(*VertexAssociationCollection).push_back(bestVertex);

   }
   outputLabel_.clear();
   iEvent.put(std::move(VertexAssociationCollection), outputLabel_);
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::unique_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(std::move(pOut));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
DisplacedSVAssociator::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
DisplacedSVAssociator::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DisplacedSVAssociator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DisplacedSVAssociator::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DisplacedSVAssociator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DisplacedSVAssociator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedSVAssociator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedSVAssociator);

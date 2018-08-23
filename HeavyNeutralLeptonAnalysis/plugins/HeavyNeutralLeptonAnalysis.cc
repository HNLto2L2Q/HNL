// -*- C++ -*-
//
// Package:    HNL/HeavyNeutralLeptonAnalysis
// Class:      HeavyNeutralLeptonAnalysis
// 
/**\class HeavyNeutralLeptonAnalysis HeavyNeutralLeptonAnalysis.cc HNL/HeavyNeutralLeptonAnalysis/plugins/HeavyNeutralLeptonAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  root
//         Created:  Fri, 23 Mar 2018 12:09:13 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// for the plotting
#include "HNL/HeavyNeutralLeptonAnalysis/interface/SmartSelectionMonitor.h"

// root includes
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include <Math/VectorUtil.h>



//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoTauTag/TauTagTools/interface/GeneratorTau.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h"
#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"
#include "HNL/DisplacedSVAssociator/interface/VertexAssociation.h"

using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class HeavyNeutralLeptonAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit HeavyNeutralLeptonAnalysis(const edm::ParameterSet&);
  ~HeavyNeutralLeptonAnalysis();
  
  vector<reco::VertexCollection> PrimaryVertex( const reco::VertexCollection &vtx);
  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle);
  double MatchGenFirstMuon(const edm::Event&,  reco::TrackRef BestTrack);
  double MatchGenSecondMuon(const edm::Event&,  reco::TrackRef BestTrack);
  double MatchGenVertex(const edm::Event& iEvent, reco::Vertex vertex);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void initialize(const edm::Event&); 
  virtual void endJob() override;
  
  

  edm::Service<TFileService> fs;
  TTree * tree_;
  BigNtuple ntuple_;

  int    debug;
  float ivfMatchingScore;
  float pvCompatibilityScore = .05;
  bool  selIVFIsPV;

  //--------------Variables------------------------
  //======================= Primary Vertex Information ========================//                     
 
  //============= Trigger Information ===========================//  
  static const Int_t MAX_TRIGGERS      = 200;
  Int_t  nTrig;
  std::vector<std::string>   triggerNames;
  Int_t         triggerPass[MAX_TRIGGERS];
  Int_t           passIsoTk18 ;
  Int_t           passIsoTk20 ;
  Int_t           passIsoTk22 ;
  Int_t           passIsoTk24 ;
  Int_t           passIsoTk27 ;
  Int_t           passIsoTk17e;
  Int_t           passIsoTk22e;

  Int_t           passIsoMu18 ;
  Int_t           passIsoMu20 ;
  Int_t           passIsoMu22 ;
  Int_t           passIsoMu24 ;
  Int_t           passIsoMu27 ;
  Int_t           passIsoMu17e;
  Int_t           passIsoMu22e;
  Int_t           passTkMu17  ;
  Int_t           passTkMu20  ;


  Int_t           passIsoMu24All;
  Int_t           passIsoMu27All;

  Int_t           passDoubleMu17TrkIsoMu8    ;
  Int_t           passDoubleMu17TrkIsoTkMu8  ;
  Int_t           passDoubleTkMu17TrkIsoTkMu8;
  
  //============= Secondary Vertex Information ===========================//

  //============= Muon Information ===========================//  

  //============= JET Information ===========================// 

  //============= MET Information ===========================// 

  // ----------  Files  ---------------------------
  
  
  //--------------template-------------------------

  //std::vector<bool>  GoodPV;
  //
  //std::vector<int> NbGoodMuons;
  
  
  
  

  // ----------member data ---------------------------
  bool   isMC;
  
  edm::EDGetTokenT         < reco::VertexCollection > vtxMiniAODToken_;
  edm::EDGetTokenT                         < double > rhoToken_;
  edm::EDGetTokenT            < pat::MuonCollection > muonsMiniAODToken_;
  edm::EDGetTokenT        < pat::ElectronCollection > electronsMiniAODToken_;
  edm::EDGetTokenT           < EcalRecHitCollection > recHitEBToken_;
  edm::EDGetTokenT           < EcalRecHitCollection > recHitEEToken_;
  edm::EDGetTokenT             < pat::TauCollection > tausMiniAODToken_;
  edm::EDGetTokenT < pat::PackedCandidateCollection > packedCandidateToken_;
  edm::EDGetTokenT             < pat::JetCollection > jetsMiniAODToken_;
  edm::EDGetTokenT             < pat::METCollection > pfMETAODToken_;
  edm::EDGetTokenT            < edm::TriggerResults > triggerResultsToken_;
  edm::EDGetTokenT            < edm::TriggerResults > metFilterResultsToken_;
  edm::EDGetTokenT    < reco::GenParticleCollection > genParticleToken_;
  edm::EDGetTokenT            < GenEventInfoProduct > genEventInfoToken_;
  edm::EDGetTokenT < std::vector<PileupSummaryInfo> > PUInfoToken_;
  edm::EDGetTokenT                < LHEEventProduct > lheEventProductToken_;
  edm::EDGetTokenT         < reco::VertexCollection > inclusiveSecondaryVertices_;
  
protected:
  edm::Handle         < reco::VertexCollection > vtxHandle;
  edm::Handle                         < double > rhoHandle;
  edm::Handle            < pat::MuonCollection > muonsHandle;
  edm::Handle        < pat::ElectronCollection > electronsHandle;
  edm::Handle           < EcalRecHitCollection > recHitCollectionEBHandle;
  edm::Handle           < EcalRecHitCollection > recHitCollectionEEHandle;
  edm::Handle             < pat::TauCollection > tausHandle;
  edm::Handle < pat::PackedCandidateCollection > pfCandidatesHandle;
  edm::Handle             < pat::JetCollection > jetsHandle;
  edm::Handle             < pat::METCollection > metsHandle;
  edm::Handle            < edm::TriggerResults > triggerResultsHandle;
  edm::Handle            < edm::TriggerResults > metFilterResultsHandle;
  edm::Handle         < reco::VertexCollection > SecondaryVertex;
  
  /* Only for MC */
  edm::Handle    < reco::GenParticleCollection > genHandle;
  edm::Handle            < GenEventInfoProduct > genEventInfoHandle;
  edm::Handle < std::vector<PileupSummaryInfo> > puInfoH;
  edm::Handle                < LHEEventProduct > lheEPHandle;
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
HeavyNeutralLeptonAnalysis::HeavyNeutralLeptonAnalysis(const edm::ParameterSet& iConfig):
  ntuple_(),
  debug(iConfig.getParameter<int>("debugLevel")),
  isMC(iConfig.getParameter<bool>("isMC")),
  vtxMiniAODToken_(mayConsume<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxSrc"))),
  rhoToken_(mayConsume<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonsMiniAODToken_(mayConsume<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  electronsMiniAODToken_(mayConsume<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  recHitEBToken_(mayConsume< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("recHitCollectionEBSrc"))),
  recHitEEToken_(mayConsume< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("recHitCollectionEESrc"))),
  tausMiniAODToken_(mayConsume< pat::TauCollection >(iConfig.getParameter<edm::InputTag>("tauSrc"))),
  packedCandidateToken_(mayConsume< pat::PackedCandidateCollection >(iConfig.getParameter<edm::InputTag>("packCandSrc"))),
  jetsMiniAODToken_(mayConsume< pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  pfMETAODToken_(mayConsume<pat::METCollection>(iConfig.getParameter<edm::InputTag>("pfMETSrc"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResultSrc"))),
  metFilterResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterResultSrc"))),
  genParticleToken_(mayConsume<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  genEventInfoToken_(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
  PUInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
  lheEventProductToken_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts"))),
  inclusiveSecondaryVertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("SecondaryVertices")))
{

  //now do what ever initialization is needed
  usesResource("TFileService");

  tree_ = fs->make<TTree>("tree_", "tree");
  ntuple_.set_evtinfo(tree_);

}


HeavyNeutralLeptonAnalysis::~HeavyNeutralLeptonAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

void HeavyNeutralLeptonAnalysis::initialize(const edm::Event& iEvent){
  iEvent.getByToken(vtxMiniAODToken_, vtxHandle);
  iEvent.getByToken(rhoToken_, rhoHandle);
  iEvent.getByToken(muonsMiniAODToken_, muonsHandle);
  iEvent.getByToken(electronsMiniAODToken_, electronsHandle);
  iEvent.getByToken(recHitEBToken_, recHitCollectionEBHandle);
  iEvent.getByToken(recHitEEToken_, recHitCollectionEEHandle);
  iEvent.getByToken(tausMiniAODToken_, tausHandle);
  iEvent.getByToken(packedCandidateToken_, pfCandidatesHandle);
  iEvent.getByToken(jetsMiniAODToken_, jetsHandle);
  iEvent.getByToken(pfMETAODToken_, metsHandle);
  iEvent.getByToken(triggerResultsToken_, triggerResultsHandle);
  iEvent.getByToken(metFilterResultsToken_, metFilterResultsHandle);
  iEvent.getByToken(inclusiveSecondaryVertices_, SecondaryVertex);

  if (isMC){
    iEvent.getByToken(genParticleToken_, genHandle);
    iEvent.getByToken(genEventInfoToken_, genEventInfoHandle);
    iEvent.getByToken(PUInfoToken_, puInfoH);
    iEvent.getByToken(lheEventProductToken_, lheEPHandle);
  }
}



//
// member functions

std::vector<reco::VertexCollection> HeavyNeutralLeptonAnalysis::PrimaryVertex( const reco::VertexCollection &vtx)
{
  std::vector<reco::VertexCollection> allPVs;

    for(reco::VertexCollection::const_iterator PV = vtx.begin(); PV!=vtx.end();++PV) 
    {
       if(!PV->isFake()) {
         if(PV->ndof() > 4 && fabs(PV->position().z()) <= 24 && fabs(PV->position().rho()) <= 2 ) allPVs.push_back(*PV);
       }
    } 
  return  allPV;
}

//======================================================================================// 
bool HeavyNeutralLeptonAnalysis::isAncestor(const reco::Candidate* ancestor,const reco::Candidate* particle)
{
  if(ancestor == particle ) return true;
  for(size_t i=0;i< particle->numberOfMothers();i++)
    {
      if(isAncestor(ancestor,particle->mother(i))) return true;
    }
  return false;
}
//===================================== First Muon Gen Match ================================================//
double HeavyNeutralLeptonAnalysis::MatchGenFirstMuon(const edm::Event& iEvent, reco::TrackRef BestTrack) {
  double minDeltaR=999;
  double recoGenDeltaR_ = 0.2;
  int genIndex=-999;
    for(size_t i = 0; i<genHandle->size(); i++){
      if(abs((*genHandle)[i].pdgId()) == 13 && abs((*genHandle)[i].mother()->pdgId()) == 24) {
	const reco::Candidate * WBoson = (*genHandle)[i].mother();
	for(size_t j = 0; j<genHandle->size(); j++){
	  const reco::Candidate * part = (*genHandle)[j].mother(0) ;
	  if(part != nullptr && isAncestor( WBoson , part) ){
	    if( (*genHandle)[j].status() == 1 && abs((*genHandle)[j].pdgId()) == 13 ){
	      double dR=deltaR(BestTrack->eta(),BestTrack->phi(),(*genHandle)[j].eta(),(*genHandle)[j].phi());
	      if (dR<minDeltaR) {
		minDeltaR = dR;
		genIndex = j;
	      }// end of DR                                                                                                                        
	    } //end of status==1 && partID==13                                                                                                     
	  }//isAncestor                                                                                                                            
	}//end of loop over gen particles                                                                                                          
      }//end of Idetitfictaion of WBoson                                                                                                           
    }// end of loop over gen particles                                                                                                             
  if (minDeltaR<recoGenDeltaR_) return genIndex;
  else return -999;
}
//===================================== Second Muon Gen Match ================================================// 
double HeavyNeutralLeptonAnalysis::MatchGenSecondMuon(const edm::Event& iEvent,reco::TrackRef BestTrack) {
  double minDeltaR=999;
  double recoGenDeltaR_ = 0.2;
  int genIndex=-999;
    for(size_t i = 0; i<genHandle->size(); i++){
      if(abs((*genHandle)[i].pdgId()) == 13 && (*genHandle)[i].mother()->pdgId() == 9900014) {
	const reco::Candidate * HNL = (*genHandle)[i].mother();
	for(size_t j = 0; j<genHandle->size(); j++){
	  const reco::Candidate * part = (*genHandle)[j].mother(0) ;
	  if(part != nullptr && isAncestor( HNL , part) ){
	    if( (*genHandle)[j].status() == 1 && abs((*genHandle)[j].pdgId()) == 13 ){
	      double dR=deltaR(BestTrack->eta(),BestTrack->phi(),(*genHandle)[j].eta(),(*genHandle)[j].phi());
	      if (dR<minDeltaR) {
		minDeltaR = dR;
		genIndex=j;
	      }// end of DR   
	    } //end of status==1 && partID==13                                                                                                     
	  }//isAncestor                                                                                                                            
	}//end of loop over gen particles                                                                                                          
      }//end of Idetitfictaion of HNL                                                                                                              
    }// end of loop over gen particles                                                                                                              
  if (minDeltaR<recoGenDeltaR_) return genIndex;
  else return -999;
}
//===================================== Displaced Vertex Gen Match ================================================// 
double HeavyNeutralLeptonAnalysis::MatchGenVertex(const edm::Event& iEvent, reco::Vertex vertex) {
  double minDVertex=999;
  double MinDistance_ = 0.2;
  int genIndex=-999;
  for(size_t i = 0; i<genHandle->size(); i++){
    if(abs((*genHandle)[i].pdgId()) == 13 && (*genHandle)[i].mother()->pdgId() == 9900014) {
      const reco::Candidate * HNL = (*genHandle)[i].mother();
      for(size_t j = 0; j<genHandle->size(); j++){
	const reco::Candidate * part = (*genHandle)[j].mother(0) ;
	if(part != nullptr && isAncestor( HNL , part) ){
	  if( (*genHandle)[j].status() == 1){
	    float vx = vertex.x(), vy = vertex.y(), vz = vertex.z();
	    float gx =  (*genHandle)[j].vx(), gy =  (*genHandle)[j].vy(), gz =  (*genHandle)[j].vz();
	    //float metric1 = std::sqrt(((gx-vx)*(gx-vx))/(gx*gx) + ((gy-vy)*(gy-vy))/(gy*gy) + ((gz-vz)*(gz-vz))/(gz*gz)); 
	    float metric = std::sqrt(((gx-vx)*(gx-vx)) + ((gy-vy)*(gy-vy)) + ((gz-vz)*(gz-vz)));
	    if (metric<minDVertex) {
	      minDVertex= metric;
	      genIndex=i;
	    }// end of DR                                                                                                                        
	  }//end of status==1 && partID==13                                                                                                     
	}//isAncestor                                                                                                                            
      }//end of loop over gen particles                                                                                                          
    }//end of loop over gen particles  
  }//end of Idetitfictaion of HNL
  if (minDVertex<MinDistance_) return genIndex;
  else return -999;
}


// ------------ method called for each event  ------------
void HeavyNeutralLeptonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  
  initialize(iEvent);
  
  ntuple_.fill_evtInfo(iEvent.id());
  
  //============================================================= 
  //
  //             Primary Vertex
  // 
  //=============================================================   
   reco::VertexCollection vtx;
  
  if(vtxHandle.isValid()){
    vtx = *vtxHandle;
    std::vector<reco::VertexCollection> pvs = PrimaryVertex(vtx);
    ntuple_.fill_evtInfo(pvs);
  }
  
  else continue;

  //=============================================================
   //
   //             Trigger 
   //    
   //=============================================================     
   const edm::TriggerResults triggerResults =  *triggerResultsHandle.product();
   const edm::TriggerNames&    trigNames  = iEvent.triggerNames(triggerResults);
   nTrig = 0;

   passIsoTk18  = 0;
   passIsoTk20  = 0;
   passIsoTk22  = 0;
   passIsoTk24  = 0;
   passIsoTk27  = 0;
   passIsoTk17e = 0;
   passIsoTk22e = 0;

   passIsoMu18  = 0;
   passIsoMu20  = 0;
   passIsoMu22  = 0;
   passIsoMu24  = 0;
   passIsoMu27  = 0;
   passIsoMu17e = 0;
   passIsoMu22e = 0;
   passTkMu17   = 0;
   passTkMu20   = 0;

   passIsoMu24All = 0;
   passIsoMu27All = 0;

   passDoubleMu17TrkIsoMu8     = 0;
   passDoubleMu17TrkIsoTkMu8   = 0;
   passDoubleTkMu17TrkIsoTkMu8 = 0;

   for (size_t i = 0; i < trigNames.size(); ++i) {
     const std::string &name = trigNames.triggerName(i);
     bool fired = triggerResults.accept(i);
     
     if(!fired) continue;
     // specific triggers iso muons                                                                                                                  
     std::size_t searchHTIsoMu18    = name.find("HLT_IsoMu18_v");
     std::size_t searchHTIsoMu20    = name.find("HLT_IsoMu20_v");
     std::size_t searchHTIsoMu22    = name.find("HLT_IsoMu22_v");                                                                                  
     std::size_t searchHTIsoMu24    = name.find("HLT_IsoMu24_v");
     std::size_t searchHTIsoMu27    = name.find("HLT_IsoMu27_v");
     std::size_t searchHTIsoMu17e   = name.find("HLT_IsoMu17_eta2p1_v");
     std::size_t searchHTIsoMu22e   = name.find("HLT_IsoMu22_eta2p1_v");                                                                           

     //ISO TRACK                                                                                                                                     
     std::size_t searchHTIsoTk18          = name.find("HLT_IsoTkMu18_v");
     std::size_t searchHTIsoTk20          = name.find("HLT_IsoTkMu20_v");
     std::size_t searchHTIsoTk22          = name.find("HLT_IsoTkMu22_v");                                                                          
     std::size_t searchHTIsoTk24          = name.find("HLT_IsoTkMu24_v");
     std::size_t searchHTIsoTk27          = name.find("HLT_IsoTkMu27_v");
     std::size_t searchHTIsoTk17e         = name.find("HLT_IsoTkMu17_eta2p1_v");
     std::size_t searchHTIsoTk22e         = name.find("HLT_IsoTkMu22_eta2p1_v");                                                                   

     std::size_t searchHTTkMu17           = name.find("HLT_TkMu17_v");
     std::size_t searchHTTkMu20           = name.find("HLT_TkMu20_v");

     std::size_t searchDoubleMu17TrkIsoMu8     = name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
     std::size_t searchDoubleMu17TrkIsoTkMu8   = name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
     std::size_t searchDoubleTkMu17TrkIsoTkMu8 = name.find("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

     bool    htist_18                 = searchHTIsoTk18               != std::string::npos ;
     bool    htist_20                 = searchHTIsoTk20               != std::string::npos ;
     bool    htist_22                 = searchHTIsoTk22               != std::string::npos ;                                                         
     bool    htist_24                 = searchHTIsoTk24               != std::string::npos ;
     bool    htist_27                 = searchHTIsoTk27               != std::string::npos ;
     bool    htist_17e                = searchHTIsoTk17e              != std::string::npos ;
     bool    htist_22e                = searchHTIsoTk22e              != std::string::npos ;                                                         

     bool    htism_18                 = searchHTIsoMu18               != std::string::npos ;
     bool    htism_20                 = searchHTIsoMu20               != std::string::npos ;
     bool    htism_22                 = searchHTIsoMu22               != std::string::npos ;                                                         
     bool    htism_24                 = searchHTIsoMu24               != std::string::npos ;
     bool    htism_27                 = searchHTIsoMu27               != std::string::npos ;
     bool    htism_17e                = searchHTIsoMu17e              != std::string::npos ;
     bool    htism_22e                = searchHTIsoMu22e              != std::string::npos ;                                                         

     bool    httkm_17                 = searchHTTkMu17                != std::string::npos ;
     bool    httkm_20                 = searchHTTkMu20                != std::string::npos ;

     bool    DoubleMu17TrkIsoMu8      = searchDoubleMu17TrkIsoMu8     != std::string::npos ;
     bool    DoubleMu17TrkIsoTkMu8    = searchDoubleMu17TrkIsoTkMu8   != std::string::npos ;
     bool    DoubleTkMu17TrkIsoTkMu8  = searchDoubleTkMu17TrkIsoTkMu8 != std::string::npos ;

     passIsoTk18  =  passIsoTk18    || htist_18;
     passIsoTk20  =  passIsoTk20    || htist_20;
     passIsoTk22  =  passIsoTk22    || htist_22;   
     passIsoTk24  =  passIsoTk24    || htist_24;
     passIsoTk27  =  passIsoTk27    || htist_27;
     passIsoTk17e =  passIsoTk17e   || htist_17e;
     passIsoTk22e =  passIsoTk22e   || htist_22e;  

     passIsoMu18  = passIsoMu18    || htism_18;
     passIsoMu20  = passIsoMu20    || htism_20;
     passIsoMu22  = passIsoMu22    || htism_22; 
     passIsoMu24  = passIsoMu24    || htism_24;
     passIsoMu27  = passIsoMu27    || htism_27;
     passIsoMu17e = passIsoMu17e   || htism_17e;
     passIsoMu22e = passIsoMu22e   || htism_22e;                             

     passTkMu17   = passTkMu17     || httkm_17;
     passTkMu20   = passTkMu20     || httkm_20;


     passIsoMu24All = passIsoMu24All   || passIsoMu24 || passIsoTk24 ;
     passIsoMu27All = passIsoMu27All   || passIsoMu27 || passIsoTk27 ;

     passDoubleMu17TrkIsoMu8     = passDoubleMu17TrkIsoMu8     ||     DoubleMu17TrkIsoMu8 ;   
     passDoubleMu17TrkIsoTkMu8   = passDoubleMu17TrkIsoTkMu8   ||     DoubleMu17TrkIsoTkMu8 ; 
     passDoubleTkMu17TrkIsoTkMu8 = passDoubleTkMu17TrkIsoTkMu8 ||     DoubleTkMu17TrkIsoTkMu8 ; 

   }

   //=============================================================
   //
   //                Secondary Vertex
   //
   //============================================================= 
    
   reco::VertexCollection sv;
   if(SecondaryVertex.isValid() && muonsHandle.isValid()){ 

     sv = *SecondaryVertex;

     std::vector<reco::Vertex> SV;
     SV.insert(SV.end(), sv.begin(), sv.end());
     
     VertexAssociation JVAIVF("IVF", pvs.front(), debug);

     reco::VertexCollection::const_iterator vtxIter = SV.begin();
     for(; vtxIter != SV.end(); ++vtxIter ) {
       JVAIVF.addVertex(*vtxIter);
     }
     
     const std::pair<reco::Vertex, float>    bestVertexPair  = JVAIVF.getBestVertex(muonsHandle, "oneOverR");
     const reco::Vertex                      bestVertex      = bestVertexPair.first;
     const float                             bestVertexScore = bestVertexPair.second;

     _ntuple.fill_svInfo(bestVertex,bestVertexScore)


       //to add to fill_svInfo
       // double GenParticleIndex = -99.99;
       //if(isMC)  GenParticleIndex = MatchGenVertex(iEvent,  bestVertex);
       //VertexMatch.push_back(GenParticleIndex);

   //=============================================================                                                                                   
   //                                                                                                                                                
   //                Method for Muon Tree                                           
   //                                                                                                                                                
   //=============================================================                                                                                   
   Muon_SecondGenMatch.clear();
  // NbGoodMuons.clear();

   int NbMuons = 0;
   std::string muon;
   pat::MuonCollection muons;
   if(muonsHandle.isValid()){ muons = *muonsHandle;
     const reco::Vertex &pv = vtxHandle->front();
     for (const pat::Muon mu : muons) {
     if( mu.pt() < 0.0 ) continue;
     if (!( fabs(mu.eta()) < 2.4 && mu.pt() > 5. )) continue;
     ++NbMuons;
       //general muons infos//
     Muon_isGlobalMuon.push_back(mu.isGlobalMuon());
     Muon_isPF.push_back(mu.isPFMuon());
     Muon_isTrackerMuon.push_back(mu.isTrackerMuon());
     Muon_isRPCMuon.push_back(mu.isRPCMuon());
     Muon_isStandAloneMuon.push_back(mu.isStandAloneMuon());
     Muon_isSoftMuon.push_back(mu.isSoftMuon(pv));
     Muon_isLoose.push_back(mu.isLooseMuon());
     Muon_isTightMuon.push_back(mu.isTightMuon(pv));
     //Muon_isGoodMuon.push_back(mu.isGood(muon));
     Muon_en.push_back(mu.energy());
     Muon_et.push_back(mu.et());
     Muon_pt.push_back(mu.pt());
     Muon_eta.push_back(mu.eta());
     Muon_phi.push_back(mu.phi());
     Muon_charge.push_back(mu.charge());
     
       // best track and high quality//
     reco::TrackRef tunePTrack = mu.muonBestTrack();
     Muon_nbMuon.push_back(NbMuons);
     Muon_ptTunePMuonBestTrack.push_back(tunePTrack->pt()); // transverse momentum                                                                
     Muon_dPToverPTTunePMuonBestTrack.push_back(tunePTrack->ptError()/tunePTrack->pt()); // error calculation of transverse momentum              
     Muon_pxTunePMuonBestTrack.push_back(tunePTrack->px()); //px component of the track                                                           
     Muon_pyTunePMuonBestTrack.push_back(tunePTrack->py()); //py component of the track                                                           
     Muon_pzTunePMuonBestTrack.push_back(tunePTrack->pz()); //pz component of the track                                                           
     Muon_pTunePMuonBestTrack.push_back(tunePTrack->p());   //magnitude of momentum vector                                                        
     Muon_etaTunePMuonBestTrack.push_back(tunePTrack->eta());
     Muon_phiTunePMuonBestTrack.push_back(tunePTrack->phi());
     Muon_thetaTunePMuonBestTrack.push_back(tunePTrack->theta());
     Muon_chargeTunePMuonBestTrack.push_back(tunePTrack->charge());
     Muon_absdxyTunePMuonBestTrack.push_back(fabs(tunePTrack->dxy(pv.position()))); //transvers  impact parameter  w.r.t. the primary vertex      
     Muon_absdxyErrorTunePMuonBestTrack.push_back(fabs(tunePTrack->dxyError())); //transvers  impact parameter  w.r.t. the primary vertex         
     Muon_absdxySigTunePMuonBestTrack.push_back(fabs(tunePTrack->dxy(pv.position()))/fabs(tunePTrack->dxyError()));
     Muon_absdzTunePMuonBestTrack.push_back(fabs(tunePTrack->dz(pv.position()))); // longitudinal impact parameter  w.r.t. the primary vertex     
     Muon_absdzErrorTunePMuonBestTrack.push_back(fabs(tunePTrack->dzError())); // longitudinal impact parameter  w.r.t. the primary vertex  
     Muon_absdzSigTunePMuonBestTrack.push_back(fabs(tunePTrack->dz(pv.position()))/fabs(tunePTrack->dzError())); 
     Muon_TrackQuality.push_back(tunePTrack->quality(reco::TrackBase::highPurity));

     int GenParticleIndex1 = -99;
     if(isMC )  GenParticleIndex1 = MatchGenFirstMuon(iEvent, tunePTrack);
     Muon_FirstGenMatch.push_back(GenParticleIndex1);

     double MatchParticleIndex = -99.99;
     if(isMC )  MatchParticleIndex = MatchGenSecondMuon(iEvent, tunePTrack);
     Muon_SecondGenMatch.push_back(MatchParticleIndex);
     

     if(mu.globalTrack().isNonnull() ) {
       Muon_normalizedChi2.push_back(mu.globalTrack()->normalizedChi2());
       Muon_numberOfValidPixelHits.push_back(mu.innerTrack()->hitPattern().numberOfValidPixelHits());
       Muon_numberOfValidMuonHits.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
       Muon_numberOftrackerLayersWithMeasurement.push_back(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
       Muon_numberOfMatchedStations.push_back(mu.numberOfMatchedStations());
       Muon_numberOfpixelLayersWithMeasurement.push_back(mu.innerTrack()->hitPattern().pixelLayersWithMeasurement());
       Muon_InnerTrackQuality.push_back(mu.innerTrack()->quality(reco::TrackBase::highPurity));
     }
     else{
       Muon_normalizedChi2.push_back(-999);
       Muon_numberOfValidPixelHits.push_back(-999);
       Muon_numberOfValidMuonHits.push_back(-999);
       Muon_numberOftrackerLayersWithMeasurement.push_back(-999);
       Muon_numberOfMatchedStations.push_back(-999);
       Muon_numberOfpixelLayersWithMeasurement.push_back(-999);
       Muon_InnerTrackQuality.push_back(-999);
     }
     
     if(mu.standAloneMuon().isNonnull() ) {
       Muon_STAnHits.push_back(mu.standAloneMuon()->numberOfValidHits());                                                                        
       Muon_STAnLost.push_back(mu.standAloneMuon()->numberOfLostHits());                                                                         
       Muon_STAnStationsWithAnyHits.push_back(mu.standAloneMuon()->hitPattern().muonStationsWithAnyHits());                                      
       Muon_STAnCscChambersWithAnyHits.push_back(mu.standAloneMuon()->hitPattern().cscStationsWithAnyHits()); //csc chambers in track fit        
       Muon_STAnDtChambersWithAnyHits.push_back(mu.standAloneMuon()->hitPattern().dtStationsWithAnyHits()); //dt chambers in track fit           
       Muon_STAnRpcChambersWithAnyHits.push_back(mu.standAloneMuon()->hitPattern().rpcStationsWithAnyHits()); //rpc chambers in track fit        
       Muon_STAinnermostStationWithAnyHits.push_back(mu.standAloneMuon()->hitPattern().innermostMuonStationWithAnyHits());                       
       Muon_STAoutermostStationWithAnyHits.push_back(mu.standAloneMuon()->hitPattern().outermostMuonStationWithAnyHits());                       
       Muon_STAnCscChambersWithValidHits.push_back(mu.standAloneMuon()->hitPattern().cscStationsWithValidHits()); //csc chambers anywhere near track
       
       Muon_STAnDtChambersWithValidHit.push_back(mu.standAloneMuon()->hitPattern().dtStationsWithValidHits()); //dt chambers anywhere near track 
       Muon_STAnRpcChambersWithValidHits.push_back(mu.standAloneMuon()->hitPattern().rpcStationsWithValidHits()); //rpc chambers anywhere near track
       Muon_STAnValidCscHits.push_back(mu.standAloneMuon()->hitPattern().numberOfValidMuonCSCHits()); //CSC hits anywhere near track             
       Muon_STAnValidDtHits.push_back(mu.standAloneMuon()->hitPattern().numberOfValidMuonDTHits()); //DT hits anywhere near track                
       Muon_STAnValidRpcHits.push_back(mu.standAloneMuon()->hitPattern().numberOfValidMuonRPCHits()); //RPC hits anywhere near track             
       Muon_STAnValidMuonHits.push_back(mu.standAloneMuon()->hitPattern().numberOfValidMuonHits()); //muon hits anywhere near track              
       Muon_STAinnermostStationWithValidHits.push_back(mu.standAloneMuon()->hitPattern().innermostMuonStationWithValidHits());                   
       Muon_STAoutermostStationWithValidHits.push_back(mu.standAloneMuon()->hitPattern().outermostMuonStationWithValidHits());                   
       Muon_STAnStationsWithValidHits.push_back(mu.standAloneMuon()->hitPattern().muonStationsWithValidHits());
     }

     else{
       Muon_STAnHits.push_back(-999);
       Muon_STAnLost.push_back(-999);
       Muon_STAnStationsWithAnyHits.push_back(-999);
       Muon_STAnCscChambersWithAnyHits.push_back(-999);
       Muon_STAnDtChambersWithAnyHits.push_back(-999);
       Muon_STAnRpcChambersWithAnyHits.push_back(-999);
       Muon_STAinnermostStationWithAnyHits.push_back(-999);
       Muon_STAoutermostStationWithAnyHits.push_back(-999);
       Muon_STAnCscChambersWithValidHits.push_back(-999);
       Muon_STAnDtChambersWithValidHit.push_back(-999);
       Muon_STAnRpcChambersWithValidHits.push_back(-999);
       Muon_STAnValidCscHits.push_back(-999);
       Muon_STAnValidDtHits.push_back(-999);
       Muon_STAnValidRpcHits.push_back(-999);
       Muon_STAnValidMuonHits.push_back(-999);
       Muon_STAinnermostStationWithValidHits.push_back(-999);
       Muon_STAoutermostStationWithValidHits.push_back(-999);
       Muon_STAnStationsWithValidHits.push_back(-999);
     }
     reco::MuonTime tofAll = mu.time();
     Muon_STATofDirection.push_back(tofAll.direction());
     Muon_STATofNDof.push_back(tofAll.nDof);
     Muon_STATofTimeAtIpInOut.push_back(tofAll.timeAtIpInOut);
     Muon_STATofTimeAtIpInOutErr.push_back(tofAll.timeAtIpInOutErr);
     Muon_STATofTimeAtIpOutIn.push_back(tofAll.timeAtIpOutIn);
     Muon_STATofTimeAtIpOutInErr.push_back(tofAll.timeAtIpOutInErr);      
     }
     //  NbGoodMuons.push_back(NbMuons);
   cout<<"nb of muons"<< NbMuons<<endl;
   }
   
   
   
   EcalRecHitCollection recHitCollectionEB;
   if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}

   EcalRecHitCollection recHitCollectionEE;
   if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}

   pat::TauCollection taus;
   if(tausHandle.isValid()){ taus = *tausHandle;}

   pat::PackedCandidateCollection pfCandidates;
   if (pfCandidatesHandle.isValid()) { pfCandidates = *pfCandidatesHandle; }

   
   pat::ElectronCollection electrons;
   if(electronsHandle.isValid()){ electrons = *electronsHandle;}

   //============================================================= 
   //
   //             Jets 
   //       
   //=============================================================

   int jetnumber=0;
   pat::JetCollection jets;

   ntuples_.ill_jetInfo(jets);
   if(jetsHandle.isValid()){ jets = *jetsHandle;
     ntuples_.ill_jetInfo(jets);
   }
   //=============================================================
   //
   //            Missing Energy 
   //     
   //============================================================= 
   
   pat::METCollection mets;
   if(metsHandle.isValid()){ mets = *metsHandle;
     const pat::MET met = mets.front();
     ntuple_.fill_metInfo(met);
     
   }
   
   }
   
}
   
// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNeutralLeptonAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNeutralLeptonAnalysis::endJob() 
{
//outputFile_->cd();
//tree_->Write();
//outputFile_->Close(); 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HeavyNeutralLeptonAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNeutralLeptonAnalysis);

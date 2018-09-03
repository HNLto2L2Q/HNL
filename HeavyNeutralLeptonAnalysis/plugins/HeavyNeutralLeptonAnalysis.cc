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

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h"
#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"

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
  
  reco::VertexCollection getMatchedVertex(const pat::Muon & mu, const reco::VertexCollection& vertexCollection);
  reco::VertexCollection PrimaryVertex( const reco::VertexCollection &vtx);
  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle);
  double MatchGenMuon(const edm::Event&,  reco::TrackRef BestTrack , int pdgId);
  double MatchGenElectron(const edm::Event& iEvent, const pat::Electron& ele_, int pdgId);
  double MatchGenVertex(const edm::Event& iEvent, reco::Vertex vertex, int pdgId);
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
  const std::vector<std::string> bDiscriminators_;
  
protected:
  edm::Handle         < reco::VertexCollection > vtxHandle;
  edm::Handle                         < double > rhoHandle;
  edm::Handle            < pat::MuonCollection > muonsHandle;
  edm::Handle       <std::vector<pat::Electron>> electronsHandle;
  edm::Handle           < EcalRecHitCollection > recHitCollectionEBHandle;
  edm::Handle           < EcalRecHitCollection > recHitCollectionEEHandle;
  edm::Handle             < pat::TauCollection > tausHandle;
  edm::Handle < pat::PackedCandidateCollection > pfCandidatesHandle;
  edm::Handle             < pat::JetCollection > jetsHandle;
  edm::Handle             < pat::METCollection > metsHandle;
  edm::Handle            < edm::TriggerResults > triggerResultsHandle;
  edm::Handle            < edm::TriggerResults > metFilterResultsHandle;
  edm::Handle         < reco::VertexCollection > secondaryVertexHandle;
  
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
  inclusiveSecondaryVertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("SecondaryVertices"))),
  bDiscriminators_(iConfig.getParameter<std::vector<std::string> >("bDiscriminators"))
{

  //now do what ever initialization is needed
  usesResource("TFileService");

  tree_ = fs->make<TTree>("tree_", "tree");
  ntuple_.set_evtInfo(tree_);
  ntuple_.set_jetInfo(tree_);
  ntuple_.set_metInfo(tree_);
  ntuple_.set_trigInfo(tree_);
  ntuple_.set_pvInfo(tree_);
  ntuple_.set_svInfo(tree_);
  ntuple_.set_muInfo(tree_);
  ntuple_.set_eleInfo(tree_);
  ntuple_.set_bjetInfo(tree_);
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
  iEvent.getByToken(inclusiveSecondaryVertices_, secondaryVertexHandle);

  if (isMC){
    iEvent.getByToken(genParticleToken_, genHandle);
    iEvent.getByToken(genEventInfoToken_, genEventInfoHandle);
    iEvent.getByToken(PUInfoToken_, puInfoH);
    iEvent.getByToken(lheEventProductToken_, lheEPHandle);
  }
}



//
// member functions

reco::VertexCollection HeavyNeutralLeptonAnalysis::getMatchedVertex(const pat::Muon & muon, const reco::VertexCollection& vertexCollection){
  reco::VertexCollection  matchedVertices;
  //cout<<muon.pfCandidateRef().isNull() << "  " <<muon.pfCandidateRef().isAvailable()<<endl;
  //cout << "Orig PTR: Num " << muon.numberOfSourceCandidatePtrs() 
  //     << " first : " << muon.sourceCandidatePtr(0).isAvailable() << " " << muon.sourceCandidatePtr(0).isNonnull() << endl;
  const pat::PackedCandidate* cand = dynamic_cast<const pat::PackedCandidate*>(muon.sourceCandidatePtr(0).get());
  //cout << "Packed: " << cand << endl;
  if(!cand) {
    cout << "THIS SHOULD NEVER HAPPEN! No packed candidated associated to muon?!" << endl;
  }
  //cout << "Packed: " << cand->hasTrackDetails() << " " << (cand->charge() != 0) << " " <<  (cand->numberOfHits() > 0) << endl;
  if(!(cand->hasTrackDetails() && cand->charge() != 0 && cand->numberOfHits() > 0)) {
    cout << "THIS SHOULD NEVER HAPPEN! Muon without track or matched to neutral?!" << endl;
  }
  //cout << cand->pseudoTrack().pt() << " " << cand->pseudoTrack().eta() << " " << cand->pseudoTrack().phi() << endl;

  for(reco::VertexCollection::const_iterator ss = vertexCollection.begin(); ss != vertexCollection.end(); ++ss) {    
    cout <<"new vertex"<<endl;
    for(reco::Vertex::trackRef_iterator tt = ss->tracks_begin(); tt != ss->tracks_end(); ++tt) {
      //cout<<"Track " << (*tt)->pt() << "  "<< (*tt)->eta()<< " " << (*tt)->phi() <<endl;
      //cout << "match " << (cand->pseudoTrack().pt() == tt->castTo<reco::TrackRef>()->pt()) << endl;
      if((cand->pseudoTrack().pt() == tt->castTo<reco::TrackRef>()->pt())) { //Options here: innerTrack, globalTrack, muonBestTrack, outerTrack, pickyTrack, track
        matchedVertices.push_back(*ss);
	break;
      }
    }
  } 
  return matchedVertices;
}

reco::VertexCollection HeavyNeutralLeptonAnalysis::PrimaryVertex( const reco::VertexCollection &vtx)
{
  reco::VertexCollection allPVs;

  for(reco::VertexCollection::const_iterator PV = vtx.begin(); PV!=vtx.end();++PV) 
    {
      if(!PV->isFake()) {
	if(PV->ndof() > 4 && fabs(PV->position().z()) <= 24 && fabs(PV->position().rho()) <= 2 ) allPVs.push_back(*PV);
      }
    } 
  return  allPVs;
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
//===================================== Muon Gen Match ================================================//
double HeavyNeutralLeptonAnalysis::MatchGenMuon(const edm::Event& iEvent, reco::TrackRef BestTrack, int pdgId) {
  double minDeltaR=999;
  double recoGenDeltaR_ = 0.2;
  int genIndex=-999;
    for(size_t i = 0; i<genHandle->size(); i++){
      if(abs((*genHandle)[i].pdgId()) == 13 && abs((*genHandle)[i].mother()->pdgId()) == pdgId) {
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
//===================================== Electron Gen Match ================================================// 
double HeavyNeutralLeptonAnalysis::MatchGenElectron(const edm::Event& iEvent, const pat::Electron& ele_, int pdgId) {
  double minDeltaR=999;
  double recoGenDeltaR_ = 0.2;
  int genIndex=-999;
  for(size_t i = 0; i<genHandle->size(); i++){
    if(abs((*genHandle)[i].pdgId()) == 11 && abs((*genHandle)[i].mother()->pdgId()) == pdgId) {
      const reco::Candidate * WBoson = (*genHandle)[i].mother();
      for(size_t j = 0; j<genHandle->size(); j++){
	const reco::Candidate * part = (*genHandle)[j].mother(0) ;
	if(part != nullptr && isAncestor( WBoson , part) ){
	  if( (*genHandle)[j].status() == 1 && abs((*genHandle)[j].pdgId()) == 11 ){
	    double dR=deltaR(ele_.eta(),ele_.phi(),(*genHandle)[j].eta(),(*genHandle)[j].phi());
	    if (dR<minDeltaR) {
	      minDeltaR = dR;
	      genIndex = j;
	    }// end of DR
	  } //end of status==1 && partID==11
	}//isAncestor
      }//end of loop over gen particles
    }//end of Idetitfictaion of WBoson  
  }// end of loop over gen particles                                            
  if (minDeltaR<recoGenDeltaR_) return genIndex;
  else return -999;
}

//===================================== Displaced Vertex Gen Match ================================================// 
double HeavyNeutralLeptonAnalysis::MatchGenVertex(const edm::Event& iEvent, reco::Vertex vertex, int pdgId) {
  double minDVertex=999;
  double MinDistance_ = 0.2;
  int genIndex=-999;
  for(size_t i = 0; i<genHandle->size(); i++){
    if(abs((*genHandle)[i].pdgId()) == pdgId && (*genHandle)[i].mother()->pdgId() == 9900012) {
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

  ntuple_.reset();
  ntuple_.fill_evtInfo(iEvent.id());
  
  //============================================================= 
  //
  //             Primary Vertex
  // 
  //=============================================================   
  if(!vtxHandle.isValid()) return;
  reco::VertexCollection pvs = PrimaryVertex(*vtxHandle);  
  if(!pvs.size()) return;
  ntuple_.fill_pvInfo(pvs);    
  //=============================================================
   //
   //             Trigger 
   //    
   //=============================================================     
   const edm::TriggerResults triggerResults =  *triggerResultsHandle.product();
   const edm::TriggerNames&    trigNames  = iEvent.triggerNames(triggerResults);
   ntuple_.fill_trigInfo(triggerResults, trigNames);    
   //=============================================================
   //
   //                Method for Muon Tree
   //
   //=============================================================
   std::string muon;
   pat::MuonCollection muons;

   vector<pat::Muon> goodMuons;
   vector<pat::Muon> looseMuons;
  
   if(muonsHandle.isValid()){ 
     muons = *muonsHandle;
     for (const pat::Muon mu : muons) {
       if( mu.pt() < 0.0 ) continue;
       if (!( fabs(mu.eta()) < 2.4 && mu.pt() > 5. )) continue;
       goodMuons.push_back(mu);
       if (!mu.isLooseMuon()) continue;
       looseMuons.push_back(mu);
     }
   }

   for (const pat::Muon mu : goodMuons){
     reco::TrackRef bestTrack = mu.muonBestTrack();
     double matching_1stmu = (isMC) ? MatchGenMuon(iEvent, bestTrack, 24) : -999;
     double matching_2ndmu = (isMC) ? MatchGenMuon(iEvent, bestTrack, 9900012) : -999;
     //added the matching
     ntuple_.fill_muInfo(mu, pvs.at(0) , matching_1stmu , matching_2ndmu);
   }

   // lambda function to sort this muons
   std::sort(looseMuons.begin(), looseMuons.end(), [](pat::Muon a, pat::Muon b) {return a.pt() < b.pt(); });
   if(!goodMuons.size()) return;
   if(!looseMuons.size()) return;
   pat::Muon muonHNL = looseMuons[0];
   //=============================================================                    
   //                     
   //                Secondary Vertex                     
   //                     
   //=============================================================                
   if(secondaryVertexHandle.isValid()){

     reco::VertexCollection bestVertices  = getMatchedVertex(muonHNL, *secondaryVertexHandle);

     // check if SV doesn't match with the PV
     for (const reco::Vertex& vtx : bestVertices){
       float x  = vtx.x(), y = vtx.y(), z = vtx.z();
       float dx = x - pvs.at(0).x() , dy = y - pvs.at(0).y(), dz = z - pvs.at(0).z();
       
       float  selIVFIsPVScore = std::sqrt((dx/x)*(dx/x) + (dy/y)*(dy/y) + (dz/z)*(dz/z));       
       if (selIVFIsPVScore < pvCompatibilityScore) continue;
       double matching_vtx = (isMC) ? MatchGenVertex(iEvent, vtx , 13) : -999;

       ntuple_.fill_svInfo(vtx, pvs.at(0), matching_vtx);	 

     }

     //////////////////////////////////////////////   
     EcalRecHitCollection recHitCollectionEB;
   if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}
   
   EcalRecHitCollection recHitCollectionEE;
   if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}
   //============================================================= 
   //
   //                Method for electrons
   //                                                   
   //============================================================= 
   // using iterator to get ref for ele
   for(auto ele = electronsHandle->begin(); ele != electronsHandle->end(); ++ele){
     if(!(!ele->gsfTrack().isNull() && ele->pt() > 5 && fabs(ele->eta()) < 2.5 ))      continue;
     
     double rho = *(rhoHandle.product());
     
     auto eleRef = edm::Ref<std::vector<pat::Electron>>(electronsHandle, (ele - electronsHandle->begin()));
     
     std::auto_ptr<EcalClusterLazyTools> recHitEcal;
     recHitEcal.reset(new EcalClusterLazyTools( iEvent, iSetup, recHitEBToken_, recHitEEToken_ ));
     
     double matching_1stele = (isMC) ? MatchGenElectron(iEvent, *ele , 24 ) : -999;
     double matching_2ndele = (isMC) ? MatchGenElectron(iEvent, *ele , 9900012) : -999;

     ntuple_.fill_eleInfo(*ele, pvs.at(0), rho , matching_1stele, matching_2ndele, recHitEcal);
   }

   pat::TauCollection taus;
   if(tausHandle.isValid()){ taus = *tausHandle;}

   pat::PackedCandidateCollection pfCandidates;
   if (pfCandidatesHandle.isValid()) { pfCandidates = *pfCandidatesHandle; }

   //============================================================= 
   //
   //             Jets 
   //       
   //=============================================================
   pat::JetCollection jets;
   if(jetsHandle.isValid()){ 
     jets = *jetsHandle;
     for (const pat::Jet jet : jets) {
       if( jet.pt() < 0.0 ) continue;
       if (!( fabs(jet.eta()) < 3 && jet.pt() > 5. )) continue;
       else ntuple_.fill_jetInfo(jet);

       int flavor = std::abs( jet.partonFlavour() );       
       for( const std::string &bDiscr : bDiscriminators_ )
	 {
	   ntuple_.fill_bjetInfo(jet, bDiscr, flavor);
	 }

     }
   }
   //=============================================================
   //
   //            Missing Energy 
   //     
   //=============================================================    
   pat::METCollection mets;
   if(metsHandle.isValid()){ 
     mets = *metsHandle;
     const pat::MET met = mets.front();
     ntuple_.fill_metInfo(met);
   }   
   }
   tree_->Fill();   
}
   
// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNeutralLeptonAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNeutralLeptonAnalysis::endJob() {
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

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
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

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
  
  reco::VertexCollection getMatchedVertex_Muon(const pat::Muon & mu, const reco::VertexCollection& vertexCollection);
  reco::VertexCollection getMatchedVertex_Electron(const pat::Electron & ele, const reco::VertexCollection& vtxCollection);
  reco::VertexCollection PrimaryVertex( const reco::VertexCollection &vtx);
  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle);
  //double MatchGenMuon(const edm::Event&,  reco::TrackRef BestTrack , int pdgId);
  //double MatchGenElectron(const edm::Event& iEvent, const pat::Electron& ele_, int pdgId);
  //double MatchGenVertex(float vgen_x, float vgen_y, float vgen_z, reco::Vertex vreco);
  std::pair<float, float> MatchGenVertex(float vgen_x, float vgen_y, float vgen_z, reco::Vertex vreco);
  std::pair<double, double> MatchGenLeptons(float v_eta, float v_phi, const reco::Candidate& lepton);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void initialize(const edm::Event&); 
  virtual void endJob() override;
    
  edm::Service<TFileService> fs;
  TTree * tree_;
  BigNtuple ntuple_;


  //--------------Variables------------------------
  int    debug;
  float pvCompatibilityScore = .05;

  float npT = -1;
  float npIT = -1;

  bool   isMC;
  bool   isMCSignal;

  SmartSelectionMonitor mon;

  //--------------Template-------------------------

  // ----------member data ---------------------------
  TH1D* prova;
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

  const std::vector        <std::string> bDiscbbToken_;
  const std::vector        <std::string> bDiscbbbToken_;
  const std::vector         <std::string> bDiscbcToken_;

  //edm::EDGetTokenT         <edm::ValueMap<float>>     eleMvaToken_;
  edm::EDGetTokenT         <edm::ValueMap<bool>>      eleVetoToken_;
  edm::EDGetTokenT         <edm::ValueMap<bool>>      eleLooseToken_;
  edm::EDGetTokenT         <edm::ValueMap<bool>>      eleMediumToken_;
  edm::EDGetTokenT         <edm::ValueMap<bool>>      eleTightToken_;
  
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
  
  //edm::Handle         <edm::ValueMap<float>> electronsMva;
  edm::Handle         <edm::ValueMap<bool>> electronsVeto;
  edm::Handle         <edm::ValueMap<bool>> electronsLoose;
  edm::Handle         <edm::ValueMap<bool>> electronsMedium;
  edm::Handle         <edm::ValueMap<bool>> electronsTight;

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
  isMCSignal(iConfig.getParameter<bool>("isMCSignal")),
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
  bDiscbbToken_(iConfig.getParameter<std::vector<std::string> >("bDiscbb")),
  bDiscbbbToken_(iConfig.getParameter<std::vector<std::string> >("bDiscbbb")),
  bDiscbcToken_(iConfig.getParameter<std::vector<std::string> >("bDiscbc")),
  //eleMvaToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("electronsMva"))),
  eleVetoToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronsVeto"))),
  eleLooseToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronsLoose"))),
  eleMediumToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronsMedium"))),
  eleTightToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronsTight")))
{

  //now do what ever initialization is needed
  usesResource("TFileService");

  //prova = fs->make<TH1D>("prova", "prova", 2000, -1001, 1000);

  tree_ = fs->make<TTree>("tree_", "tree");
  ntuple_.set_evtInfo(tree_);
  if(isMC && isMCSignal){
    ntuple_.set_pv_genInfo(tree_);
    ntuple_.set_sv_genInfo(tree_);
  }
  ntuple_.set_trigInfo(tree_);
  ntuple_.set_pvInfo(tree_);
  ntuple_.set_muInfo(tree_);
  ntuple_.set_eleInfo(tree_);
  ntuple_.set_eleIDInfo(tree_);
  ntuple_.set_sv_mu_Info(tree_);
  ntuple_.set_sv_ele_Info(tree_);
  ntuple_.set_jetInfo(tree_);
  ntuple_.set_metInfo(tree_);
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

  //iEvent.getByToken(eleMvaToken_, electronsMva);
  iEvent.getByToken(eleVetoToken_,electronsVeto);
  iEvent.getByToken(eleLooseToken_,electronsLoose);
  iEvent.getByToken(eleMediumToken_,electronsMedium);
  iEvent.getByToken(eleTightToken_,electronsTight);

  if (isMC){
    iEvent.getByToken(genParticleToken_, genHandle);
    iEvent.getByToken(genEventInfoToken_, genEventInfoHandle);
    iEvent.getByToken(PUInfoToken_, puInfoH);
    iEvent.getByToken(lheEventProductToken_, lheEPHandle);
  }
}



//
// member functions
//===================================== vertex matching mu ================================================//   
reco::VertexCollection HeavyNeutralLeptonAnalysis::getMatchedVertex_Muon(const pat::Muon & muon, const reco::VertexCollection& vertexCollection){
  reco::VertexCollection  matchedVertices;
  const pat::PackedCandidate* cand = dynamic_cast<const pat::PackedCandidate*>(muon.sourceCandidatePtr(0).get());
  if(!cand) {
    cout << "THIS SHOULD NEVER HAPPEN! No packed candidated associated to muon?!" << endl;
  }
  for(reco::VertexCollection::const_iterator ss = vertexCollection.begin(); ss != vertexCollection.end(); ++ss) {    
    for(reco::Vertex::trackRef_iterator tt = ss->tracks_begin(); tt != ss->tracks_end(); ++tt) {
      float   dpt    = fabs(cand->pseudoTrack().pt() - tt->castTo<reco::TrackRef>()->pt());
      if( (cand->pseudoTrack().pt() == tt->castTo<reco::TrackRef>()->pt()) || dpt < 0.001) {
        matchedVertices.push_back(*ss);
	break;
      }
    }
  } 
  return matchedVertices;
}
//===================================== vertex matching ele ================================================// 
reco::VertexCollection HeavyNeutralLeptonAnalysis::getMatchedVertex_Electron(const pat::Electron & ele, const reco::VertexCollection& vtxCollection){
  reco::VertexCollection  matchedVertices;
  for(reco::VertexCollection::const_iterator ss = vtxCollection.begin(); ss != vtxCollection.end(); ++ss) {
    for(reco::Vertex::trackRef_iterator tt = ss->tracks_begin(); tt != ss->tracks_end(); ++tt) {
      for(edm::Ref<pat::PackedCandidateCollection> cand : ele.associatedPackedPFCandidates()){
	float   dpt    = fabs(cand->pt() -  tt->castTo<reco::TrackRef>()->pt());
	float   deta   = fabs(cand->eta() - tt->castTo<reco::TrackRef>()->eta());
	if(dpt < 0.001 && deta < 0.001 && cand->charge() != 0){  //Je: to check! are those values ok? Don't we have a safer way to make the match as for the muons
	    matchedVertices.push_back(*ss);
	    break;
	  }
	  }
      }
    }
    return matchedVertices;
  }
//===================================== primary vertex selection ============================================//  
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
/*double HeavyNeutralLeptonAnalysis::MatchGenMuon(const edm::Event& iEvent, reco::TrackRef BestTrack, int pdgId) {
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
*/
//===================================== Electron - Muon Gen Match ================================================//
std::pair<double, double> HeavyNeutralLeptonAnalysis::MatchGenLeptons(float v_eta, float v_phi, const reco::Candidate& lepton){
  
  double recoGenDeltaR = 0.2;
  double rho = sqrt(lepton.vx()*lepton.vx() + lepton.vy()*lepton.vy());
  double dR = deltaR(lepton.eta(), lepton.phi(), v_eta, v_phi);

  if(rho < 0.5 and dR < recoGenDeltaR){ 
    return make_pair(rho, dR);
  }else{return make_pair(-999., -999.);}
}//returns the lepton's tipe (true if is muon - false if is electron) and 1 or 0 if the lepton matches with the vertex defined by v_eta and v_phi


//===================================== Displaced Vertex Gen Match ================================================// 
std::pair<float, float> HeavyNeutralLeptonAnalysis::MatchGenVertex(float vgen_x, float vgen_y, float vgen_z, reco::Vertex vreco) {
  double metric_xyz = std::sqrt(((vgen_x-vreco.x())*(vgen_x-vreco.x())) + ((vgen_y-vreco.y())*(vgen_y-vreco.y())) + ((vgen_z-vreco.z())*(vgen_z-vreco.z())));
  double metric_xy = std::sqrt(((vgen_x-vreco.x())*(vgen_x-vreco.x())) + ((vgen_y-vreco.y())*(vgen_y-vreco.y())));      
  return make_pair(metric_xyz, metric_xy);
}
  
  
// ------------ method called for each event  ------------
void HeavyNeutralLeptonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  
  initialize(iEvent);

  ntuple_.reset();
  ntuple_.fill_evtInfo(iEvent.id());
  //============================================================= 
  //
  //                 Gen level Info
  //
  //=============================================================
  if(isMCSignal){
    // mu && ele @ pv
    
    std::vector<reco::GenParticle> genParticles = *(genHandle.product());
    reco::GenParticle majN;

    for(auto& genPart: genParticles){
      if ((genPart.pdgId() == 9900012 or genPart.pdgId() == 9900014 or genPart.pdgId() == 9900016) and genPart.isLastCopy()){
	majN = genPart;
	ntuple_.fill_pv_genInfo(genPart,genParticles);
	break;
      }
    }
    
    set<const reco::GenParticle*> finalParticles;

    for(auto& genPart: genParticles){
      if (genPart.status()==1 and genPart.isLastCopy()){
	if (isAncestor(&majN, &genPart)){
	  const reco::GenParticle* mother = static_cast<const reco::GenParticle*>(genPart.mother(0));
	  if (mother->pdgId()==111 and (mother->isLastCopy())){
	    finalParticles.insert(static_cast<const reco::GenParticle*>(genPart.mother(0)));
	  }
	  else{
    	    finalParticles.insert(&genPart);
   	  }
    	}
      }
    }
    //set 2 vector
    vector<reco::GenParticle> final_particles;
    for(auto& gp : finalParticles) {
      final_particles.push_back(*gp);
    }
    ntuple_.fill_sv_genInfo(majN, final_particles);
  }


  //============================================================= 
  //
  //                 Primary Vertex
  // 
  //=============================================================   
  if(!vtxHandle.isValid()) return;
  reco::VertexCollection pvs = PrimaryVertex(*vtxHandle);  
  if(!pvs.size()) return;
  ntuple_.fill_pvInfo(pvs);    
  //=============================================================
  //
  //             Trigger Info
  //    
  //=============================================================     
   const edm::TriggerResults triggerResults =  *triggerResultsHandle.product();
   const edm::TriggerNames&    trigNames  = iEvent.triggerNames(triggerResults);
   ntuple_.fill_trigInfo(triggerResults, trigNames);    

   //=============================================================
   //
   //                Select Muon To Fill
   //
   //=============================================================
   std::string muon;
   pat::MuonCollection muons;
   vector<pat::Muon> looseMuons;
   if(muonsHandle.isValid()){ 
     muons = *muonsHandle;
     for(auto& mu : muons){
     if (!( fabs(mu.eta()) < 2.4 && mu.pt() > 5. )) continue;
     if (!mu.isLooseMuon()) continue;
     looseMuons.push_back(mu);
     }
   }
   // lambda function to sort this muons
   std::sort(looseMuons.begin(), looseMuons.end(), [](pat::Muon a, pat::Muon b) {return a.p() > b.p(); }); //p ordering

   //============================================================= 
   //
   //                Method for electrons
   //                                                   
   //============================================================= 

   EcalRecHitCollection recHitCollectionEB;
   if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}

   EcalRecHitCollection recHitCollectionEE;
   if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}

   vector<pat::Electron>  looseElectrons;
   // using iterator to get ref for ele
   for(auto ele = electronsHandle->begin(); ele != electronsHandle->end(); ++ele){
     if(ele->gsfTrack().isNull() || ele->pt() < 5 || fabs(ele->eta()) > 2.5 )      continue;
    
     double rho = *(rhoHandle.product());
     
     auto eleRef = edm::Ref<std::vector<pat::Electron>>(electronsHandle, (ele - electronsHandle->begin()));
     
     std::auto_ptr<EcalClusterLazyTools> recHitEcal;
     recHitEcal.reset(new EcalClusterLazyTools( iEvent, iSetup, recHitEBToken_, recHitEEToken_ ));
     
     std::pair<double, double> matching_1stele = (isMC && isMCSignal) ? MatchGenLeptons(ntuple_.get_lep1_eta(), ntuple_.get_lep1_phi(), *ele) : -999;
     std::pair<double, double> matching_2ndele = (isMC && isMCSignal) ? MatchGenLeptons(ntuple_.get_lep2_eta(), ntuple_.get_lep2_phi(), *ele) : -999;
     ntuple_.fill_eleInfo(*ele, pvs.at(0), rho , matching_1stele, matching_2ndele, recHitEcal);
     
     if(ele->full5x5_sigmaIetaIeta() <  0.036 && ele->passConversionVeto() == 1) looseElectrons.push_back(*ele); 

     pat::Electron eleMva = *ele;
     float  ele_Mva_   = eleMva.electronID("mvaEleID-Fall17-iso-V1-wp90");
       //static_cast<float *>(pat::Electron::userFloat("<ElectronMVAEstimatorRun2Fall17IsoV1>Values"));
       //pat::Electron::electronID("mvaEleID-Fall17-iso-V1-wp90");
       //static_cast<float *>(pat::Electron::userFloat("<ElectronMVAEstimatorRun2Fall17IsoV1>Values"));
     bool  ele_Veto_   = ((*electronsVeto)[eleRef]);
     bool  ele_Loose_  = ((*electronsLoose)[eleRef]);
     bool  ele_Medium_ = ((*electronsMedium)[eleRef]);
     bool  ele_Tight_  = ((*electronsTight)[eleRef]);

     ntuple_.fill_eleIDInfo(ele_Mva_, ele_Veto_, ele_Loose_ , ele_Medium_, ele_Tight_);
   }
   
   // lambda function to sort this electrons
   std::sort(looseElectrons.begin(), looseElectrons.end(), [](pat::Electron a, pat::Electron b) {return a.pt() > b.pt(); });
   //=============================================================                    
   //                     
   //                Secondary Vertex                     
   //                     
   //============================================================= 
   vector<reco::VertexCollection> sv;
   if(secondaryVertexHandle.isValid()){
     //sv due to muon
     if(looseMuons.size() > 1){
       pat::Muon muonHNL = looseMuons[1];
       reco::VertexCollection bestVertices_mu  = getMatchedVertex_Muon(muonHNL, *secondaryVertexHandle);
       sv.push_back(bestVertices_mu);
       // check if SV doesn't match with the PV
       for (const reco::Vertex& vtx_mu : bestVertices_mu){
	 float x  = vtx_mu.x(), y = vtx_mu.y(), z = vtx_mu.z();
	 float dx = x - pvs.at(0).x() , dy = y - pvs.at(0).y(), dz = z - pvs.at(0).z();	 
         float  selIVFIsPVScore = std::sqrt((dx*dx) + (dy*dy) + (dz*dz));
	 if (selIVFIsPVScore < pvCompatibilityScore) continue;

	 std::pair<float, float> matching_vtx = (isMC && isMCSignal) ? MatchGenVertex(ntuple_.get_sv_x(), ntuple_.get_sv_y(), ntuple_.get_sv_z(), vtx_mu) : make_pair(static_cast<float>(-999.), static_cast<float>(-999.));
	 ntuple_.fill_sv_Info(vtx_mu, pvs.at(0), matching_vtx, true);	 
       }
     }
     //file muon branches with events contain sv
     if(sv.size()){
       for (const pat::Muon mu : looseMuons){
	 double rho = *(rhoHandle.product());
	 reco::TrackRef bestTrack = mu.muonBestTrack();

	   std::pair<double, double> matching_1stmu = (isMC && isMCSignal) ? MatchGenLeptons(ntuple_.get_lep1_eta(), ntuple_.get_lep1_phi(), mu) : -999;
	   std::pair<double, double> matching_2ndmu = (isMC && isMCSignal) ? MatchGenLeptons(ntuple_.get_lep2_eta(), ntuple_.get_lep2_phi(), mu) : -999;
	 //added the matching
	 ntuple_.fill_muInfo(mu, pvs.at(0), rho, matching_1stmu, matching_2ndmu);

       }
     }
     
     //sv due to electron
     if(looseElectrons.size() > 1){
       pat::Electron electronHNL = looseElectrons[1];       
       reco::VertexCollection bestVertices_ele  = getMatchedVertex_Electron(electronHNL, *secondaryVertexHandle);
       for (const reco::Vertex& vtx_ele : bestVertices_ele){
	 float x  = vtx_ele.x(), y = vtx_ele.y(), z = vtx_ele.z();
	 float dx = x - pvs.at(0).x() , dy = y - pvs.at(0).y(), dz = z - pvs.at(0).z();
         float  selIVFIsPVScore = std::sqrt((dx*dx) + (dy*dy) + (dz*dz));
	 if (selIVFIsPVScore < pvCompatibilityScore) continue;
	 std::pair<float, float> matching_vtx = (isMC && isMCSignal) ? MatchGenVertex(ntuple_.get_sv_x(), ntuple_.get_sv_y(), ntuple_.get_sv_z(), vtx_mu) : make_pair(static_cast<float>(-999.), static_cast<float>(-999.));
	 ntuple_.fill_sv_Info(vtx_ele, pvs.at(0), matching_vtx, false);
       }
     }


   }
   //============================================================= 
   //
   //             Jets 
   //       
   //=============================================================
   pat::JetCollection jets;
   if(jetsHandle.isValid() && sv.size()){ 
     jets = *jetsHandle;
     for (const pat::Jet jet : jets) {
       if( jet.pt() < 0.0 ) continue;
       if (!( fabs(jet.eta()) < 3 && jet.pt() > 5. )) continue;
       else ntuple_.fill_jetInfo(jet);
       int flavor = std::abs( jet.partonFlavour() );       
       for( const std::string &bDiscr : bDiscbbToken_ )
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
   if(metsHandle.isValid() && sv.size()){ 
     mets = *metsHandle;
     const pat::MET met = mets.front();
     ntuple_.fill_metInfo(met);
   }   
   pat::TauCollection taus;
   if(tausHandle.isValid()){ taus = *tausHandle;}

   pat::PackedCandidateCollection pfCandidates;
   if (pfCandidatesHandle.isValid()) { pfCandidates = *pfCandidatesHandle; }


   tree_->Fill();   
}
   
// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNeutralLeptonAnalysis::beginJob(){
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

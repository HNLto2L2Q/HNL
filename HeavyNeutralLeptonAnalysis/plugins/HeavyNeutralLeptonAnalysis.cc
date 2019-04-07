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
#include "TVector3.h"
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

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"

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
  std::pair<double, double> MatchGenVertex(float vgen_x, float vgen_y, float vgen_z, reco::Vertex vreco);
  std::pair<double, double> MatchGenLeptons(float v_eta, float v_phi, const reco::Candidate& lepton);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  //transverse mass
  inline double MT(TLorentzVector *l, TLorentzVector *met) {
    return sqrt(pow(l->Pt() + met->Pt(), 2) - pow(l->Px() + met->Px(), 2) - pow(l->Py() + met->Py(), 2));
  }
  

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
  edm::EDGetTokenT            <std::vector<pat::Jet>> jetSmearedToken_;
  edm::EDGetTokenT            <std::vector<pat::Jet>> jetSmearedUpToken_;
  edm::EDGetTokenT            <std::vector<pat::Jet>> jetSmearedDownToken_;
  edm::EDGetTokenT             < pat::METCollection > pfMETAODToken_;
  edm::EDGetTokenT            < edm::TriggerResults > triggerResultsToken_;
  edm::EDGetTokenT            < edm::TriggerResults > metFilterResultsToken_;
  edm::EDGetTokenT    < reco::GenParticleCollection > genParticleToken_;
  edm::EDGetTokenT            < GenEventInfoProduct > genEventInfoToken_;
  edm::EDGetTokenT < std::vector<PileupSummaryInfo> > PUInfoToken_;
  edm::EDGetTokenT                < LHEEventProduct > lheEventProductToken_;
  edm::EDGetTokenT         < reco::VertexCollection > inclusiveSecondaryVertices_;

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
  edm::Handle            <std::vector<pat::Jet>> jetsSmeared; 
  edm::Handle            <std::vector<pat::Jet>> jetsSmearedUp; 
  edm::Handle            <std::vector<pat::Jet>> jetsSmearedDown; 
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

  //Prefiring stuff
  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;


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
  jetSmearedToken_( consumes<std::vector<pat::Jet>>( iConfig.getParameter<edm::InputTag>("jetsSmeared"))),
  jetSmearedUpToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetsSmearedUp"))),
  jetSmearedDownToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetsSmearedDown"))),
  pfMETAODToken_(mayConsume<pat::METCollection>(iConfig.getParameter<edm::InputTag>("pfMETSrc"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResultSrc"))),
  metFilterResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterResultSrc"))),
  genParticleToken_(mayConsume<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  genEventInfoToken_(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
  PUInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
  lheEventProductToken_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts"))),
  inclusiveSecondaryVertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("SecondaryVertices"))),
  eleVetoToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronsVeto"))),
  eleLooseToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronsLoose"))),
  eleMediumToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronsMedium"))),
  eleTightToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("electronsTight")))
  //prefweight_token(consumes< double >(iConfig.getParameter<edm::InputTag>("prefiringweight:NonPrefiringProb"))),
  //prefweightup_token(consumes< double >(iConfig.getParameter<edm::InputTag>("prefiringweight:NonPrefiringProbUp"))),
  //prefweightdown_token(consumes< double >(iConfig.getParameter<edm::InputTag>("prefiringweight:NonPrefiringProbDown")))
{

  prefweight_token = consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProb"));
  prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProbUp"));
  prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:NonPrefiringProbDown"));
  //now do what ever initialization is needed
  usesResource("TFileService");

  //prova = fs->make<TH1D>("prova", "prova", 2000, -1001, 1000);

  tree_ = fs->make<TTree>("tree_", "tree");
  ntuple_.set_evtInfo(tree_);
  if(isMC && isMCSignal){
    ntuple_.set_pv_genInfo(tree_);
    ntuple_.set_sv_genInfo(tree_);
  }
  ntuple_.set_prefiring(tree_);
  ntuple_.set_trigInfo(tree_);
  ntuple_.set_pileupInfo(tree_);
  ntuple_.set_pvInfo(tree_);
  ntuple_.set_muInfo(tree_);
  ntuple_.set_eleInfo(tree_);
  ntuple_.set_eleIDInfo(tree_);
  ntuple_.set_sv_Info(tree_);
  ntuple_.set_jetInfo(tree_);
  ntuple_.set_metInfo(tree_);
  ntuple_.set_transverseMassInfo(tree_);
  ntuple_.set_massCorrection(tree_);  
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
  iEvent.getByToken(jetSmearedToken_, jetsSmeared);
  iEvent.getByToken(jetSmearedUpToken_, jetsSmearedUp);
  iEvent.getByToken(jetSmearedDownToken_, jetsSmearedDown);
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
  if(ancestor->pdgId() == particle->pdgId() ) return true;
  for(size_t i=0;i< particle->numberOfMothers();i++)
    {
      if(isAncestor(ancestor,particle->mother(i))) return true;
    }
  return false;
}

//===================================== Electron - Muon Gen Match ================================================//
std::pair<double, double> HeavyNeutralLeptonAnalysis::MatchGenLeptons(float v_eta, float v_phi, const reco::Candidate& lepton){
  
  double recoGenDeltaR = 0.2;
  double rho = sqrt(lepton.vx()*lepton.vx() + lepton.vy()*lepton.vy());
  double dR = deltaR(lepton.eta(), lepton.phi(), v_eta, v_phi);
  
  if(rho < 100 and dR < recoGenDeltaR){ 
    return make_pair(rho, dR);
  }else{return make_pair(-999., -999.);}
}//returns the lepton's tipe (true if is muon - false if is electron) and 1 or 0 if the lepton matches with the vertex defined by v_eta and v_phi


//===================================== Displaced Vertex Gen Match ================================================// 
std::pair<double, double> HeavyNeutralLeptonAnalysis::MatchGenVertex(float vgen_x, float vgen_y, float vgen_z, reco::Vertex vreco) {
  double metric_xyz = std::sqrt(((vgen_x-vreco.x())*(vgen_x-vreco.x())) + ((vgen_y-vreco.y())*(vgen_y-vreco.y())) + ((vgen_z-vreco.z())*(vgen_z-vreco.z())));
  double metric_xy = std::sqrt(((vgen_x-vreco.x())*(vgen_x-vreco.x())) + ((vgen_y-vreco.y())*(vgen_y-vreco.y())));      
  if(metric_xyz < 0.2){
  return make_pair(metric_xyz, metric_xy);
  }else{return make_pair(-999., -999.);}
}
  
  
// ------------ method called for each event  ------------
void HeavyNeutralLeptonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  
  initialize(iEvent);

  ntuple_.reset();
  ntuple_.fill_evtInfo(iEvent.id());

  edm::Handle< double > theprefweight;
  iEvent.getByToken(prefweight_token, theprefweight ) ;
  float _prefiringweight =(*theprefweight);

  edm::Handle< double > theprefweightup;
  iEvent.getByToken(prefweightup_token, theprefweightup ) ;
  float _prefiringweightup =(*theprefweightup);

  edm::Handle< double > theprefweightdown;
  iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
  float _prefiringweightdown =(*theprefweightdown);
  ntuple_.fill_prefiring(_prefiringweight, _prefiringweightup, _prefiringweightdown);

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
      if ((genPart.pdgId() == 9990012 or genPart.pdgId() == 9900012 or genPart.pdgId() == 9900014 or genPart.pdgId() == 9900016) and  abs(genPart.mother()->pdgId()) == 24){	
	majN = genPart;
	  ntuple_.fill_pv_genInfo(genPart,genParticles);
	break;
      }
    }
    
    /*
    set<const reco::GenParticle*> finalParticles;
    // vector<reco::GenParticle> finalParticles;    
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
    */
    //set 2 vector
    vector<reco::GenParticle> final_particles;
    for(auto& gp: genParticles) {
      if (gp.status()==1 and gp.isLastCopy()){
	const reco::Candidate * mom = gp.mother(0) ;
	if (isAncestor(&majN, mom)  ){
	  final_particles.push_back(gp);
	  ntuple_.fill_sv_genInfo(majN, final_particles);
	}
      }
    }    
    
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
   //               Pile Up Info
   //
   //=============================================================
   if(isMC){     
     std::vector<PileupSummaryInfo>::const_iterator PVI;
     for(PVI = puInfoH->begin(); PVI != puInfoH->end(); ++PVI) {
       int BX = PVI->getBunchCrossing();
       if(BX == 0) {
	 npT = PVI->getTrueNumInteractions();
	 npIT = PVI->getPU_NumInteractions();
       }
     }
     ntuple_.fill_pileupInfo(npT, npIT );
   }    
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
   
   //fill muon branches with events at least 2 loose muons                                                                                                                                                
   if(looseMuons.size() > 1){
     for (const pat::Muon mu : looseMuons){
       double rho = *(rhoHandle.product());
       reco::TrackRef bestTrack = mu.muonBestTrack();
       std::pair<double, double> matching_1stmu = (isMC && isMCSignal) ? MatchGenLeptons(ntuple_.get_lep1_eta(), ntuple_.get_lep1_phi(), mu) : make_pair(-999., -999.);
       std::pair<double, double> matching_2ndmu = (isMC && isMCSignal) ? MatchGenLeptons(ntuple_.get_lep2_eta(), ntuple_.get_lep2_phi(), mu) : make_pair(-999., -999.);
       //added the matching                                                                                                                                                                               
       ntuple_.fill_muInfo(mu, pvs.at(0), rho, matching_1stmu, matching_2ndmu);
     }
   }

   //============================================================= 
   //
   //                Method for electrons
   //                                                   
   //============================================================= 

   EcalRecHitCollection recHitCollectionEB;
   if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}

   EcalRecHitCollection recHitCollectionEE;
   if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}

   //vector<pat::Electron>  looseElectrons;
   std::vector<pat::Electron> looseElectrons;

   // using iterator to get ref for ele
   for(auto ele = electronsHandle->begin(); ele != electronsHandle->end(); ++ele){
     if(ele->gsfTrack().isNull() || ele->pt() < 5 || fabs(ele->eta()) > 2.5 )      continue;
     if(ele->full5x5_sigmaIetaIeta() <  0.036 && ele->passConversionVeto() == 1) looseElectrons.push_back(*ele);
   }

   if(looseElectrons.size() > 1){
     for(auto ele = electronsHandle->begin(); ele != electronsHandle->end(); ++ele){
       if(ele->gsfTrack().isNull() || ele->pt() < 5 || fabs(ele->eta()) > 2.5 )      continue;
       
       double rho = *(rhoHandle.product());
       
       auto eleRef = edm::Ref<std::vector<pat::Electron>>(electronsHandle, (ele - electronsHandle->begin()));
       
       std::auto_ptr<EcalClusterLazyTools> recHitEcal;
       recHitEcal.reset(new EcalClusterLazyTools( iEvent, iSetup, recHitEBToken_, recHitEEToken_ ));
       
       std::pair<double, double> matching_1stele = (isMC && isMCSignal) ? MatchGenLeptons(ntuple_.get_lep1_eta(), ntuple_.get_lep1_phi(), *ele) : make_pair(-999., -999.);
       std::pair<double, double> matching_2ndele = (isMC && isMCSignal) ? MatchGenLeptons(ntuple_.get_lep2_eta(), ntuple_.get_lep2_phi(), *ele) : make_pair(-999., -999.);
       ntuple_.fill_eleInfo(*ele, pvs.at(0), rho , matching_1stele, matching_2ndele, recHitEcal);
       
       pat::Electron eleMva = *ele;
       float  ele_Mva_   = eleMva.electronID("mvaEleID-Spring16-GeneralPurpose-V1-wp80");
       //float  ele_Mva_   = eleMva.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90"); to try with 2016
       bool  ele_Veto_   = ((*electronsVeto)[eleRef]);
       bool  ele_Loose_  = ((*electronsLoose)[eleRef]);
       bool  ele_Medium_ = ((*electronsMedium)[eleRef]);
       bool  ele_Tight_  = ((*electronsTight)[eleRef]);
       
       ntuple_.fill_eleIDInfo(ele_Mva_, ele_Veto_, ele_Loose_ , ele_Medium_, ele_Tight_);
     }
   }
   
   // lambda function to sort this electrons
   std::sort(looseElectrons.begin(), looseElectrons.end(), [](pat::Electron a, pat::Electron b) {return a.pt() > b.pt(); });
   //=============================================================                    
   //                     
   //                Secondary Vertex                     
   //                     
   //============================================================= 
   if(secondaryVertexHandle.isValid()){
     //sv due to muon
     if(looseMuons.size() > 1){
       pat::Muon muonHNL = looseMuons[1];
       reco::VertexCollection bestVertices_mu  = getMatchedVertex_Muon(muonHNL, *secondaryVertexHandle);
       // check if SV doesn't match with the PV
       for (const reco::Vertex& vtx_mu : bestVertices_mu){
	 float x  = vtx_mu.x(), y = vtx_mu.y(), z = vtx_mu.z();
	 float dx = x - pvs.at(0).x() , dy = y - pvs.at(0).y(), dz = z - pvs.at(0).z();	 
         float  selIVFIsPVScore = std::sqrt((dx*dx) + (dy*dy) + (dz*dz));
	 if (selIVFIsPVScore < pvCompatibilityScore) continue;
	 std::pair<double, double> matching_vtx = (isMC && isMCSignal) ? MatchGenVertex(ntuple_.get_sv_x(), ntuple_.get_sv_y(), ntuple_.get_sv_z(), vtx_mu) : make_pair(-999.,-999.);
	 ntuple_.fill_sv_Info(vtx_mu, pvs.at(0), matching_vtx, true);	 
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
	 //std::pair<float, float> matching_vtx = (isMC && isMCSignal) ? MatchGenVertex(ntuple_.get_sv_x(), ntuple_.get_sv_y(), ntuple_.get_sv_z(), vtx_ele) : make_pair(static_cast<float>(-999.), static_cast<float>(-999.));
	 std::pair<double, double> matching_vtx = (isMC && isMCSignal) ? MatchGenVertex(ntuple_.get_sv_x(), ntuple_.get_sv_y(), ntuple_.get_sv_z(), vtx_ele) : make_pair(-999.,-999.);
	 ntuple_.fill_sv_Info(vtx_ele, pvs.at(0), matching_vtx, false);
       }
     }


   }
   //============================================================= 
   //
   //             Jets 
   //       
   //=============================================================

   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
   iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
   JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
   JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);


   pat::JetCollection jets;
   //std::vector<tuple < std::vector<std::string>, std::vector<std::string>, std::vector<std::string> > > bDisc_info = make_tuple(bDiscbbToken_, bDiscbbbToken_, bDiscbcToken_); 
   if(jetsHandle.isValid() && ( looseMuons.size() > 1 || looseElectrons.size() > 1 )  ){ 
     jets = *jetsHandle;
     for (const pat::Jet jet : jets) {
       if (!( fabs(jet.eta()) < 3 && jet.pt() > 5. )) continue;

       auto jetSmearedIt = jetsSmeared->begin();
       for(auto j = jetsSmeared->cbegin(); j != jetsSmeared->cend(); ++j){
	 if(reco::deltaR(jet, *j) < reco::deltaR(jet, *jetSmearedIt)) jetSmearedIt = j;
       }
       
       auto jetSmearedUpIt = jetsSmearedUp->begin();
       for(auto j = jetsSmearedUp->cbegin(); j != jetsSmearedUp->cend(); ++j){
	 if(reco::deltaR(jet, *j) < reco::deltaR(jet, *jetSmearedUpIt))  jetSmearedUpIt = j;
       }

       auto jetSmearedDownIt = jetsSmearedDown->begin();
       for(auto j = jetsSmearedDown->cbegin(); j != jetsSmearedDown->cend(); ++j){
	 if(reco::deltaR(jet, *j) < reco::deltaR(jet, *jetSmearedDownIt))  jetSmearedDownIt = j;
       }

       jecUnc->setJetEta(jet.eta());
       jecUnc->setJetPt(jet.pt());
       double unc = jecUnc->getUncertainty(true);

       //double unc = 9 ;

       jecUnc->setJetEta(jetSmearedIt->eta());
       jecUnc->setJetPt(jetSmearedIt->pt());
       double uncSmeared = jecUnc->getUncertainty(true);

       //double uncSmeared = 9;

       float jetSmearedPt         = jetSmearedIt->pt();
       float jetSmearedPt_JERDown = jetSmearedDownIt->pt();
       float jetSmearedPt_JERUp   = jetSmearedUpIt->pt();

       ntuple_.fill_jetInfo(jet ,jetSmearedPt ,jetSmearedPt_JERUp ,jetSmearedPt_JERDown ,unc ,uncSmeared );

     }
   }
   //=============================================================
   //
   //            Missing Energy 
   //     
   //=============================================================    
   pat::METCollection mets;

   if(metsHandle.isValid() && ( looseMuons.size() > 1 || looseElectrons.size() > 1 ) ){ 
     mets = *metsHandle;
     const pat::MET met = mets.front();
     ntuple_.fill_metInfo(met);
   }

   pat::TauCollection taus;
   if(tausHandle.isValid()){ taus = *tausHandle;}
   
   pat::PackedCandidateCollection pfCandidates;
   if (pfCandidatesHandle.isValid()) { pfCandidates = *pfCandidatesHandle; }

   TLorentzVector ivf_, prompt_lep;
   ivf_.SetPxPyPzE(ntuple_.get_sv_px(), ntuple_.get_sv_py(), ntuple_.get_sv_pz(), ntuple_.get_sv_en());
   prompt_lep.SetPtEtaPhiE(ntuple_.get_lep1_pt(), ntuple_.get_lep1_eta(), ntuple_.get_lep1_phi(), ntuple_.get_lep1_en());
   TLorentzVector ivfPlus_lep1 = ivf_ + prompt_lep;

   if(ntuple_.get_met_px() > 0){ 

   float transverse_mass_ivf = sqrt(pow(ivf_.Pt() + ntuple_.get_met_pt(), 2) - pow(ivf_.Px() + ntuple_.get_met_px(), 2) - pow(ivf_.Py() + ntuple_.get_met_py(), 2));//ivf transverse mass
   float transverse_mass_lep1 = sqrt(pow(prompt_lep.Pt() + ntuple_.get_met_pt(), 2) - pow(prompt_lep.Px() + ntuple_.get_met_px(), 2) - pow(prompt_lep.Py() + ntuple_.get_met_py(), 2));//prompt lepton transverse mass
   float transverse_mass_ivfPlus_lep1 = sqrt(pow(ivfPlus_lep1.Pt() + ntuple_.get_met_pt(), 2) - pow(ivfPlus_lep1.Px() + ntuple_.get_met_px(), 2) - pow(ivfPlus_lep1.Py() + ntuple_.get_met_py(), 2));//ivf + prompt lepton transverse mass

   ntuple_.fill_transverseMassInfo(transverse_mass_ivf, transverse_mass_lep1, transverse_mass_ivfPlus_lep1); //ivf transverse mass

   //mass correction
   TVector3 dir_sv, ivf3D_;
   dir_sv.SetMagThetaPhi(ntuple_.get_pvTosv_rho(), ntuple_.get_pvTosv_phi(), ntuple_.get_pvTosv_theta());
   ivf3D_.SetXYZ(ntuple_.get_sv_recox(), ntuple_.get_sv_recoy(), ntuple_.get_sv_recoz());
   
   double ivf_mass = ivf_.M();
   double vertexPt2 = dir_sv.Cross(ivf3D_).Mag2() / dir_sv.Mag2();
   double mass_corrected = std::sqrt(ivf_mass * ivf_mass + vertexPt2) + std::sqrt(vertexPt2);
   
   ntuple_.fill_massCorrection(mass_corrected);
   }
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

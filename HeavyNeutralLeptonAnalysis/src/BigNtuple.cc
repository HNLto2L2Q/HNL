#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "TVector3.h"
using namespace std;

void BigNtuple::set_evtInfo(TTree* tree) {
  tree->Branch("run" , &run_, "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt" , &evt_, "evt/i");
}

void BigNtuple::fill_evtInfo(const edm::EventID& id) {
  lumi_ = id.run();
  run_  = id.luminosityBlock();
  evt_  = id.event();
}


void BigNtuple::set_prefiring(TTree* tree){
  tree->Branch("prefire_weight" , &prefire_weight_);
  tree->Branch("prefire_weightup" , &prefire_weightup_);
  tree->Branch("prefire_weightdown" , &prefire_weightdown_);
}

void BigNtuple::fill_prefiring(float weight, float weightup, float weightdown){
  prefire_weight_.push_back(weight);
  prefire_weightup_.push_back(weightup);
  prefire_weightdown_.push_back(weightdown);
}

void BigNtuple::set_mc_genInfo(TTree* tree) {
  tree->Branch("lep_PID",   &lep_PID_);
  tree->Branch("lep_Charge",&lep_Charge_);
  tree->Branch("lep_Pt",    &lep_Pt_);
  tree->Branch("lep_Eta",   &lep_Eta_);
  tree->Branch("lep_Phi",   &lep_Phi_);
  tree->Branch("lep_En",    &lep_En_);
  tree->Branch("lep_motherID",  &lep_motherID_);
  tree->Branch("lep_DecayChain",&lep_DecayChain_);
  tree->Branch("lep_VX",    &lep_VX_);
  tree->Branch("lep_VY",    &lep_VY_);
  tree->Branch("lep_VZ",    &lep_VZ_);
  tree->Branch("lep_status",&lep_status_);
}
void  BigNtuple::fill_mc_genInfo(const reco::Candidate* genPart , float decay){

    lep_PID_.push_back(genPart->pdgId());
    lep_Charge_.push_back(genPart->charge());
    lep_Pt_.push_back(genPart->pt());
    lep_Eta_.push_back(genPart->eta());
    lep_Phi_.push_back(genPart->phi());
    lep_En_.push_back(genPart->energy());
    lep_motherID_.push_back(genPart->mother()->pdgId());
    lep_DecayChain_.push_back(decay);
    lep_VX_.push_back(genPart->vertex().x());
    lep_VY_.push_back(genPart->vertex().y());
    lep_VZ_.push_back(genPart->vertex().z());
    lep_status_.push_back(genPart->status());
}

void BigNtuple::set_pv_genInfo(TTree* tree) {
  
  tree->Branch("lep1_gen_PID" , &lep1_gen_PID_, "lep1_gen_PID/I");
  tree->Branch("lep1_gen_Charge",&lep1_gen_Charge_, "lep1_gen_Charge/I");
  tree->Branch("lep1_gen_Pt",&lep1_gen_Pt_,"lep1_gen_Pt/F");
  tree->Branch("lep1_gen_Eta",&lep1_gen_Eta_,"lep1_gen_Eta/F");
  tree->Branch("lep1_gen_Phi",&lep1_gen_Phi_,"lep1_gen_Phi/F");
  tree->Branch("lep1_gen_En",&lep1_gen_En_,"lep1_gen_En/F");
  tree->Branch("lep1_gen_vx",&lep1_gen_vx_,"lep1_gen_vx/F");
  tree->Branch("lep1_gen_vy",&lep1_gen_vy_,"lep1_gen_vy/F");
  tree->Branch("lep1_gen_vz",&lep1_gen_vz_,"lep1_gen_vz/F");
  tree->Branch("lep1_gen_Lxy",&lep1_gen_Lxy_,"lep1_gen_Lxy/F");
  tree->Branch("lep1_gen_Lxyz",&lep1_gen_Lxyz_,"lep1_gen_Lxyz/F");
  tree->Branch("HNL_gen_PID",&HNL_gen_PID_,"HNL_gen_PID/I");
  tree->Branch("HNL_gen_Mass",&HNL_gen_Mass_,"HNL_gen_Mass/F");
  tree->Branch("HNL_gen_Charge",&HNL_gen_Charge_,"HNL_gen_Charge/I");
  tree->Branch("HNL_gen_Pt",&HNL_gen_Pt_,"HNL_gen_Pt/F");
  tree->Branch("HNL_gen_Eta",&HNL_gen_Eta_,"HNL_gen_Eta/F");
  tree->Branch("HNL_gen_Phi",&HNL_gen_Phi_,"HNL_gen_Phi/F");

}


void  BigNtuple::fill_pv_genInfo(const reco::GenParticle prt ,std::vector<reco::GenParticle> prtCollection ){
  float vx = prt.vx(), vy = prt.vy(), vz = prt.vz();  
  reco::GenParticle lep1;
  auto pOrdering = [&](reco::GenParticle first, reco::GenParticle second){return first.pt() > second.pt();};
  std::sort(prtCollection.begin(), prtCollection.end(), pOrdering);
  for (auto genPart: prtCollection){
    if( abs(genPart.pdgId()) == 13 ) {
      lep1 = genPart;
      break;
    }
  }//loop on daughter of HNL to find the lept2 with highest p                                                                                                                                              
  //cout<<"the pdgId of first muon  = "<<lep1.pdgId()<<" the pT of first lep 1 = "<<lep1.pt()<<" the mother of first muon = "<<lep1.mother()->pdgId() <<endl;
  lep1_gen_PID_ = lep1.pdgId();
  // genkinematics
  lep1_gen_Pt_     = lep1.pt();
  lep1_gen_Eta_    = lep1.eta();
  lep1_gen_Phi_    = lep1.phi();
  lep1_gen_En_     = lep1.energy();
  lep1_gen_Charge_ = lep1.charge();
  // vertexposition
  lep1_gen_vx_ = vx;
  lep1_gen_vy_ = vy;
  lep1_gen_vz_ = vz;
  // vertexflight
  lep1_gen_Lxy_  = std::sqrt( vx * vx + vy * vy );
  lep1_gen_Lxyz_ = std::sqrt( vx * vx + vy * vy + vz * vz);
  // mother beta, gamma,ctau
  // HNL infos
  // gen id and status
  HNL_gen_PID_    = prt.pdgId();
  HNL_gen_Charge_ = prt.charge();
  HNL_gen_Mass_   = prt.mass();
  HNL_gen_Pt_     = prt.pt();
  HNL_gen_Eta_    = prt.eta();
  HNL_gen_Phi_    = prt.phi();
}//filling the information of the HNL (mjN) particle

void BigNtuple::set_sv_genInfo(TTree* tree) {
  
  tree->Branch("lep2_gen_PID" , &lep2_gen_PID_, "lep2_gen_PID/I");
  tree->Branch("lep2_gen_Charge",&lep2_gen_Charge_,"lep2_gen_Charge/I");
  tree->Branch("lep2_gen_Pt",&lep2_gen_Pt_,"lep2_gen_Pt/F");
  tree->Branch("lep2_gen_Eta",&lep2_gen_Eta_,"lep2_gen_Eta/F");
  tree->Branch("lep2_gen_Phi",&lep2_gen_Phi_,"lep2_gen_Phi/F");
  tree->Branch("lep2_gen_vx",&lep2_gen_vx_,"lep2_gen_vx/F");
  tree->Branch("lep2_gen_vy",&lep2_gen_vy_,"lep2_gen_vy/F");
  tree->Branch("lep2_gen_vz",&lep2_gen_vz_,"lep2_gen_vz/F");
  tree->Branch("lep2_gen_MomLxyz",&lep2_gen_MomLxyz_,"lep2_gen_MomLxyz/F");
  tree->Branch("lep2_gen_MomLz",&lep2_gen_MomLz_,"lep2_gen_MomLz/F");
  tree->Branch("lep2_gen_MomLxy",&lep2_gen_MomLxy_,"lep2_gen_MomLxy/F");
  tree->Branch("lep2_gen_MomCTau0" ,&lep2_gen_MomCTau0_,"lep2_gen_MomCTau0/F");

  tree->Branch("daugh_gen_PID"   , &daugh_gen_PID_);
  tree->Branch("daugh_gen_Charge", &daugh_gen_Charge_);
  tree->Branch("daugh_gen_Pt"    , &daugh_gen_Pt_);
  tree->Branch("daugh_gen_Eta"   , &daugh_gen_Eta_);
  tree->Branch("daugh_gen_Phi"   , &daugh_gen_Phi_);
  tree->Branch("daugh_gen_Mass"  , &daugh_gen_Mass_);

}

void  BigNtuple::fill_sv_genInfo(const reco::GenParticle hnl , std::vector<reco::GenParticle> prtCollection ){ 
  reco::GenParticle lep2;
  float vx = hnl.vx(), vy = hnl.vy(), vz = hnl.vz();
  auto pOrdering = [&](reco::GenParticle first, reco::GenParticle second){return first.pt() > second.pt();};
  std::sort(prtCollection.begin(), prtCollection.end(), pOrdering);
  for (auto genPart: prtCollection){
    if( abs(genPart.pdgId()) == 13 ) {
      lep2 = genPart;
      break;
    }
  }//loop on daughter of HNL to find the lept2 with highest p
  lep2_gen_PID_ = lep2.pdgId();
  // genkinematics
  lep2_gen_Pt_ = lep2.pt();
  lep2_gen_Eta_ = lep2.eta();
  lep2_gen_Phi_ = lep2.phi();
  lep2_gen_Charge_ = lep2.charge();
  // vertexposition
  lep2_gen_vx_ = lep2.vx();
  lep2_gen_vy_ = lep2.vy();
  lep2_gen_vz_ = lep2.vz();
  //set_sv_x(lep2.vx());
  //set_sv_y(lep2.vy());
  //set_sv_z(lep2.vz());
  
  // Vertexflight ()
  float dx = vx - lep2.vx(), dy = vy - lep2.vy(), dz = vz - lep2.vz();
  float beta_hnl  = hnl.p() / hnl.energy();
  float gamma_hnl = hnl.energy() / hnl.mass();
  
  lep2_gen_MomLxyz_ = std::sqrt(dx*dx + dy*dy + dz*dz);
  lep2_gen_MomLz_ = dz;
  lep2_gen_MomLxy_ = std::sqrt(dx*dx + dy*dy);
  lep2_gen_MomCTau0_ = std::sqrt(dx*dx + dy*dy + dz*dz) / (beta_hnl * gamma_hnl);
  
  for(auto genPart: prtCollection){
    daugh_gen_PID_.push_back(genPart.pdgId());
    daugh_gen_Pt_.push_back(genPart.pt());
    daugh_gen_Eta_.push_back(genPart.eta());
    daugh_gen_Phi_.push_back(genPart.phi());
    daugh_gen_Charge_.push_back(genPart.charge());
    daugh_gen_Mass_.push_back(genPart.mass());
  }//loop on dougheter
}//store info of HNL and all daughters


void BigNtuple::set_pvInfo(TTree* tree){
       tree->Branch("pvX" , &pvX_, "pvX/F");
       tree->Branch("pvY" , &pvY_, "pvY/F");
       tree->Branch("pvZ" , &pvZ_, "pvZ/F");
       tree->Branch("pvXErr" , &pvXErr_, "pvXErr/F");
       tree->Branch("pvYErr" , &pvYErr_, "pvYErr/F");
       tree->Branch("pvZErr" , &pvZErr_, "pvZErr/F");
       tree->Branch("pvMass" , &pvMass_, "pvMass/F");
       tree->Branch("pvLxy" , &pvLxy_, "pvLxy/F");
       tree->Branch("pvLxyz" , &pvLxyz_, "pvLxyz/F");
       tree->Branch("pvLxySigma" , &pvLxySigma_, "pvLxySigma/F");
       tree->Branch("pvLxyzSigma" , &pvLxyzSigma_, "pvLxyzSigma/F");
       tree->Branch("pvChi2" , &pvChi2_, "pvChi2/F");
       tree->Branch("pvNTrack" , &pvNTrack_, "pvNTrack/i");
       tree->Branch("pvSumPtSq" , &pvSumPtSq_, "pvSumPtSq/F");
       tree->Branch("numberPV" , &numberPV_, "numberPV/i");

}

void BigNtuple::fill_pvInfo(const reco::VertexCollection& pvs){
  numberPV_ = pvs.size();
  const reco::Vertex& pv = pvs.front();

  float x  = pv.x(), y = pv.y(), z = pv.z();
  float xE = pv.xError(), yE = pv.yError(), zE = pv.zError();
  
  pvX_ = x;
  pvY_ = y;
  pvZ_ = z;
  pvXErr_ = xE;
  pvYErr_ = yE;
  pvZErr_ = zE;
  pvMass_ = 0;
  pvLxy_ = std::sqrt( x * x + y * y );
  pvLxyz_ = std::sqrt( x * x + y * y + z * z);
  pvLxySigma_ = std::sqrt( x * x + y * y ) / std::sqrt(xE * xE + yE * yE);
  pvLxyzSigma_ = std::sqrt( x * x + y * y + z * z )/ std::sqrt(xE * xE + yE * yE + zE * zE);
  pvChi2_ = pv.chi2();
  
  reco::Vertex::trackRef_iterator vtxIter = pv.tracks_begin();
  float  SumPtSq =  0;
  int NTrack = 0;
  for(; vtxIter != pv.tracks_end(); ++vtxIter) {
    NTrack++;
    SumPtSq += (*vtxIter)->pt() * (*vtxIter)->pt();
  }
  pvNTrack_ = NTrack;
  pvSumPtSq_ = SumPtSq;
}//fill primary vertex reco info


void BigNtuple::set_trigInfo(TTree* tree){

    tree->Branch("passMu3_PFJet40"   , &passMu3_PFJet40_   , "passMu3_PFJet40/O");
    tree->Branch("passMu8_TrkIsoVVL" , &passMu8_TrkIsoVVL_ , "passMu8_TrkIsoVVL/O");
    tree->Branch("passMu17_TrkIsoVVL", &passMu17_TrkIsoVVL_, "passMu17_TrkIsoVVL/O");
    tree->Branch("passIsoMuTk18" , &passIsoMuTk18_, "passIsoMuTk18/O");
    tree->Branch("passIsoMuTk20" , &passIsoMuTk20_, "passIsoMuTk20/O");
    tree->Branch("passIsoMuTk22" , &passIsoMuTk22_, "passIsoMuTk22/O");
    tree->Branch("passIsoMuTk24" , &passIsoMuTk24_, "passIsoMuTk24/O");
    tree->Branch("passIsoMuTk27" , &passIsoMuTk27_, "passIsoMuTk27/O");
    tree->Branch("passIsoMuTk17e" , &passIsoMuTk17e_, "passIsoMuTk17e/O");
    tree->Branch("passIsoMuTk22e" , &passIsoMuTk22e_, "passIsoMuTk22e/O");
    tree->Branch("passIsoMu18" , &passIsoMu18_, "passIsoMu18/O");
    tree->Branch("passIsoMu20" , &passIsoMu20_, "passIsoMu20/O");
    tree->Branch("passIsoMu22" , &passIsoMu22_, "passIsoMu22/O");
    tree->Branch("passIsoMu24" , &passIsoMu24_, "passIsoMu24/O");
    tree->Branch("passIsoMu27" , &passIsoMu27_, "passIsoMu27/O");
    tree->Branch("passIsoMu17e" , &passIsoMu17e_, "passIsoMu17e/O");
    tree->Branch("passIsoMu22e" , &passIsoMu22e_, "passIsoMu22e/O");
    tree->Branch("passTkMu17" , &passTkMu17_, "passTkMu17/O");
    tree->Branch("passTkMu20" , &passTkMu20_, "passTkMu20/O");
    tree->Branch("passIsoMu24All" , &passIsoMu24All_, "passIsoMu24All/O");
    tree->Branch("passIsoMu27All" , &passIsoMu27All_, "passIsoMu27All/O");
    tree->Branch("passDoubleMu17TrkIsoMu8" , &passDoubleMu17TrkIsoMu8_, "passDoubleMu17TrkIsoMu8/O");
    tree->Branch("passDoubleMu17TrkIsoTkMu8" , &passDoubleMu17TrkIsoTkMu8_, "passDoubleMu17TrkIsoTkMu8/O");
    tree->Branch("passDoubleTkMu17TrkIsoTkMu8" , &passDoubleTkMu17TrkIsoTkMu8_, "passDoubleTkMu17TrkIsoTkMu8/O");
    tree->Branch("passIsoEle27", &passIsoEle27_,"passIsoEle27/O");
    tree->Branch("passNonIsoEle115", &passNonIsoEle115_,"passNonIsoEle115/O");
    tree->Branch("passDoubleEle23andEle12DZ", &passDoubleEle23andEle12DZ_,"passDoubleEle23andEle12DZ/O");
    tree->Branch("passDoubleEle23andEle12", &passDoubleEle23andEle12_,"passDoubleEle23andEle12/O");
    tree->Branch("passDoubleEle33TrkMW", &passDoubleEle33TrkMW_,"passDoubleEle33TrkMW/O");
    tree->Branch("passDoubleEle33MW", &passDoubleEle33MW_,"passDoubleEle33MW/O");
    tree->Branch("passDoubleEle33", &passDoubleEle33_,"passDoubleEle33/O");
    tree->Branch("passDoubleMu33Ele33", &passDoubleMu33Ele33_,"passDoubleMu33Ele33/O");
    tree->Branch("passEle32_WPTight_Gsf", &passEle32_WPTight_Gsf_, "passEle32_WPTight_Gsf/O");

    tree->Branch("pass_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",  & pass_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_ , "pass_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30/O");
    tree->Branch("pass_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", & pass_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_ , "pass_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30/O");
    tree->Branch("pass_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30", & pass_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_ , "pass_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30/O");
    tree->Branch("pass_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", & pass_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_ , "pass_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30/O");

}

void BigNtuple::fill_trigInfo(const edm::TriggerResults& triggerResults, const edm::TriggerNames& trigNames){

  for (size_t i = 0; i < trigNames.size(); ++i) {
    const std::string &name = trigNames.triggerName(i);
    bool fired = triggerResults.accept(i);
    if(!fired) continue;
    
    passEle32_WPTight_Gsf_ |= name.find("HLT_Ele32_WPTight_Gsf_v") != std::string::npos;
    passMu3_PFJet40_    |= name.find("HLT_Mu3_PFJet40_v")    != std::string::npos;
    passMu8_TrkIsoVVL_  |= name.find("HLT_Mu8_TrkIsoVVL_v")  != std::string::npos;
    passMu17_TrkIsoVVL_ |= name.find("HLT_Mu17_TrkIsoVVL_v") != std::string::npos;
    
    passIsoMuTk18_  |=  name.find("HLT_IsoTkMu18_v") != std::string::npos;
    passIsoMuTk20_  |=  name.find("HLT_IsoTkMu20_v") != std::string::npos;
    passIsoMuTk22_  |=  name.find("HLT_IsoTkMu22_v") != std::string::npos;
    passIsoMuTk24_  |=  name.find("HLT_IsoTkMu24_v") != std::string::npos;
    passIsoMuTk27_  |=  name.find("HLT_IsoTkMu27_v") != std::string::npos;
    passIsoMuTk17e_  |=  name.find("HLT_IsoTkMu17_eta2p1_v") != std::string::npos;
    passIsoMuTk22e_  |=  name.find("HLT_IsoTkMu22_eta2p1_v") != std::string::npos;
    passIsoMu18_  |=  name.find("HLT_IsoMu18_v") != std::string::npos;
    passIsoMu20_  |=  name.find("HLT_IsoMu20_v") != std::string::npos;
    passIsoMu22_  |=  name.find("HLT_IsoMu22_v") != std::string::npos;
    passIsoMu24_  |=  name.find("HLT_IsoMu24_v") != std::string::npos;
    passIsoMu27_  |=  name.find("HLT_IsoMu27_v") != std::string::npos;

    passIsoMu17e_ |=  name.find("HLT_IsoTkMu17_eta2p1_v") != std::string::npos;
    passIsoMu22e_ |=  name.find("HLT_IsoTkMu22_eta2p1_v") != std::string::npos;

    passTkMu17_   |=  name.find("HLT_TkMu17_v") != std::string::npos;
    passTkMu20_   |=  name.find("HLT_TkMu20_v") != std::string::npos;

    passDoubleMu17TrkIsoMu8_     |=  name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") != std::string::npos;
    passDoubleMu17TrkIsoTkMu8_   |=  name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos;
    passDoubleTkMu17TrkIsoTkMu8_ |=  name.find("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos;

    passIsoEle27_                |=  name.find("HLT_Ele27_WPTight_Gsf_v") != std::string::npos;
    passNonIsoEle115_            |=  name.find("HLT_Ele115_CaloIdVT_GsfTrkIdT_v") != std::string::npos;

    passDoubleEle23andEle12DZ_   |=  name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos;
    passDoubleEle23andEle12_     |=  name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos;

    passDoubleEle33TrkMW_        |=  name.find("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v") != std::string::npos;
    passDoubleEle33MW_           |=  name.find("HLT_DoubleEle33_CaloIdL_MW_v") != std::string::npos;
    passDoubleEle33_             |=  name.find("HLT_DoubleEle33_CaloIdL_v") != std::string::npos;

    passDoubleMu33Ele33_         |=  name.find("HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v") != std::string::npos;

    pass_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_   |=  name.find("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") != std::string::npos;
    pass_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_  |=  name.find("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") != std::string::npos;
    pass_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_  |=  name.find("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v") != std::string::npos;
    pass_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_  |=  name.find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") != std::string::npos;


  }

  passIsoMu24All_ = passIsoMu24All_   || passIsoMu24_ || passIsoMuTk24_ ;
  passIsoMu27All_ = passIsoMu27All_   || passIsoMu27_ || passIsoMuTk27_ ;

}

void BigNtuple::set_pileupInfo(TTree* tree){
  tree->Branch("npT" , &npT_);
  tree->Branch("npIT" , &npIT_);
}

void BigNtuple::fill_pileupInfo( float npt, float npit){
  npT_.push_back(npt);
  npIT_.push_back(npit);
}


 void BigNtuple::set_muInfo(TTree* tree){
   tree->Branch("mu_en" , &mu_en_);
   tree->Branch("mu_pt" , &mu_pt_);
   tree->Branch("mu_eta" , &mu_eta_);
   tree->Branch("mu_phi" , &mu_phi_);
   tree->Branch("mu_et" , &mu_et_);
   tree->Branch("mu_rhoIso",&mu_rhoIso_);
   tree->Branch("mu_charge" , &mu_charge_);
   tree->Branch("mu_trackiso" , &mu_trackiso_);
   tree->Branch("mu_pfSumChargedHadronPt" , &mu_pfSumChargedHadronPt_);
   tree->Branch("mu_pfSumNeutralHadronEt" , &mu_pfSumNeutralHadronEt_);
   tree->Branch("mu_PFSumPhotonEt" , &mu_PFSumPhotonEt_);
   tree->Branch("mu_pfSumPUPt" , &mu_pfSumPUPt_);
   tree->Branch("mu_pfSumPUPt04" ,&mu_pfSumPUPt04_);
   tree->Branch("mu_PFSumPhotonEt04",&mu_PFSumPhotonEt04_);
   tree->Branch("mu_pfSumChargedHadronPt04",&mu_pfSumChargedHadronPt04_);
   tree->Branch("mu_pfSumNeutralHadronEt04",&mu_pfSumNeutralHadronEt04_);
   tree->Branch("mu_numberOfValidMuonHits" , &mu_numberOfValidMuonHits_);
   tree->Branch("mu_emIso" , &mu_emIso_);
   tree->Branch("mu_hadIso" , &mu_hadIso_);
   tree->Branch("mu_segmentCompatibilityMuonBestTrack" ,&mu_segmentCompatibilityMuonBestTrack_);
   tree->Branch("mu_trkKinkMuonBestTrack" ,&mu_trkKinkMuonBestTrack_);
   tree->Branch("mu_chi2LocalPositionMuonBestTrack" ,&mu_chi2LocalPositionMuonBestTrack_);
   tree->Branch("mu_normalizedChi2" , &mu_normalizedChi2_);
   tree->Branch("mu_numberOfMatchedStations" , &mu_numberOfMatchedStations_);
   tree->Branch("mu_numberOfValidPixelHits" , &mu_numberOfValidPixelHits_);
   tree->Branch("mu_numberOftrackerLayersWithMeasurement" , &mu_numberOftrackerLayersWithMeasurement_);
   tree->Branch("mu_numberOfpixelLayersWithMeasurement" , &mu_numberOfpixelLayersWithMeasurement_);
   tree->Branch("mu_TrackQuality" , &mu_TrackQuality_);
   tree->Branch("mu_InnerTrackQuality" , &mu_InnerTrackQuality_);
   tree->Branch("mu_InnerTrackValidFraction" ,& mu_InnerTrackValidFraction_);
   tree->Branch("mu_pxTunePMuonBestTrack" , &mu_pxTunePMuonBestTrack_);
   tree->Branch("mu_pyTunePMuonBestTrack" , &mu_pyTunePMuonBestTrack_);
   tree->Branch("mu_pzTunePMuonBestTrack" , &mu_pzTunePMuonBestTrack_);
   tree->Branch("mu_pTunePMuonBestTrack" , &mu_pTunePMuonBestTrack_);
   tree->Branch("mu_etaTunePMuonBestTrack" , &mu_etaTunePMuonBestTrack_);
   tree->Branch("mu_ptTunePMuonBestTrack" , &mu_ptTunePMuonBestTrack_);
   tree->Branch("mu_phiTunePMuonBestTrack" , &mu_phiTunePMuonBestTrack_);
   tree->Branch("mu_thetaTunePMuonBestTrack" , &mu_thetaTunePMuonBestTrack_);
   tree->Branch("mu_chargeTunePMuonBestTrack" , &mu_chargeTunePMuonBestTrack_);
   tree->Branch("mu_dPToverPTTunePMuonBestTrack" , &mu_dPToverPTTunePMuonBestTrack_);
   tree->Branch("mu_absdxyTunePMuonBestTrack" , &mu_absdxyTunePMuonBestTrack_);
   tree->Branch("mu_absdxyErrorTunePMuonBestTrack" , &mu_absdxyErrorTunePMuonBestTrack_);
   tree->Branch("mu_absdxySigTunePMuonBestTrack" , &mu_absdxySigTunePMuonBestTrack_);
   tree->Branch("mu_absdzTunePMuonBestTrack" , &mu_absdzTunePMuonBestTrack_);
   tree->Branch("mu_absdzErrorTunePMuonBestTrack" , &mu_absdzErrorTunePMuonBestTrack_);
   tree->Branch("mu_absdzSigTunePMuonBestTrack" , &mu_absdzSigTunePMuonBestTrack_);
   tree->Branch("mu_3dIP",& mu_3dIP);
   tree->Branch("mu_3dIPSig",& mu_3dIPSig);
   tree->Branch("mu_2dIP",& mu_2dIP );
   tree->Branch("mu_2dIPSig",& mu_2dIPSig);

   tree->Branch("mu_recoDeltaBeta" , &mu_recoDeltaBeta_);
   tree->Branch("mu_recoiso" , &mu_recoiso_);
   tree->Branch("mu_isGlobalMuon" , &mu_isGlobalMuon_);
   tree->Branch("mu_isStandAloneMuon" , &mu_isStandAloneMuon_);
   tree->Branch("mu_isPFMuon" , &mu_isPFMuon_);
   tree->Branch("mu_isRPCMuon" , &mu_isRPCMuon_);
   tree->Branch("mu_isTrackerMuon" , &mu_isTrackerMuon_);
   tree->Branch("mu_isGoodMuon" , &mu_isGoodMuon_);
   tree->Branch("mu_isSoftMuon" , &mu_isSoftMuon_);
   tree->Branch("mu_isLooseMuon" , &mu_isLooseMuon_);
   tree->Branch("mu_isTightMuon" , &mu_isTightMuon_);
   /*
   tree->Branch("mu_STAnHits" , &mu_STAnHits_);
   tree->Branch("mu_STAnLost" , &mu_STAnLost_);
   tree->Branch("mu_STAnStationsWithAnyHits" , &mu_STAnStationsWithAnyHits_);
   tree->Branch("mu_STAnCscChambersWithAnyHits" , &mu_STAnCscChambersWithAnyHits_);
   tree->Branch("mu_STAnDtChambersWithAnyHits" , &mu_STAnDtChambersWithAnyHits_);
   tree->Branch("mu_STAnRpcChambersWithAnyHits" , &mu_STAnRpcChambersWithAnyHits_);
   tree->Branch("mu_STAinnermostStationWithAnyHits" , &mu_STAinnermostStationWithAnyHits_);
   tree->Branch("mu_STAoutermostStationWithAnyHits" , &mu_STAoutermostStationWithAnyHits_);
   tree->Branch("mu_STAnStationsWithValidHits" , &mu_STAnStationsWithValidHits_);
   tree->Branch("mu_STAnCscChambersWithValidHits" , &mu_STAnCscChambersWithValidHits_);
   tree->Branch("mu_STAnDtChambersWithValidHit" , &mu_STAnDtChambersWithValidHit_);
   tree->Branch("mu_STAnRpcChambersWithValidHits" , &mu_STAnRpcChambersWithValidHits_);
   tree->Branch("mu_STAnValidMuonHits" , &mu_STAnValidMuonHits_);
   tree->Branch("mu_STAnValidCscHits" , &mu_STAnValidCscHits_);
   tree->Branch("mu_STAnValidDtHits" , &mu_STAnValidDtHits_);
   tree->Branch("mu_STAnValidRpcHits" , &mu_STAnValidRpcHits_);
   tree->Branch("mu_STAinnermostStationWithValidHits" , &mu_STAinnermostStationWithValidHits_);
   tree->Branch("mu_STAoutermostStationWithValidHits" , &mu_STAoutermostStationWithValidHits_);
   */
   tree->Branch("mu_STATofDirection" , &mu_STATofDirection_);
   tree->Branch("mu_STATofNDof" , &mu_STATofNDof_);
   tree->Branch("mu_STATofTimeAtIpInOut" , &mu_STATofTimeAtIpInOut_);
   tree->Branch("mu_STATofTimeAtIpInOutErr" , &mu_STATofTimeAtIpInOutErr_);
   tree->Branch("mu_STATofTimeAtIpOutIn" , &mu_STATofTimeAtIpOutIn_);
   tree->Branch("mu_STATofTimeAtIpOutInErr" , &mu_STATofTimeAtIpOutInErr_);
   tree->Branch("mu_FirstGenMatch" , &mu_FirstGenMatch_);
   tree->Branch("mu_SecondGenMatch" , &mu_SecondGenMatch_);
   tree->Branch("mu_DecayChain" , &mu_GenMatchTest_);
   tree->Branch("mu_RPCTofDirection" , &mu_RPCTofDirection_);
   tree->Branch("mu_RPCTofNDof" , &mu_RPCTofNDof_);
   tree->Branch("mu_RPCTofTimeAtIpInOut" , &mu_RPCTofTimeAtIpInOut_);
   tree->Branch("mu_RPCTofTimeAtIpInOutErr" , &mu_RPCTofTimeAtIpInOutErr_);
   tree->Branch("mu_RPCTofTimeAtIpOutIn" , &mu_RPCTofTimeAtIpOutIn_);
   tree->Branch("mu_RPCTofTimeAtIpOutInErr" , &mu_RPCTofTimeAtIpOutInErr_);
}



void BigNtuple::fill_muInfo(const pat::Muon& mu, const reco::Vertex& pv, double Rho,  double match1 ,  double match2 , double match3){

  mu_isGlobalMuon_.push_back(mu.isGlobalMuon());
  mu_isPFMuon_.push_back(mu.isPFMuon());
  mu_isTrackerMuon_.push_back(mu.isTrackerMuon());
  mu_isRPCMuon_.push_back(mu.isRPCMuon());
  mu_isStandAloneMuon_.push_back(mu.isStandAloneMuon());
  mu_isSoftMuon_.push_back(mu.isSoftMuon(pv));
  mu_isLooseMuon_.push_back(mu.isLooseMuon());
  mu_isTightMuon_.push_back(mu.isTightMuon(pv));
  mu_en_.push_back(mu.energy());
  mu_et_.push_back(mu.et());
  mu_pt_.push_back(mu.pt());
  mu_eta_.push_back(mu.eta());
  mu_phi_.push_back(mu.phi());
  mu_charge_.push_back(mu.charge());
  mu_FirstGenMatch_.push_back(match1);
  mu_SecondGenMatch_.push_back(match2);
  mu_GenMatchTest_.push_back(match3);

  reco::TrackRef tunePTrack = mu.muonBestTrack();

  mu_ptTunePMuonBestTrack_.push_back(tunePTrack->pt()); // transverse momentum
  mu_dPToverPTTunePMuonBestTrack_.push_back(tunePTrack->ptError()/tunePTrack->pt()); // error calculation of transverse momentum
  mu_pxTunePMuonBestTrack_.push_back(tunePTrack->px()); //px component of the track
  mu_pyTunePMuonBestTrack_.push_back(tunePTrack->py()); //py component of the track
  mu_pzTunePMuonBestTrack_.push_back(tunePTrack->pz()); //pz component of the track
  mu_pTunePMuonBestTrack_.push_back(tunePTrack->p());   //magnitude of momentum vector
  mu_etaTunePMuonBestTrack_.push_back(tunePTrack->eta());
  mu_phiTunePMuonBestTrack_.push_back(tunePTrack->phi());
  mu_thetaTunePMuonBestTrack_.push_back(tunePTrack->theta());
  mu_chargeTunePMuonBestTrack_.push_back(tunePTrack->charge());
  mu_absdxyTunePMuonBestTrack_.push_back(fabs(tunePTrack->dxy(pv.position()))); //transvers  impact parameter  w.r.t. the primary vertex
  mu_absdxyErrorTunePMuonBestTrack_.push_back(fabs(tunePTrack->dxyError())); //transvers  impact parameter  w.r.t. the primary vertex
  mu_absdxySigTunePMuonBestTrack_.push_back(fabs(tunePTrack->dxy(pv.position()))/fabs(tunePTrack->dxyError()));
  mu_absdzTunePMuonBestTrack_.push_back(fabs(tunePTrack->dz(pv.position()))); // longitudinal impact parameter  w.r.t. the primary vertex
  mu_absdzErrorTunePMuonBestTrack_.push_back(fabs(tunePTrack->dzError())); // longitudinal impact parameter  w.r.t. the primary vertex
  mu_absdzSigTunePMuonBestTrack_.push_back(fabs(tunePTrack->dz(pv.position()))/fabs(tunePTrack->dzError()));
  mu_TrackQuality_.push_back(tunePTrack->quality(reco::TrackBase::highPurity));

  mu_3dIP.push_back(mu.dB(pat::Muon::PV3D));
  mu_3dIPSig.push_back(mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D));
  mu_2dIP.push_back(mu.dB());
  mu_2dIPSig.push_back(mu.dB()/mu.edB());


  mu_rhoIso_.push_back(Rho); //transverse momentum per unit area

  if(mu.globalTrack().isNonnull() ) {
    mu_normalizedChi2_.push_back(mu.globalTrack()->normalizedChi2());
    mu_numberOfValidPixelHits_.push_back(mu.innerTrack()->hitPattern().numberOfValidPixelHits());
    mu_numberOfValidMuonHits_.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
    mu_numberOftrackerLayersWithMeasurement_.push_back(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
    mu_numberOfMatchedStations_.push_back(mu.numberOfMatchedStations());
    mu_numberOfpixelLayersWithMeasurement_.push_back(mu.innerTrack()->hitPattern().pixelLayersWithMeasurement());
    mu_InnerTrackQuality_.push_back(mu.innerTrack()->quality(reco::TrackBase::highPurity));
    mu_InnerTrackValidFraction_.push_back(mu.innerTrack()->validFraction());
  }
  else{
    mu_normalizedChi2_.push_back(-999);
    mu_numberOfValidPixelHits_.push_back(-999);
    mu_numberOfValidMuonHits_.push_back(-999);
    mu_numberOftrackerLayersWithMeasurement_.push_back(-999);
    mu_numberOfMatchedStations_.push_back(-999);
    mu_numberOfpixelLayersWithMeasurement_.push_back(-999);
    mu_InnerTrackQuality_.push_back(-999);
    mu_InnerTrackValidFraction_.push_back(-999);
  }
  /*
  if(mu.standAloneMuon().isNonnull() ) {
    mu_STAnHits_.push_back(mu.standAloneMuon()->numberOfValidHits());
    mu_STAnLost_.push_back(mu.standAloneMuon()->numberOfLostHits());
    mu_STAnStationsWithAnyHits_.push_back(mu.standAloneMuon()->hitPattern().muonStationsWithAnyHits());
    mu_STAnCscChambersWithAnyHits_.push_back(mu.standAloneMuon()->hitPattern().cscStationsWithAnyHits()); //csc chambers in track fit
    mu_STAnDtChambersWithAnyHits_.push_back(mu.standAloneMuon()->hitPattern().dtStationsWithAnyHits()); //dt chambers in track fit
    mu_STAnRpcChambersWithAnyHits_.push_back(mu.standAloneMuon()->hitPattern().rpcStationsWithAnyHits()); //rpc chambers in track fit
    mu_STAinnermostStationWithAnyHits_.push_back(mu.standAloneMuon()->hitPattern().innermostMuonStationWithAnyHits());
    mu_STAoutermostStationWithAnyHits_.push_back(mu.standAloneMuon()->hitPattern().outermostMuonStationWithAnyHits());
    mu_STAnCscChambersWithValidHits_.push_back(mu.standAloneMuon()->hitPattern().cscStationsWithValidHits()); //csc chambers anywhere near track
    mu_STAnDtChambersWithValidHit_.push_back(mu.standAloneMuon()->hitPattern().dtStationsWithValidHits()); //dt chambers anywhere near track
    mu_STAnRpcChambersWithValidHits_.push_back(mu.standAloneMuon()->hitPattern().rpcStationsWithValidHits()); //rpc chambers anywhere near track
    mu_STAnValidCscHits_.push_back(mu.standAloneMuon()->hitPattern().numberOfValidMuonCSCHits()); //CSC hits anywhere near track
    mu_STAnValidDtHits_.push_back(mu.standAloneMuon()->hitPattern().numberOfValidMuonDTHits()); //DT hits anywhere near track
    mu_STAnValidRpcHits_.push_back(mu.standAloneMuon()->hitPattern().numberOfValidMuonRPCHits()); //RPC hits anywhere near track
    mu_STAnValidMuonHits_.push_back(mu.standAloneMuon()->hitPattern().numberOfValidMuonHits()); //muon hits anywhere near track
    mu_STAinnermostStationWithValidHits_.push_back(mu.standAloneMuon()->hitPattern().innermostMuonStationWithValidHits());
    mu_STAoutermostStationWithValidHits_.push_back(mu.standAloneMuon()->hitPattern().outermostMuonStationWithValidHits());
    mu_STAnStationsWithValidHits_.push_back(mu.standAloneMuon()->hitPattern().muonStationsWithValidHits());
  }
  else{
    mu_STAnHits_.push_back(-999);
    mu_STAnLost_.push_back(-999);
    mu_STAnStationsWithAnyHits_.push_back(-999);
    mu_STAnCscChambersWithAnyHits_.push_back(-999);
    mu_STAnDtChambersWithAnyHits_.push_back(-999);
    mu_STAnRpcChambersWithAnyHits_.push_back(-999);
    mu_STAinnermostStationWithAnyHits_.push_back(-999);
    mu_STAoutermostStationWithAnyHits_.push_back(-999);
    mu_STAnCscChambersWithValidHits_.push_back(-999);
    mu_STAnDtChambersWithValidHit_.push_back(-999);
    mu_STAnRpcChambersWithValidHits_.push_back(-999);
    mu_STAnValidCscHits_.push_back(-999);
    mu_STAnValidDtHits_.push_back(-999);
    mu_STAnValidRpcHits_.push_back(-999);
    mu_STAnValidMuonHits_.push_back(-999);
    mu_STAinnermostStationWithValidHits_.push_back(-999);
    mu_STAoutermostStationWithValidHits_.push_back(-999);
    mu_STAnStationsWithValidHits_.push_back(-999);
  }
  */
  //time info
  reco::MuonTime tofAll = mu.time();
  reco::MuonTime tofRPC = mu.rpcTime();
  mu_STATofDirection_.push_back(tofAll.direction());
  mu_STATofNDof_.push_back(tofAll.nDof);
  mu_STATofTimeAtIpInOut_.push_back(tofAll.timeAtIpInOut);
  mu_STATofTimeAtIpInOutErr_.push_back(tofAll.timeAtIpInOutErr);
  mu_STATofTimeAtIpOutIn_.push_back(tofAll.timeAtIpOutIn);
  mu_STATofTimeAtIpOutInErr_.push_back(tofAll.timeAtIpOutInErr);

  mu_RPCTofDirection_.push_back(tofRPC.direction());
  mu_RPCTofNDof_.push_back(tofRPC.nDof);
  mu_RPCTofTimeAtIpInOut_.push_back(tofRPC.timeAtIpInOut);
  mu_RPCTofTimeAtIpInOutErr_.push_back(tofRPC.timeAtIpInOutErr);
  mu_RPCTofTimeAtIpOutIn_.push_back(tofRPC.timeAtIpOutIn);
  mu_RPCTofTimeAtIpOutInErr_.push_back(tofRPC.timeAtIpOutInErr);

  //============= Parameters related to detector isolation =====================
  double charged   = mu.pfIsolationR04().sumChargedHadronPt;
  double neutral   = mu.pfIsolationR04().sumNeutralHadronEt;
  double pileup    = mu.pfIsolationR04().sumPUPt;
  double sumPhotonEt = mu.pfIsolationR04().sumPhotonEt; //Sum Et of PF photonds
  double Mu_iso = 1.0*(charged  +  neutral + sumPhotonEt )/mu.pt(); //recommended this be < 0.20 (loose) or < 0.12 (tight)
  double deltaBeta = (charged + std::max(0.0, neutral+sumPhotonEt-0.5*pileup))/mu.pt();
  mu_recoDeltaBeta_.push_back(deltaBeta); //Delta Beta
  mu_recoiso_.push_back(Mu_iso);
  mu_emIso_.push_back(mu.isolationR03().emEt);
  mu_hadIso_.push_back(mu.isolationR03().hadEt);
  mu_trackiso_.push_back(mu.isolationR03().sumPt);
  mu_segmentCompatibilityMuonBestTrack_.push_back(mu.segmentCompatibility());
  mu_trkKinkMuonBestTrack_.push_back(mu.combinedQuality().trkKink);
  mu_chi2LocalPositionMuonBestTrack_.push_back(mu.combinedQuality().chi2LocalPosition);
  //============= Parameters related to PF isolation =====================
  mu_pfSumPUPt_.push_back(mu.pfIsolationR03().sumPUPt);
  mu_PFSumPhotonEt_.push_back(mu.pfIsolationR03().sumPhotonEt);
  mu_pfSumChargedHadronPt_.push_back(mu.pfIsolationR03().sumChargedHadronPt);
  mu_pfSumNeutralHadronEt_.push_back(mu.pfIsolationR03().sumNeutralHadronEt);

  mu_pfSumPUPt04_.push_back(mu.pfIsolationR04().sumPUPt);
  mu_PFSumPhotonEt04_.push_back(mu.pfIsolationR04().sumPhotonEt);
  mu_pfSumChargedHadronPt04_.push_back(mu.pfIsolationR04().sumChargedHadronPt);
  mu_pfSumNeutralHadronEt04_.push_back(mu.pfIsolationR04().sumNeutralHadronEt);

}

  void BigNtuple::set_sv_Info(TTree* tree){

    tree->Branch("sv_munTracks" , &sv_numTracks_);
    tree->Branch("sv_X" , &sv_x_);
    tree->Branch("sv_Y" , &sv_y_);
    tree->Branch("sv_Z" , &sv_z_);
    tree->Branch("sv_xErr" , &sv_xErr_);
    tree->Branch("sv_yErr" , &sv_yErr_);
    tree->Branch("sv_zErr" , &sv_zErr_);
    tree->Branch("sv_LxySig" , &sv_LxySig_);
    tree->Branch("sv_LxyzSig" , &sv_LxyzSig_);
    tree->Branch("sv_Lxy" , &sv_Lxy_);
    tree->Branch("sv_Lxyz" , &sv_Lxyz_);
    tree->Branch("sv_mass" , &sv_mass_);
    tree->Branch("sv_charge" , &sv_charge_);
    tree->Branch("sv_eta" , &sv_eta_);
    tree->Branch("sv_phi" , &sv_phi_);
    tree->Branch("sv_pt" , &sv_pt_);
    tree->Branch("sv_p" , &sv_p_);
    tree->Branch("sv_px" , &sv_px_);
    tree->Branch("sv_py" , &sv_py_);
    tree->Branch("sv_pz" , &sv_pz_);
    tree->Branch("sv_energy" , &sv_energy_);
    tree->Branch("sv_Beta" , &sv_Beta_);
    tree->Branch("sv_Gamma" , &sv_Gamma_);
    tree->Branch("sv_CTau0" , &sv_CTau0_);
    tree->Branch("sv_NDof" , &sv_NDof_);
    tree->Branch("sv_Chi2" , &sv_Chi2_);
    tree->Branch("sv_Angle3D" , &sv_Angle3D_);
    tree->Branch("sv_Angle2D" , &sv_Angle2D_);
    tree->Branch("sv_tracks_charge" , &sv_tracks_charge_);
    tree->Branch("sv_tracks_eta" , &sv_tracks_eta_);
    tree->Branch("sv_tracks_phi" , &sv_tracks_phi_);
    tree->Branch("sv_tracks_pt" , &sv_tracks_pt_);
    tree->Branch("sv_tracks_p" , &sv_tracks_p_);
    tree->Branch("sv_tracks_dxySig" , &sv_tracks_dxySig_);
    tree->Branch("sv_tracks_dxy" , &sv_tracks_dxy_);
    tree->Branch("sv_tracks_dxyz" , &sv_tracks_dxyz_);
    tree->Branch("sv_tracks_Sumcharge" , &sv_tracks_Sumcharge_);
    tree->Branch("sv_tracks_Sumpt" , &sv_tracks_Sumpt_);
    tree->Branch("sv_match_dxyz" , &sv_match_dxyz_);
    tree->Branch("sv_match_dxy" , &sv_match_dxy_);
    tree->Branch("sv_lx" , & sv_lx_);
    tree->Branch("sv_ly" , & sv_ly_);
    tree->Branch("sv_lz" , & sv_lz_);
    tree->Branch("sv_hasMuon", & sv_hasMuon_);
  }

void BigNtuple::fill_sv_Info(const reco::Vertex& bestVertex, const reco::Vertex& pv , std::pair<float, float> match, bool lept){

  float svChi2 = bestVertex.chi2();
  float svNDof = bestVertex.ndof();

  //flight distance from the firstPV
  float x  = bestVertex.x(), y = bestVertex.y(), z = bestVertex.z();
  float dx = x - pv.x() , dy = y - pv.y(), dz = z - pv.z();

  //build the total error
  float sv_xErr = bestVertex.xError(), sv_yErr = bestVertex.yError(), sv_zErr = bestVertex.zError();
  float pv_xErr = pv.xError(), pv_yErr = pv.yError(), pv_zErr = pv.zError();
  float l_xErr   = std::sqrt(sv_xErr * sv_xErr + pv_xErr * pv_xErr), l_yErr = std::sqrt(sv_yErr * sv_yErr + pv_yErr * pv_yErr), l_zErr = std::sqrt(sv_zErr * sv_zErr + pv_zErr * pv_zErr);

  // mother beta, gamma, ctau
  float beta_hnl  = bestVertex.p4().P() / bestVertex.p4().energy();
  float gamma_hnl = bestVertex.p4().energy() / bestVertex.p4().mass();

  TVector3 pvVector3D(pv.x(), pv.y(), pv.z());
  TVector3 pvVector2D(pv.x(), pv.y(), 0);
  TVector3 svVector3D(x, y, z);
  TVector3 svVector2D(x, y, 0);

  //projection in x, y (and z) of the sum of the tracks which form the vertex
  TVector3 sv_momentum_3D( bestVertex.p4().x(), bestVertex.p4().y(), bestVertex.p4().z());
  TVector3 sv_momentum_2D( bestVertex.p4().x(), bestVertex.p4().y(), 0); 

  float sign2D =  (sv_momentum_2D * (svVector2D - pvVector2D)) > 0 ? -1: 1;
  float sign3D =  (sv_momentum_3D * (svVector3D - pvVector3D)) > 0 ? -1: 1;

  TVector3 pvToVertex3D( sign3D * dx, sign3D * dy, sign3D * dz);
  TVector3 pvToVertex2D( sign2D * dx, sign2D * dy, 0);

  pvTosv_rho_ = pvToVertex3D.Mag();
  pvTosv_phi_ = pvToVertex3D.Phi();
  pvTosv_theta_ = pvToVertex3D.Theta();

  float svAngle3D = pvToVertex3D.Angle(sv_momentum_3D);
  float svAngle2D = pvToVertex2D.Angle(sv_momentum_2D);
  
  //bestVertex info
  best_sv_px_ = bestVertex.p4().px();
  best_sv_py_ = bestVertex.p4().py();
  best_sv_pz_ = bestVertex.p4().pz();
  best_sv_pt_ = bestVertex.p4().pt();
  best_sv_energy_ = bestVertex.p4().energy();

  best_sv_recox_ = bestVertex.x();
  best_sv_recoy_ = bestVertex.y();
  best_sv_recoz_ = bestVertex.z();
    

  sv_lx_.push_back(dx);
  sv_ly_.push_back(dy);
  sv_lz_.push_back(dz);
  
  sv_numTracks_.push_back(bestVertex.nTracks());
  sv_x_.push_back(x);
  sv_y_.push_back(y);
  sv_z_.push_back(z);
  sv_xErr_.push_back(sv_xErr);
  sv_yErr_.push_back(sv_yErr);
  sv_zErr_.push_back(sv_zErr);
  sv_Lxy_.push_back(std::sqrt( dx * dx + dy * dy ));
  sv_Lxyz_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz ));
  sv_LxySig_.push_back(std::sqrt( dx * dx + dy * dy ) / std::sqrt(l_xErr * l_xErr + l_yErr * l_yErr));
  sv_LxyzSig_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz) / std::sqrt(l_xErr * l_xErr + l_yErr * l_yErr + l_zErr * l_zErr));
  sv_mass_.push_back(bestVertex.p4().mass());
  sv_eta_.push_back(bestVertex.p4().eta());
  sv_phi_.push_back(bestVertex.p4().phi());
  sv_pt_.push_back(bestVertex.p4().pt());
  sv_p_.push_back(bestVertex.p4().P());
  sv_px_.push_back(bestVertex.p4().px());
  sv_py_.push_back(bestVertex.p4().py());
  sv_pz_.push_back(bestVertex.p4().pz());
  sv_energy_.push_back(bestVertex.p4().energy());
  sv_Beta_.push_back(beta_hnl);
  sv_Gamma_.push_back(gamma_hnl);
  sv_CTau0_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz) / (beta_hnl * gamma_hnl));
  sv_NDof_.push_back(svNDof);
  sv_Chi2_.push_back(svChi2);
  sv_Angle3D_.push_back(svAngle3D);
  sv_Angle2D_.push_back(svAngle2D);
  sv_hasMuon_.push_back(lept);


  int ch = 0;
  float pt = 0;
  sv_tracks_charge_.emplace_back();
  sv_tracks_eta_.emplace_back();
  sv_tracks_phi_.emplace_back();
  sv_tracks_pt_.emplace_back();
  sv_tracks_p_.emplace_back();
  sv_tracks_dxySig_.emplace_back();
  sv_tracks_dxy_.emplace_back();
  sv_tracks_dxyz_.emplace_back();

  reco::Vertex::trackRef_iterator tt = bestVertex.tracks_begin();
  for(; tt != bestVertex.tracks_end(); ++tt) {

    sv_tracks_charge_.back().push_back((*tt)->charge());
    sv_tracks_eta_.back().push_back((*tt)->eta());
    sv_tracks_phi_.back().push_back((*tt)->phi());
    sv_tracks_pt_.back().push_back((*tt)->pt());
    sv_tracks_p_.back().push_back((*tt)->p());
    sv_tracks_dxySig_.back().push_back(fabs((*tt)->dxy(pv.position()))/fabs((*tt)->dxyError()));
    sv_tracks_dxy_.back().push_back((*tt)->dxy(pv.position()));

    ROOT::Math::SVector<double, 3> lxyz1((*tt)->vx()-pv.position().x(), (*tt)->vy()-pv.position().y(), (*tt)->vz()-pv.position().z());
    float dxyz = (float)ROOT::Math::Mag(lxyz1); // magntude of the vector

    sv_tracks_dxyz_.back().push_back(dxyz);

    ch+=(*tt)->charge();
    pt+=(*tt)->pt();
  }

  sv_tracks_Sumcharge_.push_back(ch);
  sv_tracks_Sumpt_.push_back(pt);
  sv_match_dxyz_.push_back(match.first);
  sv_match_dxy_.push_back(match.second);

}

double BigNtuple::reducedPdgId( int pdgId ){
  static const std::map< unsigned, double > pdgIdMap = {
    { 0, 0.},
    { 1, 0.125},
    { 2, 0.25},
    { 11, 0.375},
    { 13, 0.5},
    { 22, 0.625},
    { 130, 0.75},
    { 211, 0.875}
  };
  auto entry = pdgIdMap.find( fabs( pdgId ) );
  if( entry != pdgIdMap.cend() ){
    return entry->second;
  } else {
    return 1;
  }
}

double BigNtuple::catchNanOrInf( double value ){
  if( std::isnan( value ) ){
    return -1;
  } else if( std::isinf( value ) ){
    return -1;
  } else{
    return value;
  }
}

void BigNtuple::set_jetInfo(TTree* tree){

  tree->Branch("jetPt_JECUp", &jetPt_JECUp_);
  tree->Branch("jetPt_JECDown", &jetPt_JECDown_);
  tree->Branch("jetSmearedPt", &jetSmearedPt_);
  tree->Branch("jetSmearedPt_JERUp", &jetSmearedPt_JERUp_);
  tree->Branch("jetSmearedPt_JERDown", &jetSmearedPt_JERDown_);
  tree->Branch("jetSmearedPt_unUp", &jetSmearedPt_unUp_);
  tree->Branch("jetSmearedPt_unDown", &jetSmearedPt_unDown_);
  tree->Branch("jet_charge" , &jet_charge_);
  tree->Branch("jet_et" , &jet_et_);
  tree->Branch("jet_pt" , &jet_pt_);
  tree->Branch("jet_px" , &jet_px_);
  tree->Branch("jet_py" , &jet_py_);
  tree->Branch("jet_pz" , &jet_pz_);
  tree->Branch("jet_eta" , &jet_eta_);
  tree->Branch("jet_phi" , &jet_phi_);
  tree->Branch("jet_theta" , &jet_theta_);
  tree->Branch("jet_en" , &jet_en_);
  tree->Branch("jet_chargedEmEnergy" , &jet_chargedEmEnergy_);
  tree->Branch("jet_chargedEmEnergyFraction" ,&jet_chargedEmEnergyFraction_);
  tree->Branch("jet_neutralEmEnergyFraction" , &jet_neutralEmEnergyFraction_);
  tree->Branch("jet_chargedHadronEnergy" , &jet_chargedHadronEnergy_);
  tree->Branch("jet_chargedHadronEnergyFraction" , &jet_chargedHadronEnergyFraction_);
  tree->Branch("jet_neutralHadronEnergyFraction" , &jet_neutralHadronEnergyFraction_);
  tree->Branch("jet_chargedMuEnergy" , &jet_chargedMuEnergy_);
  tree->Branch("jet_chargedMuEnergyFraction" , &jet_chargedMuEnergyFraction_);
  tree->Branch("jet_numberOfDaughters" , &jet_numberOfDaughters_);
  tree->Branch("jet_muonEnergy" , &jet_muonEnergy_);
  tree->Branch("jet_muonEnergyFraction" , &jet_muonEnergyFraction_);
  tree->Branch("jet_muonMultiplicity" , &jet_muonMultiplicity_);
  tree->Branch("jet_neutralEmEnergy" , &jet_neutralEmEnergy_);
  tree->Branch("jet_neutralHadronEnergy" , &jet_neutralHadronEnergy_);
  tree->Branch("jet_neutralHadronMultiplicity" , &jet_neutralHadronMultiplicity_);
  tree->Branch("jet_neutralMultiplicity" , &jet_neutralMultiplicity_);
  tree->Branch("jet_chargedMultiplicity" , &jet_chargedMultiplicity_);
  tree->Branch("jet_pileUpid", &jet_pileUpid_);
  tree->Branch("jet_ptuncorrected", &jet_ptuncorrected_);
  tree->Branch("jet_L1ptcorrection", &jet_L1ptcorrection_);
  tree->Branch("jet_L2ptcorrection", &jet_L2ptcorrection_);
  tree->Branch("jet_L3ptcorrection", &jet_L3ptcorrection_);
  tree->Branch("jet_CsvV2", &jet_CsvV2_);
  tree->Branch("jet_DeepCsv_udsg", &jet_DeepCsv_udsg_);
  tree->Branch("jet_DeepCsv_b", &jet_DeepCsv_b_);
  tree->Branch("jet_DeepCsv_c", &jet_DeepCsv_c_);
  tree->Branch("jet_DeepCsv_bb", &jet_DeepCsv_bb_);
  tree->Branch("jet_HadronFlavor", &jet_HadronFlavor_);

  tree->Branch("jet_daughter_Pt",&jet_daughter_Pt_);
  tree->Branch("jet_daughter_Eta",&jet_daughter_Eta_);
  tree->Branch("jet_daughter_Phi",&jet_daughter_Phi_);
  tree->Branch("jet_daughter_Mass",&jet_daughter_Mass_);
  tree->Branch("jet_daughter_PdgId",&jet_daughter_PdgId_);
  tree->Branch("jet_daughter_PdgIdReduced",&jet_daughter_PdgIdReduced_);
  tree->Branch("jet_daughter_Charge",&jet_daughter_Charge_);
  tree->Branch("jet_daughter_dxy",&jet_daughter_dxy_);
  tree->Branch("jet_daughter_dz",&jet_daughter_dz_);
  tree->Branch("jet_daughter_dxyErr",&jet_daughter_dxyErr_);
  tree->Branch("jet_daughter_dzErr",&jet_daughter_dzErr_);
  tree->Branch("jet_daughter_NumberOfHits",&jet_daughter_NumberOfHits_);
  tree->Branch("jet_daughter_NumberOfPixelHits",&jet_daughter_NumberOfPixelHits_);
  tree->Branch("jet_daughter_HasTrack",&jet_daughter_HasTrack_);

}



void BigNtuple::fill_jetInfo(const pat::Jet& jet, float smeared,float smearedUp ,float smearedDown ,double un , double unSmeared ){

  /*
  JME::JetResolution resolution;
  JME::JetParameters parameters_1;
  JME::JetResolutionScaleFactor res_sf;

  parameters_1.setJetPt(jet.pt());
  parameters_1.setJetEta(jet.eta());

  float r = resolution.getResolution(parameters_1);

  JME::JetParameters parameters = {{JME::Binning::JetEta, jet.eta()}, {JME::Binning::Rho, 2}};

  float sf = res_sf.getScaleFactor(parameters);
  float sf_up = res_sf.getScaleFactor(parameters, Variation::UP);
  float sf_down = res_sf.getScaleFactor(parameters, Variation::DOWN);
  */

  jet_pt_.push_back(jet.pt());
  jet_px_.push_back(jet.px());
  jet_py_.push_back(jet.py());
  jet_pz_.push_back(jet.pz());
  jetPt_JECUp_.push_back(jet.pt()*(1 + un));
  jetPt_JECDown_.push_back(jet.pt()*(1 - un));
  jetSmearedPt_.push_back(smeared);
  jetSmearedPt_JERUp_.push_back(smearedUp);
  jetSmearedPt_JERDown_.push_back(smearedDown);
  jetSmearedPt_unUp_.push_back(jet.pt()*(1.+unSmeared));
  jetSmearedPt_unDown_.push_back(jet.pt()*(1.-unSmeared));

  jet_charge_.push_back(jet.charge());
  jet_et_.push_back(jet.et());
  jet_eta_.push_back(jet.eta());
  jet_phi_.push_back(jet.phi());
  jet_theta_.push_back(jet.theta());
  jet_en_.push_back(jet.energy());
  jet_chargedEmEnergy_.push_back(jet.chargedEmEnergy());
  jet_chargedEmEnergyFraction_.push_back(jet.chargedEmEnergyFraction());
  jet_neutralEmEnergyFraction_.push_back(jet.neutralEmEnergyFraction());
  jet_chargedHadronEnergy_.push_back(jet.chargedHadronEnergy());
  jet_chargedHadronEnergyFraction_.push_back(jet.chargedHadronEnergyFraction());
  jet_neutralHadronEnergyFraction_.push_back(jet.neutralHadronEnergyFraction());
  jet_chargedMuEnergy_.push_back(jet.chargedMuEnergy());
  jet_chargedMuEnergyFraction_.push_back(jet.chargedMuEnergyFraction());
  jet_chargedMultiplicity_.push_back(jet.chargedMultiplicity());
  jet_numberOfDaughters_.push_back(jet.numberOfDaughters());
  jet_muonEnergy_.push_back(jet.muonEnergy());
  jet_muonEnergyFraction_.push_back(jet.muonEnergyFraction());
  jet_muonMultiplicity_.push_back(jet.muonMultiplicity());
  jet_neutralEmEnergy_.push_back(jet.neutralEmEnergy());
  jet_neutralHadronEnergy_.push_back(jet.neutralHadronEnergy());
  jet_neutralHadronMultiplicity_.push_back(jet.neutralHadronMultiplicity());
  jet_neutralMultiplicity_.push_back(jet.neutralMultiplicity());
  jet_pileUpid_.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
  jet_L1ptcorrection_.push_back(jet.correctedP4("L3Absolute").Pt());
  jet_L2ptcorrection_.push_back(jet.correctedP4("L3Absolute").Pt());
  jet_L3ptcorrection_.push_back(jet.correctedP4("L3Absolute").Pt());
  jet_ptuncorrected_.push_back(jet.correctedP4("Uncorrected").Pt());

  float CsvV2 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
  float DeepCsv_udsg = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");
  float DeepCsv_b = jet.bDiscriminator("pfDeepCSVJetTags:probb");
  float DeepCsv_c = jet.bDiscriminator("pfDeepCSVJetTags:probc");
  float DeepCsv_bb = jet.bDiscriminator("pfDeepCSVJetTags:probbb");

  jet_CsvV2_.push_back(CsvV2);
  jet_DeepCsv_udsg_.push_back(DeepCsv_udsg);
  jet_DeepCsv_b_.push_back(DeepCsv_b);
  jet_DeepCsv_c_.push_back(DeepCsv_c);
  jet_DeepCsv_bb_.push_back(DeepCsv_bb);
  jet_HadronFlavor_.push_back(jet.hadronFlavour());

  jet_daughter_Pt_.emplace_back();
  jet_daughter_Eta_.emplace_back();
  jet_daughter_Phi_.emplace_back();
  jet_daughter_Mass_.emplace_back();
  jet_daughter_PdgId_.emplace_back();
  jet_daughter_PdgIdReduced_.emplace_back();
  jet_daughter_Charge_.emplace_back();
  jet_daughter_dxy_.emplace_back();
  jet_daughter_dz_.emplace_back();
  jet_daughter_dxyErr_.emplace_back();
  jet_daughter_dzErr_.emplace_back();
  jet_daughter_NumberOfHits_.emplace_back();
  jet_daughter_NumberOfPixelHits_.emplace_back();
  jet_daughter_HasTrack_.emplace_back();

  for(unsigned d = 0; d < jet.numberOfDaughters(); ++d){
    const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet.daughter(d);

    jet_daughter_Pt_.back().push_back(daughter->pt());
    jet_daughter_Eta_.back().push_back(daughter->eta());
    jet_daughter_Phi_.back().push_back(daughter->phi());
    jet_daughter_Mass_.back().push_back(daughter->mass());
    jet_daughter_PdgId_.back().push_back(daughter->pdgId());
    jet_daughter_PdgIdReduced_.back().push_back(reducedPdgId(daughter->pdgId()));
    jet_daughter_Charge_.back().push_back(daughter->charge());

    if( daughter->hasTrackDetails() ){
      jet_daughter_dxy_.back().push_back(fabs(daughter->dxy()));
      jet_daughter_dz_.back().push_back(fabs(daughter->dz()));
      jet_daughter_dxyErr_.back().push_back(catchNanOrInf(fabs(daughter->dxyError())));
      jet_daughter_dzErr_.back().push_back(catchNanOrInf(fabs(daughter->dzError())));
      jet_daughter_NumberOfHits_.back().push_back(daughter->numberOfHits());
      jet_daughter_NumberOfPixelHits_.back().push_back(daughter->numberOfPixelHits());
      jet_daughter_HasTrack_.back().push_back(true);
    } 
  }


}

void BigNtuple::set_eleInfo(TTree* tree){

  tree->Branch("ele_Et",&ele_Et_);
  tree->Branch("ele_EtFromCaloEn",&ele_EtFromCaloEn_);
  tree->Branch("ele_pt",&ele_pt_);
  tree->Branch("ele_px",&ele_px_);
  tree->Branch("ele_py",&ele_py_);
  tree->Branch("ele_pz",&ele_pz_);
  tree->Branch("ele_etaSC",&ele_etaSC_);
  tree->Branch("ele_phiSC",&ele_phiSC_);
  tree->Branch("ele_phiWidth",&ele_phiWidth_);
  tree->Branch("ele_etaWidth",&ele_etaWidth_);
  tree->Branch("ele_energySC",&ele_energySC_);
  tree->Branch("ele_thetaSC",&ele_thetaSC_);
  tree->Branch("ele_preshowerEnergySC",&ele_preshowerEnergySC_);
  tree->Branch("ele_etaTrack",&ele_etaTrack_);
  tree->Branch("ele_phiTrack",&ele_phiTrack_);
  tree->Branch("ele_thetaTrack",&ele_thetaTrack_);
  tree->Branch("ele_eta",&ele_eta_);
  tree->Branch("ele_phi",&ele_phi_);
  tree->Branch("ele_theta",&ele_theta_);
  tree->Branch("ele_energy",&ele_energy_);
  tree->Branch("ele_x",&ele_x_);
  tree->Branch("ele_y",&ele_y_);
  tree->Branch("ele_z",&ele_z_);
  tree->Branch("ele_e2x5Max",&ele_e2x5Max_);
  tree->Branch("ele_e1x5",&ele_e1x5_);
  tree->Branch("ele_e5x5",&ele_e5x5_);
  tree->Branch("ele_e2x5MaxOver5x5",&ele_e2x5MaxOver5x5_);
  tree->Branch("ele_e1x5Over5x5",&ele_e1x5Over5x5_);
  tree->Branch("ele_sigmaIetaIetaFull5x5",&ele_sigmaIetaIetaFull5x5_);
  tree->Branch("ele_e2x5MaxFull5x5",&ele_e2x5MaxFull5x5_);
  tree->Branch("ele_e1x5Full5x5",&ele_e1x5Full5x5_);
  tree->Branch("ele_e5x5Full5x5",&ele_e5x5Full5x5_);
  tree->Branch("ele_e2x5MaxOver5x5Full5x5",&ele_e2x5MaxOver5x5Full5x5_);
  tree->Branch("ele_e1x5Over5x5Full5x5",&ele_e1x5Over5x5Full5x5_);
  tree->Branch("ele_zTrackPositionAtVtx",&ele_zTrackPositionAtVtx_);
  tree->Branch("ele_hadronicOverEm",&ele_hadronicOverEm_);
  tree->Branch("ele_deltaEtaInSC",&ele_deltaEtaInSC_);
  tree->Branch("ele_deltaPhiInSC",&ele_deltaPhiInSC_);
  tree->Branch("ele_deltaEtaInSeedCluster",&ele_deltaEtaInSeedCluster_);
  tree->Branch("ele_deltaPhiInSeedCluster",&ele_deltaPhiInSeedCluster_);
  tree->Branch("ele_sigmaIetaIeta",&ele_sigmaIetaIeta_);
  tree->Branch("ele_rawId",&ele_rawId_);
  tree->Branch("ele_ieta",&ele_ieta_);
  tree->Branch("ele_e2x5Right",&ele_e2x5Right_);
  tree->Branch("ele_e2x5Left",&ele_e2x5Left_);
  tree->Branch("ele_e2x5Top",&ele_e2x5Top_);
  tree->Branch("ele_e2x5Bottom",&ele_e2x5Bottom_);
  tree->Branch("ele_eMax",&ele_eMax_);
  tree->Branch("ele_eRight",&ele_eRight_);
  tree->Branch("ele_eLeft",&ele_eLeft_);
  tree->Branch("ele_eTop",&ele_eTop_);
  tree->Branch("ele_eBottom",&ele_eBottom_);
  tree->Branch("ele_e3x3",&ele_e3x3_);
  tree->Branch("ele_frac51",&ele_frac51_);
  tree->Branch("ele_frac15",&ele_frac15_);
  tree->Branch("ele_dxy",&ele_dxy_);
  tree->Branch("ele_dz",&ele_dz_);
  tree->Branch("ele_3dIP",& ele_3dIP);
  tree->Branch("ele_3dIPSig",& ele_3dIPSig); 
  tree->Branch("ele_2dIP",& ele_2dIP );
  tree->Branch("ele_2dIPSig",& ele_2dIPSig);

  tree->Branch("ele_isEcalDrivenSeed",&ele_isEcalDrivenSeed_);
  tree->Branch("ele_isPassConversionVeto",&ele_isPassConversionVeto_);
  tree->Branch("ele_charge",&ele_charge_);
  tree->Branch("ele_rhoIso",&ele_rhoIso_);
  tree->Branch("ele_nbOfMissingHits",&ele_nbOfMissingHits_);
  tree->Branch("ele_fbrem",&ele_fbrem_);
  tree->Branch("ele_EoverP",&ele_EoverP_);
  tree->Branch("ele_Xposition",&ele_Xposition_);
  tree->Branch("ele_Yposition",&ele_Yposition_);
  tree->Branch("ele_dr03TkSumPt",&ele_dr03TkSumPt_);
  tree->Branch("ele_hcalDepth1OverEcal",&ele_hcalDepth1OverEcal_);
  tree->Branch("ele_hcalDepth2OverEcal",&ele_hcalDepth2OverEcal_);
  tree->Branch("ele_dr03HcalDepth2TowerSumEt",&ele_dr03HcalDepth2TowerSumEt_);
  tree->Branch("ele_hcalDepth2TowerSumEtNoVeto",&ele_hcalDepth2TowerSumEtNoVeto_);
  tree->Branch("ele_hcalDepth1TowerSumEtNoVeto",&ele_hcalDepth1TowerSumEtNoVeto_);
  tree->Branch("ele_EcalPlusHcald1iso",&ele_EcalPlusHcald1iso_);
  tree->Branch("ele_dr03EcalRecHitSumEt",&ele_dr03EcalRecHitSumEt_);
  tree->Branch("ele_dr03HcalDepth1TowerSumEt",&ele_dr03HcalDepth1TowerSumEt_);
  tree->Branch("ele_dr03HcalDepth1TowerSumEtBc",&ele_dr03HcalDepth1TowerSumEtBc_);
  tree->Branch("ele_pfSumPhotonEt",&ele_pfSumPhotonEt_);
  tree->Branch("ele_pfSumChargedHadronPt",&ele_pfSumChargedHadronPt_);
  tree->Branch("ele_pfSumNeutralHadronEt",&ele_pfSumNeutralHadronEt_);
  tree->Branch("ele_pfSumPUPt",&ele_pfSumPUPt_);
  tree->Branch("ele_pfDeltaBeta",&ele_pfDeltaBeta_);
  tree->Branch("ele_FirstGenMatch",&ele_FirstGenMatch_);
  tree->Branch("ele_SecondGenMatch", &ele_SecondGenMatch_);
  tree->Branch("ele_DecayChain", &ele_GenMatchTest_);
  tree->Branch("ele_isEB" ,&ele_isEB_);
  tree->Branch("ele_isEE" ,&ele_isEE_);
  tree->Branch("ele_eSuperClusterOverP" ,&ele_eSuperClusterOverP_);
  tree->Branch("ele_ecalEnergy" ,&ele_ecalEnergy_);
  tree->Branch("ele_dEtaInSeed" ,&ele_dEtaInSeed_);
  tree->Branch("ele_InvMinusPInv" ,&ele_InvMinusPInv_);
  tree->Branch("ele_PtCorr",&ele_PtCorr_ );
  tree->Branch("ele_PtScaleUp",& ele_PtScaleUp_);
  tree->Branch("ele_PtScaleDown",& ele_PtScaleDown_);
  tree->Branch("ele_PtResUp",& ele_PtResUp_);
  tree->Branch("ele_PtResDown",& ele_PtResDown_);
  tree->Branch("ele_ECorr",& ele_ECorr_);
  tree->Branch("ele_EScaleUp",& ele_EScaleUp_);
  tree->Branch("ele_EScaleDown",& ele_EScaleDown_);
  tree->Branch("ele_EResUp",& ele_EResUp_);
  tree->Branch("ele_EResDown",& ele_EResDown_);
  tree->Branch("ele_IsoEffArea",& ele_IsoEffArea_);

}

void BigNtuple::fill_eleInfo(const pat::Electron& ele_, const reco::Vertex& pv, double Rho,  double match1, double match2,  std::auto_ptr<EcalClusterLazyTools> recHitEcal, double Iso, double match3){

  float dEtaInSeed;

  if(ele_.superCluster().isNonnull() and ele_.superCluster()->seed().isNonnull())
    dEtaInSeed = ele_.deltaEtaSuperClusterTrackAtVtx() - ele_.superCluster()->eta() + ele_.superCluster()->seed()->eta();
  
  else dEtaInSeed =  std::numeric_limits<float>::max();

  ele_FirstGenMatch_.push_back(match1);
  ele_SecondGenMatch_.push_back(match2);
  ele_GenMatchTest_.push_back(match3);

  ele_Et_.push_back(ele_.superCluster()->energy() * sin(ele_.p4().theta()));
  ele_EtFromCaloEn_.push_back(ele_.caloEnergy() * sin(ele_.p4().theta()));
  
  ele_pt_.push_back(ele_.pt());
  ele_px_.push_back(ele_.px());
  ele_py_.push_back(ele_.py());
  ele_pz_.push_back(ele_.pz());
  ele_etaSC_.push_back(ele_.superCluster()->eta());    //eta SC
  ele_phiSC_.push_back(ele_.superCluster()->phi());    //phi SC
  ele_phiWidth_.push_back(ele_.superCluster()->phiWidth());
  ele_etaWidth_.push_back(ele_.superCluster()->etaWidth());
  ele_energySC_.push_back(ele_.superCluster()->energy()); //energy SC
  ele_thetaSC_.push_back(ele_.caloPosition().theta()); //theta SC
  ele_preshowerEnergySC_.push_back(ele_.superCluster()->preshowerEnergy());

  ele_etaTrack_.push_back(ele_.p4().eta());     //eta track
  ele_phiTrack_.push_back(ele_.p4().phi());     //phi track
  ele_thetaTrack_.push_back(ele_.p4().theta()); //theta track

  ele_eta_.push_back(ele_.eta());
  ele_phi_.push_back(ele_.phi());
  ele_theta_.push_back(ele_.theta());
  ele_energy_.push_back(ele_.energy());

  ele_isEB_.push_back(ele_.isEB());
  ele_isEE_.push_back(ele_.isEE());
  ele_eSuperClusterOverP_.push_back(ele_.eSuperClusterOverP());
  ele_ecalEnergy_.push_back(ele_.ecalEnergy());
  ele_dEtaInSeed_.push_back(std::abs(dEtaInSeed));
  ele_InvMinusPInv_.push_back(std::abs(1.0 - ele_.eSuperClusterOverP())/ele_.ecalEnergy());

  ele_x_.push_back(ele_.p4().x());
  ele_y_.push_back(ele_.p4().y());
  ele_z_.push_back(ele_.p4().z());

  ele_e2x5Max_.push_back(ele_.e2x5Max());
  ele_e1x5_.push_back(ele_.e1x5());
  ele_e5x5_.push_back(ele_.e5x5());
  ele_e2x5MaxOver5x5_.push_back(ele_.e2x5Max()/ele_.e5x5());
  ele_e1x5Over5x5_.push_back(ele_.e1x5()/ele_.e5x5());
  ele_sigmaIetaIetaFull5x5_.push_back(ele_.full5x5_sigmaIetaIeta());
  ele_e2x5MaxFull5x5_.push_back(ele_.full5x5_e2x5Max());
  ele_e1x5Full5x5_.push_back(ele_.full5x5_e1x5());
  ele_e5x5Full5x5_.push_back(ele_.full5x5_e5x5());
  ele_e2x5MaxOver5x5Full5x5_.push_back(ele_.full5x5_e2x5Max()/ele_.full5x5_e5x5());
  ele_e1x5Over5x5Full5x5_.push_back(ele_.full5x5_e1x5()/ele_.full5x5_e5x5());

  ele_zTrackPositionAtVtx_.push_back(ele_.TrackPositionAtVtx().Z());
  ele_hadronicOverEm_.push_back(ele_.hadronicOverEm());
  ele_deltaEtaInSC_.push_back(ele_.deltaEtaSuperClusterTrackAtVtx());
  ele_deltaPhiInSC_.push_back(ele_.deltaPhiSuperClusterTrackAtVtx());
  ele_deltaEtaInSeedCluster_.push_back(ele_.deltaEtaSeedClusterTrackAtVtx());
  ele_deltaPhiInSeedCluster_.push_back(ele_.deltaPhiSeedClusterTrackAtCalo());
  ele_sigmaIetaIeta_.push_back(ele_.sigmaIetaIeta());

  EBDetId BarrelId = ele_.superCluster()->seed()->seed();

  ele_rawId_.push_back(BarrelId.rawId());
  ele_ieta_.push_back(BarrelId.ieta());

  ele_e2x5Right_.push_back(recHitEcal->e2x5Right(*(ele_.superCluster()->seed())));
  ele_e2x5Left_.push_back(recHitEcal->e2x5Left(*(ele_.superCluster()->seed())));
  ele_e2x5Top_.push_back(recHitEcal->e2x5Top(*(ele_.superCluster()->seed())));
  ele_e2x5Bottom_.push_back(recHitEcal->e2x5Bottom(*(ele_.superCluster()->seed())));
  ele_eMax_.push_back(recHitEcal->eMax(*(ele_.superCluster()->seed())));
  ele_eRight_.push_back(recHitEcal->eRight(*(ele_.superCluster()->seed())));
  ele_eLeft_.push_back(recHitEcal->eLeft(*(ele_.superCluster()->seed())));
  ele_eTop_.push_back(recHitEcal->eTop(*(ele_.superCluster()->seed())));
  ele_eBottom_.push_back(recHitEcal->eBottom(*(ele_.superCluster()->seed())));
  ele_e3x3_.push_back(recHitEcal->e3x3(*(ele_.superCluster()->seed())));
  ele_frac51_.push_back( recHitEcal->e5x1(*(ele_.superCluster()->seed()))/ele_.full5x5_e5x5() );
  ele_frac15_.push_back( recHitEcal->e1x5(*(ele_.superCluster()->seed()))/ele_.full5x5_e5x5() );

  ele_dxy_.push_back(ele_.gsfTrack()->dxy(pv.position()));   //GSF -> Gaussian Sum Filter
  ele_dz_.push_back(ele_.gsfTrack()->dz(pv.position()));

  ele_3dIP.push_back(ele_.dB(pat::Electron::PV3D));
  ele_3dIPSig.push_back(ele_.dB(pat::Electron::PV3D)/ele_.edB(pat::Electron::PV3D));
  ele_2dIP.push_back(ele_.dB());
  ele_2dIPSig.push_back(ele_.dB()/ele_.edB());

  ele_isEcalDrivenSeed_.push_back(ele_.ecalDrivenSeed());
  ele_isPassConversionVeto_.push_back(ele_.passConversionVeto());
  ele_charge_.push_back(ele_.gsfTrack()->charge());
  ele_rhoIso_.push_back(Rho); //transverse momentum per unit area
  ele_nbOfMissingHits_.push_back(ele_.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
  ele_fbrem_.push_back(ele_.fbrem());
  ele_EoverP_.push_back(ele_.eSeedClusterOverP());
  ele_Xposition_.push_back(ele_.caloPosition().x());
  ele_Yposition_.push_back(ele_.caloPosition().y());

  // electron correction
  ele_PtCorr_.push_back(ele_.pt()*ele_.userFloat("ecalTrkEnergyPostCorr")/ele_.energy());
  ele_PtScaleUp_.push_back(ele_.pt()*ele_.userFloat("energyScaleUp")/ele_.energy());
  ele_PtScaleDown_.push_back(ele_.pt()*ele_.userFloat("energyScaleDown")/ele_.energy());
  ele_PtResUp_.push_back(ele_.pt()*ele_.userFloat("energySigmaUp")/ele_.energy());
  ele_PtResDown_.push_back(ele_.pt()*ele_.userFloat("energySigmaDown")/ele_.energy());
  ele_ECorr_.push_back(ele_.userFloat("ecalTrkEnergyPostCorr"));
  ele_EScaleUp_.push_back(ele_.userFloat("energyScaleUp"));
  ele_EScaleDown_.push_back(ele_.userFloat("energyScaleDown"));
  ele_EResUp_.push_back(ele_.userFloat("energySigmaUp"));
  ele_EResDown_.push_back(ele_.userFloat("energySigmaDown"));

    //tracker isolation
  ele_dr03TkSumPt_.push_back(ele_.dr03TkSumPt());

    //------------- detector isolation -------------------------
  ele_hcalDepth1OverEcal_.push_back(ele_.hcalDepth1OverEcal());
  ele_hcalDepth2OverEcal_.push_back(ele_.hcalDepth2OverEcal());
  ele_dr03HcalDepth2TowerSumEt_.push_back(ele_.dr03HcalDepth2TowerSumEt());
  ele_hcalDepth2TowerSumEtNoVeto_.push_back(ele_.isolationVariables03().hcalDepth2TowerSumEt);// hcaldepht2 iso deposit with
  // electron footprint removed
  ele_hcalDepth1TowerSumEtNoVeto_.push_back(ele_.isolationVariables03().hcalDepth1TowerSumEt);// hcaldepht1 iso deposit with
  // electron footprint removed
  ele_EcalPlusHcald1iso_.push_back(ele_.dr03EcalRecHitSumEt() + ele_.dr03HcalDepth1TowerSumEt());
  ele_dr03EcalRecHitSumEt_.push_back(ele_.dr03EcalRecHitSumEt());
  ele_dr03HcalDepth1TowerSumEt_.push_back(ele_.dr03HcalDepth1TowerSumEt());
  ele_dr03HcalDepth1TowerSumEtBc_.push_back(ele_.dr03HcalDepth1TowerSumEtBc());
  //------------- PF isolation from pat::electron -------------------------
  ele_pfSumPhotonEt_.push_back(ele_.pfIsolationVariables().sumPhotonEt);
  ele_pfSumChargedHadronPt_.push_back(ele_.pfIsolationVariables().sumChargedHadronPt);
  ele_pfSumNeutralHadronEt_.push_back(ele_.pfIsolationVariables().sumNeutralHadronEt);
  ele_pfSumPUPt_.push_back(ele_.pfIsolationVariables().sumPUPt);
  // deltaBeta
  double charged   = ele_.pfIsolationVariables().sumPhotonEt;
  double neutral   = ele_.pfIsolationVariables().sumNeutralHadronEt;
  double pileup    = ele_.pfIsolationVariables().sumPUPt;
  double deltaBeta = charged + std::max(0.0, neutral-0.5*pileup)/ele_.pt();
  ele_pfDeltaBeta_.push_back(deltaBeta);
  ele_IsoEffArea_.push_back(Iso);

}
void BigNtuple::set_eleIDInfo(TTree* tree){

  tree->Branch("ele_Mva2016" , &ele_Mva2016_);
  tree->Branch("ele_CutVeto" , &ele_CutVeto_);
  tree->Branch("ele_CutLoose" , &ele_CutLoose_);
  tree->Branch("ele_CutMedium" , &ele_CutMedium_);
  tree->Branch("ele_CutTight" , &ele_CutTight_);

}

void BigNtuple::fill_eleIDInfo(float ele_mva , bool ele_veto , bool ele_loose , bool ele_medium ,bool ele_tight){

  ele_Mva2016_.push_back(ele_mva);
  ele_CutVeto_.push_back(ele_veto);
  ele_CutLoose_.push_back(ele_loose);
  ele_CutMedium_.push_back(ele_medium);
  ele_CutTight_.push_back(ele_tight);

}

void BigNtuple::set_metInfo(TTree* tree){

  tree->Branch("pfMet_et" , &pfMet_et_, "pfMet_et/F");
  tree->Branch("pfMet_pt" , &pfMet_pt_, "pfMet_pt/F");
  tree->Branch("pfMet_phi" ,&pfMet_phi_,"pfMet_phi/F");
  tree->Branch("pfMet_en" , &pfMet_en_, "pfMet_en/F");
  tree->Branch("pfMet_px" , &pfMet_px_, "pfMet_px/F");
  tree->Branch("pfMet_py" , &pfMet_py_, "pfMet_py/F");
  tree->Branch("pfMet_pz" , &pfMet_pz_, "pfMet_pz/F");
  tree->Branch("pfMet_sumEt" , &pfMet_sumEt_ ,"pfMet_sumEt/F");
  tree->Branch("caloMet_pt"  , &caloMet_pt_  ,"caloMet_pt/F");
  tree->Branch("caloMet_phi" , &caloMet_phi_ ,"caloMet_phi/F");
  tree->Branch("metJECDown"  , &metJECDown_  ,"metJECDown/F");
  tree->Branch("metJECUp"    , &metJECUp_    ,"metJECUp/F");
  tree->Branch("metUnclDown" , &metUnclDown_ ,"metUnclDow/F");
  tree->Branch("metUnclUp"   , &metUnclUp_   ,"metUnclUp/F");
  tree->Branch("metPhiJECDown" , &metPhiJECDown_ ,"metPhiJECDown/F");
  tree->Branch("metPhiUnclUp"  , &metPhiUnclUp_  ,"metPhiUnclUp/F");
  tree->Branch("metPhiUnclDown", &metPhiUnclDown_,"metPhiUnclDown/F");
  tree->Branch("pfmet_Rawpt"   , &pfmet_Rawpt_   ,"pfmet_Rawpt/F");
  tree->Branch("pfmet_RawPhi"  , &pfmet_RawPhi_  ,"pfmet_RawPhi/F");
}



void BigNtuple::fill_metInfo(const pat::MET& met){

  pfmet_Rawpt_        = met.uncorPt();
  pfmet_RawPhi_       = met.uncorPhi();


  pfMet_et_  = met.et();
  pfMet_pt_  = met.pt();
  pfMet_phi_ = met.phi();
  pfMet_en_  = met.energy();
  pfMet_px_  = met.px();
  pfMet_py_  = met.py();
  pfMet_pz_  = met.pz();
  pfMet_sumEt_ = met.sumEt();

  caloMet_pt_  = met.caloMETPt();
  caloMet_phi_ = met.caloMETPhi();

  metJECDown_      = met.shiftedPt(pat::MET::JetEnDown);
  metJECUp_        = met.shiftedPt(pat::MET::JetEnUp);
  metUnclDown_     = met.shiftedPt(pat::MET::UnclusteredEnDown);
  metUnclUp_       = met.shiftedPt(pat::MET::UnclusteredEnUp);
  metPhiJECDown_   = met.shiftedPhi(pat::MET::JetEnDown);
  metPhiJECUp_     = met.shiftedPhi(pat::MET::JetEnUp);
  metPhiUnclUp_    = met.shiftedPhi(pat::MET::UnclusteredEnUp);
  metPhiUnclDown_  = met.shiftedPhi(pat::MET::UnclusteredEnDown);



}

//transverse mass
void BigNtuple::set_transverseMassInfo(TTree* tree){
  tree->Branch("tranvsverseMass_ivf" , &tranvsverseMass_ivf_, "tranvsverseMass_ivf/F");
  tree->Branch("tranvsverseMass_lep1" , &tranvsverseMass_lep1_, "tranvsverseMass_lep1/F");
  tree->Branch("tranvsverseMass_ivfPluslep1" , &tranvsverseMass_ivfPluslep1_, "tranvsverseMass_ivfPluslep1/F");
}

void BigNtuple::fill_transverseMassInfo(float tr_mass_ivf, float tr_mass_lep1, float tr_mass_ivfPluslep1){
  tranvsverseMass_ivf_ = tr_mass_ivf;
  tranvsverseMass_lep1_ = tr_mass_lep1;
  tranvsverseMass_ivfPluslep1_ = tr_mass_ivfPluslep1;
}

void BigNtuple::set_massCorrection(TTree *tree){
  //tree->Branch("sv_rho", sv_rho_, "sv_rho/F");
  //tree->Branch("sv_phi", _sv_phi_, "sv_phi/F");
  //tree->Branch("sv_theta", sv_theta_, "sv_theta/F");
  tree->Branch("sv_mass_corr", &sv_mass_corr_, "sv_mass_corr/F");
}

void BigNtuple::fill_massCorrection(double mass_corr){
  sv_mass_corr_ = mass_corr;
}


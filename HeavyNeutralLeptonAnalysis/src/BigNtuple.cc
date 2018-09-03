#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"
#include "TVector3.h"

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
       tree->Branch("pvLxySig" , &pvLxySig_, "pvLxySig/F");
       tree->Branch("pvLxyzSig" , &pvLxyzSig_, "pvLxyzSig/F");
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
       pvLxySig_ = std::sqrt( x * x + y * y ) / std::sqrt(xE * xE + yE * yE);
       pvLxyzSig_ = std::sqrt( x * x + y * y + z * z )/ std::sqrt(xE * xE + yE * yE + zE * zE);
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
}


void BigNtuple::set_trigInfo(TTree* tree){
  
    tree->Branch("passIsoTk18" , &passIsoTk18_, "passIsoTk18/O");
    tree->Branch("passIsoTk20" , &passIsoTk20_, "passIsoTk20/O");
    tree->Branch("passIsoTk22" , &passIsoTk22_, "passIsoTk22/O");
    tree->Branch("passIsoTk24" , &passIsoTk24_, "passIsoTk24/O");
    tree->Branch("passIsoTk27" , &passIsoTk27_, "passIsoTk27/O");
    tree->Branch("passIsoTk17e" , &passIsoTk17e_, "passIsoTk17e/O");
    tree->Branch("passIsoTk22e" , &passIsoTk22e_, "passIsoTk22e/O");
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


}

void BigNtuple::fill_trigInfo(const edm::TriggerResults& triggerResults, const edm::TriggerNames& trigNames){

  for (size_t i = 0; i < trigNames.size(); ++i) {
    const std::string &name = trigNames.triggerName(i);
    bool fired = triggerResults.accept(i);
    if(!fired) continue;

    passIsoTk18_  |=  name.find("HLT_IsoTkMu18_v") != std::string::npos;
    passIsoTk20_  |=  name.find("HLT_IsoTkMu20_v") != std::string::npos;
    passIsoTk22_  |=  name.find("HLT_IsoTkMu22_v") != std::string::npos;
    passIsoTk24_  |=  name.find("HLT_IsoTkMu24_v") != std::string::npos;
    passIsoTk27_  |=  name.find("HLT_IsoTkMu27_v") != std::string::npos;
    passIsoTk17e_  |=  name.find("HLT_IsoTkMu17_eta2p1_v") != std::string::npos;
    passIsoTk22e_  |=  name.find("HLT_IsoTkMu22_eta2p1_v") != std::string::npos;
    passIsoMu18_  |=  name.find("HLT_IsoMu18_v") != std::string::npos;
    passIsoMu20_  |=  name.find("HLT_IsoMu20_v") != std::string::npos;
    passIsoMu22_  |=  name.find("HLT_IsoMu22_v") != std::string::npos;
    passIsoMu24_  |=  name.find("HLT_IsoMu24_v") != std::string::npos;
    passIsoMu27_  |=  name.find("HLT_IsoMu27_v") != std::string::npos;

    passIsoMu17e_ |=  name.find("HLT_IsoTkMu17_eta2p1_v") != std::string::npos;
    passIsoMu22e_ |=  name.find("HLT_IsoTkMu22_eta2p1_v") != std::string::npos;

    passTkMu17_   |=  name.find("HLT_TkMu17_v") != std::string::npos;
    passTkMu20_   |=  name.find("HLT_TkMu20_v") != std::string::npos;

    passDoubleMu17TrkIsoMu8_  |=  name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") != std::string::npos;
    passDoubleMu17TrkIsoTkMu8_ |=  name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos;
    passDoubleTkMu17TrkIsoTkMu8_ |=  name.find("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos;

  }

  passIsoMu24All_ = passIsoMu24All_   || passIsoMu24_ || passIsoTk24_ ;
  passIsoMu27All_ = passIsoMu27All_   || passIsoMu27_ || passIsoTk27_ ;

}

  

void BigNtuple::set_svInfo(TTree* tree){

    tree->Branch("sv_TrackSize" , &sv_TrackSize_);
    tree->Branch("sv_LXYSig" , &sv_LXYSig_);
    tree->Branch("sv_LXYZSig" , &sv_LXYZSig_);
    tree->Branch("sv_LXY" , &sv_LXY_);
    tree->Branch("sv_LXYZ" , &sv_LXYZ_);
    tree->Branch("sv_mass" , &sv_mass_);
    tree->Branch("sv_charge" , &sv_charge_);
    tree->Branch("sv_eta" , &sv_eta_);
    tree->Branch("sv_phi" , &sv_phi_);
    tree->Branch("sv_pt" , &sv_pt_);
    tree->Branch("sv_p" , &sv_p_);
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
    tree->Branch("sv_tracks_dxySig" , &sv_tracks_dxySig_);
    tree->Branch("sv_tracks_dxy" , &sv_tracks_dxy_);
    tree->Branch("sv_tracks_dxyz" , &sv_tracks_dxyz_);
    tree->Branch("sv_tracks_Sumcharge" , &sv_tracks_Sumcharge_);
    tree->Branch("sv_tracks_Sumpt" , &sv_tracks_Sumpt_);
    tree->Branch("sv_match" , &sv_match_);
}


void BigNtuple::fill_svInfo(const reco::Vertex& bestVertex, const reco::Vertex& pv , double match){

  float  svChi2 = bestVertex.chi2();
  float  svNDof = bestVertex.ndof();

  //flight distance from the firstPV                                                                                                                                  
  float x  = bestVertex.x(), y = bestVertex.y(), z = bestVertex.z();
  float dx = x - pv.x() , dy = y - pv.y(), dz = z - pv.z();

  //build the total error                                                                                                                                                  
  float svxE = bestVertex.xError(), svyE = bestVertex.yError(), svzE = bestVertex.zError();
  float pvxE = pv.xError(), pvyE = pv.yError(), pvzE = pv.zError();
  float xE   = std::sqrt(svxE * svxE + pvxE * pvxE), yE = std::sqrt(svyE * svyE + pvyE * pvyE), zE = std::sqrt(svzE * svzE + pvzE * pvzE);

  // mother beta, gamma, ctau                                                                                                                                              
  float   beta_mom  = bestVertex.p4().P() / bestVertex.p4().energy();
  float   gamma_mom = bestVertex.p4().energy() / bestVertex.p4().mass();

  TVector3 pvVector3D(pv.x(), pv.y(), pv.z());
  TVector3 pvVector2D(pv.x(), pv.y(), 0);
  TVector3 svVector3D(x, y, z);
  TVector3 svVector2D(x, y, 0);

  // line pointing form the primary vertex through the sceondary vertex                                                                                                 
                
  TVector3 svMom3D( bestVertex.p4().x(), bestVertex.p4().y(), bestVertex.p4().z());
  TVector3 svMom2D( bestVertex.p4().x(), bestVertex.p4().y(), 0);


  // you want the negative when the momentum and sv are in the same                                                                                                     
  // direction relative to the PV                                                                                                                                       
  // this makes sure the angle is not pi when the vertex is fit behind                                                                                                  
  // the primary vertex                                                                                                                                                  
                
  float sign2D =  (svMom2D * (svVector2D - pvVector2D)) > 0 ? -1: 1;
  float sign3D =  (svMom3D * (svVector3D - pvVector3D)) > 0 ? -1: 1;

  TVector3 pvToVertex3D( sign3D * dx, sign3D * dy, sign3D * dz);
  TVector3 pvToVertex2D( sign2D * dx, sign2D * dy, 0);

  float  svAngle3D = pvToVertex3D.Angle(svMom3D);
  float  svAngle2D = pvToVertex2D.Angle(svMom2D);

  sv_TrackSize_.push_back(bestVertex.nTracks());
  sv_LXY_.push_back(std::sqrt( dx * dx + dy * dy ));
  sv_LXYZ_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz ));
  sv_LXYSig_.push_back(std::sqrt( dx * dx + dy * dy ) / std::sqrt(xE * xE + yE * yE));
  sv_LXYZSig_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz) / std::sqrt(xE * xE + yE * yE + zE * zE));
  sv_mass_.push_back(bestVertex.p4().mass());
  sv_eta_.push_back(bestVertex.p4().eta());
  sv_phi_.push_back(bestVertex.p4().phi());
  sv_pt_.push_back(bestVertex.p4().pt());
  sv_p_.push_back(bestVertex.p4().P());
  sv_Beta_.push_back(beta_mom);
  sv_Gamma_.push_back(gamma_mom);
  sv_CTau0_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz) / (beta_mom * gamma_mom));
  sv_NDof_.push_back(svNDof);
  sv_Chi2_.push_back(svChi2);
  sv_Angle3D_.push_back(svAngle3D);
  sv_Angle2D_.push_back(svAngle2D);

  int ch = 0;
  float pt = 0;

  // What is emplace_back(); ?
  sv_tracks_charge_.emplace_back();
  sv_tracks_eta_.emplace_back();
  sv_tracks_phi_.emplace_back();
  sv_tracks_pt_.emplace_back();
  sv_tracks_dxySig_.emplace_back();
  sv_tracks_dxy_.emplace_back();
  sv_tracks_dxyz_.emplace_back();

  reco::Vertex::trackRef_iterator tt = bestVertex.tracks_begin();
  for(; tt != bestVertex.tracks_end(); ++tt) {
    
    sv_tracks_charge_.back().push_back((*tt)->charge());
    sv_tracks_eta_.back().push_back((*tt)->eta());
    sv_tracks_phi_.back().push_back((*tt)->phi());
    sv_tracks_pt_.back().push_back((*tt)->pt());
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
  sv_match_.push_back(match);
}


 void BigNtuple::set_muInfo(TTree* tree){
   tree->Branch("mu_en" , &mu_en_);
   tree->Branch("mu_pt" , &mu_pt_);
   tree->Branch("mu_eta" , &mu_eta_);
   tree->Branch("mu_phi" , &mu_phi_);
   tree->Branch("mu_et" , &mu_et_);
   tree->Branch("mu_charge" , &mu_charge_);
   tree->Branch("mu_trackiso" , &mu_trackiso_);
   tree->Branch("mu_pfSumChargedHadronPt" , &mu_pfSumChargedHadronPt_);
   tree->Branch("mu_pfSumNeutralHadronEt" , &mu_pfSumNeutralHadronEt_);
   tree->Branch("mu_PFSumPhotonEt" , &mu_PFSumPhotonEt_);
   tree->Branch("mu_pfSumPUPt" , &mu_pfSumPUPt_);
   tree->Branch("mu_numberOfValidMuonHits" , &mu_numberOfValidMuonHits_);
   tree->Branch("mu_emIso" , &mu_emIso_);
   tree->Branch("mu_hadIso" , &mu_hadIso_);
   tree->Branch("mu_normalizedChi2" , &mu_normalizedChi2_);
   tree->Branch("mu_numberOfMatchedStations" , &mu_numberOfMatchedStations_);
   tree->Branch("mu_numberOfValidPixelHits" , &mu_numberOfValidPixelHits_);
   tree->Branch("mu_numberOftrackerLayersWithMeasurement" , &mu_numberOftrackerLayersWithMeasurement_);
   tree->Branch("mu_numberOfpixelLayersWithMeasurement" , &mu_numberOfpixelLayersWithMeasurement_);
   tree->Branch("mu_TrackQuality" , &mu_TrackQuality_);
   tree->Branch("mu_InnerTrackQuality" , &mu_InnerTrackQuality_);
   tree->Branch("mu_pxTunePMuonBestTrack" , &mu_pxTunePMuonBestTrack_);
   tree->Branch("mu_pyTunePMuonBestTrack" , &mu_pyTunePMuonBestTrack_);
   tree->Branch("mu_pzTunePMuonBestTrack" , &mu_pzTunePMuonBestTrack_);
   tree->Branch("mu_pTunePMuonBestTrack" , &mu_pTunePMuonBestTrack_);
   tree->Branch("mu_etaTunePMuonBestTrack" , &mu_etaTunePMuonBestTrack_);
   tree->Branch("mu_LXYZ" , &mu_LXYZ_);
   tree->Branch("mu_LXY" , &mu_LXY_);
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
   tree->Branch("mu_recoDeltaBeta" , &mu_recoDeltaBeta_);
   tree->Branch("mu_recoiso" , &mu_recoiso_);
   tree->Branch("mu_isGlobalMuon" , &mu_isGlobalMuon_);
   tree->Branch("mu_isStandAloneMuon" , &mu_isStandAloneMuon_);
   tree->Branch("mu_isPF" , &mu_isPF_);
   tree->Branch("mu_isRPCMuon" , &mu_isRPCMuon_);
   tree->Branch("mu_isTrackerMuon" , &mu_isTrackerMuon_);
   tree->Branch("mu_isGoodMuon" , &mu_isGoodMuon_);
   tree->Branch("mu_isSoftMuon" , &mu_isSoftMuon_);
   tree->Branch("mu_isLoose" , &mu_isLoose_);
   tree->Branch("mu_isTightMuon" , &mu_isTightMuon_);
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
   tree->Branch("mu_STATofDirection" , &mu_STATofDirection_);
   tree->Branch("mu_STATofNDof" , &mu_STATofNDof_);
   tree->Branch("mu_STATofTimeAtIpInOut" , &mu_STATofTimeAtIpInOut_);
   tree->Branch("mu_STATofTimeAtIpInOutErr" , &mu_STATofTimeAtIpInOutErr_);
   tree->Branch("mu_STATofTimeAtIpOutIn" , &mu_STATofTimeAtIpOutIn_);
   tree->Branch("mu_STATofTimeAtIpOutInErr" , &mu_STATofTimeAtIpOutInErr_);
   tree->Branch("mu_FirstGenMatch" , &mu_FirstGenMatch_);
   tree->Branch("mu_SecondGenMatch" , &mu_SecondGenMatch_);

 }



void BigNtuple::fill_muInfo(const pat::Muon& mu, const reco::Vertex& pv , double match1 , double match2){

  mu_isGlobalMuon_.push_back(mu.isGlobalMuon());
  mu_isPF_.push_back(mu.isPFMuon());
  mu_isTrackerMuon_.push_back(mu.isTrackerMuon());
  mu_isRPCMuon_.push_back(mu.isRPCMuon());
  mu_isStandAloneMuon_.push_back(mu.isStandAloneMuon());
  mu_isSoftMuon_.push_back(mu.isSoftMuon(pv));
  mu_isLoose_.push_back(mu.isLooseMuon());
  mu_isTightMuon_.push_back(mu.isTightMuon(pv));
  mu_en_.push_back(mu.energy());
  mu_et_.push_back(mu.et());
  mu_pt_.push_back(mu.pt());
  mu_eta_.push_back(mu.eta());
  mu_phi_.push_back(mu.phi());
  mu_charge_.push_back(mu.charge());
  mu_FirstGenMatch_.push_back(match1);
  mu_SecondGenMatch_.push_back(match2);

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

  if(mu.globalTrack().isNonnull() ) {
    mu_normalizedChi2_.push_back(mu.globalTrack()->normalizedChi2());
    mu_numberOfValidPixelHits_.push_back(mu.innerTrack()->hitPattern().numberOfValidPixelHits());
    mu_numberOfValidMuonHits_.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
    mu_numberOftrackerLayersWithMeasurement_.push_back(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
    mu_numberOfMatchedStations_.push_back(mu.numberOfMatchedStations());
    mu_numberOfpixelLayersWithMeasurement_.push_back(mu.innerTrack()->hitPattern().pixelLayersWithMeasurement());
    mu_InnerTrackQuality_.push_back(mu.innerTrack()->quality(reco::TrackBase::highPurity));
  }

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

  reco::MuonTime tofAll = mu.time();
  mu_STATofDirection_.push_back(tofAll.direction());
  mu_STATofNDof_.push_back(tofAll.nDof);
  mu_STATofTimeAtIpInOut_.push_back(tofAll.timeAtIpInOut);
  mu_STATofTimeAtIpInOutErr_.push_back(tofAll.timeAtIpInOutErr);
  mu_STATofTimeAtIpOutIn_.push_back(tofAll.timeAtIpOutIn);
  mu_STATofTimeAtIpOutInErr_.push_back(tofAll.timeAtIpOutInErr);
}



void BigNtuple::set_jetInfo(TTree* tree){

  tree->Branch("jet_charge" , &jet_charge_);
  tree->Branch("jet_et" , &jet_et_);
  tree->Branch("jet_pt" , &jet_pt_);
  tree->Branch("jet_eta" , &jet_eta_);
  tree->Branch("jet_phi" , &jet_phi_);
  tree->Branch("jet_theta" , &jet_theta_);
  tree->Branch("jet_en" , &jet_en_);
  tree->Branch("jet_chargedEmEnergy" , &jet_chargedEmEnergy_);
  tree->Branch("jet_neutralEmEnergyFraction" , &jet_neutralEmEnergyFraction_);
  tree->Branch("jet_chargedHadronEnergy" , &jet_chargedHadronEnergy_);
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

}



void BigNtuple::fill_jetInfo(const pat::Jet& jet){

  jet_charge_.push_back(jet.charge());
  jet_et_.push_back(jet.et());
  jet_pt_.push_back(jet.pt());
  jet_eta_.push_back(jet.eta());
  jet_phi_.push_back(jet.phi());
  jet_theta_.push_back(jet.theta());
  jet_en_.push_back(jet.energy());
  jet_chargedEmEnergy_.push_back(jet.chargedEmEnergy());
  jet_neutralEmEnergyFraction_.push_back(jet.neutralEmEnergyFraction());
  jet_chargedHadronEnergy_.push_back(jet.chargedHadronEnergy());
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

}

void BigNtuple::set_eleInfo(TTree* tree){
  /*
  tree->Branch("ele_Mva",&ele_Mva_);
  tree->Branch("ele_MvaFall17Iso",&ele_MvaFall17Iso_);
  tree->Branch("ele_MvaFall17NoIso",&ele_MvaFall17NoIso_);
  tree->Branch("ele_CutBasedVeto",&ele_CutBasedVeto_);
  tree->Branch("ele_CutBasedLoose",&ele_CutBasedLoose_);
  tree->Branch("ele_CutBasedMedium",&ele_CutBasedMedium_);
  tree->Branch("ele_CutBasedTight",&ele_CutBasedTight_);
  */
  tree->Branch("ele_Et",&ele_Et_);
  tree->Branch("ele_EtFromCaloEn",&ele_EtFromCaloEn_);    
  tree->Branch("ele_pt",&ele_pt_); 
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
}

void BigNtuple::fill_eleInfo(const pat::Electron& ele_, const reco::Vertex& pv, double Rho, double match1, double match2,  std::auto_ptr<EcalClusterLazyTools> recHitEcal){


  ele_FirstGenMatch_.push_back(match1);
  ele_SecondGenMatch_.push_back(match2);
  ele_Et_.push_back(ele_.superCluster()->energy() * sin(ele_.p4().theta()));
  ele_EtFromCaloEn_.push_back(ele_.caloEnergy() * sin(ele_.p4().theta()));
    
  ele_pt_.push_back(ele_.pt()); 
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

  ele_isEcalDrivenSeed_.push_back(ele_.ecalDrivenSeed());
  ele_isPassConversionVeto_.push_back(ele_.passConversionVeto());
  ele_charge_.push_back(ele_.gsfTrack()->charge());
  ele_rhoIso_.push_back(Rho); //transverse momentum per unit area
  ele_nbOfMissingHits_.push_back(ele_.gsfTrack()->numberOfLostHits()); 
  ele_fbrem_.push_back(ele_.fbrem());
  ele_EoverP_.push_back(ele_.eSeedClusterOverP());
  ele_Xposition_.push_back(ele_.caloPosition().x());   
  ele_Yposition_.push_back(ele_.caloPosition().y()); 

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
  double deltaBeta = charged + std::max(0.0, neutral-0.5*pileup);
  ele_pfDeltaBeta_.push_back(deltaBeta);

}

void BigNtuple::set_metInfo(TTree* tree){

  tree->Branch("pfMet_et" , &pfMet_et_, "pfMet_et/F");
  tree->Branch("pfMet_pt" , &pfMet_pt_, "pfMet_pt/F");
  tree->Branch("pfMet_phi" , &pfMet_phi_, "pfMet_phi/F");
  tree->Branch("pfMet_en" , &pfMet_en_, "pfMet_en/F");
  tree->Branch("pfMet_px" , &pfMet_px_, "pfMet_px/F");
  tree->Branch("pfMet_py" , &pfMet_py_, "pfMet_py/F");
  tree->Branch("pfMet_pz" , &pfMet_pz_, "pfMet_pz/F");
  tree->Branch("pfMet_sumEt" , &pfMet_sumEt_, "pfMet_sumEt/F");
  tree->Branch("caloMet_pt" , &caloMet_pt_, "caloMet_pt/F");
  tree->Branch("caloMet_phi" , &caloMet_phi_, "caloMet_phi/F");   
  
}



void BigNtuple::fill_metInfo(const pat::MET& met){

  pfMet_et_ = met.et();
  pfMet_pt_ = met.pt();
  pfMet_phi_ = met.phi();
  pfMet_en_ = met.energy();
  pfMet_px_ = met.px();
  pfMet_py_ = met.py();
  pfMet_pz_ = met.pz();
  pfMet_sumEt_ = met.sumEt();

  caloMet_pt_ = met.caloMETPt();
  caloMet_phi_ = met.caloMETPhi();


}

void BigNtuple::set_bjetInfo(TTree* tree){
  tree->Branch("jet_btag_pt",&jet_btag_pt);
  tree->Branch("jet_btag_eta",&jet_btag_eta);
  tree->Branch("jet_btag_phi",&jet_btag_phi);
  tree->Branch("jet_btag_flavor",&jet_btag_flavor);
  tree->Branch("jet_btag_pfCSVv2IVF_discriminator",&jet_btag_pfCSVv2IVF_discriminator);
}

void BigNtuple::fill_bjetInfo(const pat::Jet& jet,  const std::string & bDiscr, int flavor){

  jet_btag_pt.push_back(jet.pt());
  jet_btag_eta.push_back(jet.eta());
  jet_btag_phi.push_back(jet.phi());
  jet_btag_flavor.push_back(flavor);
  jet_btag_pfCSVv2IVF_discriminator.push_back(jet.bDiscriminator(bDiscr));

}

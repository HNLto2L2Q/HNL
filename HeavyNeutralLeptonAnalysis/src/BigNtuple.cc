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

void BigNtuple::set_pv_genInfo(TTree* tree) {

  tree->Branch("mu_gen_PID1" , &mu_gen_PID1_);
  tree->Branch("mu_gen_Status1",&mu_gen_Status1_);
  tree->Branch("mu_gen_Charge1",&mu_gen_Charge1_);
  tree->Branch("mu_gen_Pt1",&mu_gen_Pt1_);
  tree->Branch("mu_gen_Eta1",&mu_gen_Eta1_);
  tree->Branch("mu_gen_Phi1",&mu_gen_Phi1_);
  tree->Branch("mu_gen_VX1",&mu_gen_VX1_);
  tree->Branch("mu_gen_VY1",&mu_gen_VY1_);
  tree->Branch("mu_gen_VZ1",&mu_gen_VZ1_);
  tree->Branch("mu_gen_PartVLxy1",&mu_gen_Lxy1_);
  tree->Branch("mu_gen_PartVLxyz1",&mu_gen_Lxyz1_);
  tree->Branch("mu_gen_MomPID1",&mu_gen_MomPID1_);
  tree->Branch("mu_gen_MomStatus1",&mu_gen_MomStatus1_);
  tree->Branch("mu_gen_MomMass1",&mu_gen_MomMass1_);
  tree->Branch("mu_gen_MomCharge1",&mu_gen_MomCharge1_);
  tree->Branch("mu_gen_MomPt1",&mu_gen_MomPt1_);
  tree->Branch("mu_gen_MomEta1",&mu_gen_MomEta1_);
  tree->Branch("mu_gen_MomPhi1",&mu_gen_MomPhi1_);
  tree->Branch("mu_gen_MomBeta1",&mu_gen_MomBeta1_);
  tree->Branch("mu_gen_MomGamma1",&mu_gen_MomGamma1_);
  tree->Branch("mu_gen_MomLxyz1",&mu_gen_MomLxyz1_);
  tree->Branch("mu_gen_MomLz1",&mu_gen_MomLz1_);
  tree->Branch("mu_gen_MomLxy1",&mu_gen_MomLxy1_);

  tree->Branch("ele_gen_PID1" , &ele_gen_PID1_);
  tree->Branch("ele_gen_Status1",&ele_gen_Status1_);
  tree->Branch("ele_gen_Charge1",&ele_gen_Charge1_);
  tree->Branch("ele_gen_Pt1",&ele_gen_Pt1_);
  tree->Branch("ele_gen_Eta1",&ele_gen_Eta1_);
  tree->Branch("ele_gen_Phi1",&ele_gen_Phi1_);
  tree->Branch("ele_gen_VX1",&ele_gen_VX1_);
  tree->Branch("ele_gen_VY1",&ele_gen_VY1_);
  tree->Branch("ele_gen_VZ1",&ele_gen_VZ1_);
  tree->Branch("ele_gen_PartVLxy1",&ele_gen_Lxy1_);
  tree->Branch("ele_gen_PartVLxyz1",&ele_gen_Lxyz1_);
  tree->Branch("ele_gen_MomPID1",&ele_gen_MomPID1_);
  tree->Branch("ele_gen_MomStatus1",&ele_gen_MomStatus1_);
  tree->Branch("ele_gen_MomMass1",&ele_gen_MomMass1_);
  tree->Branch("ele_gen_MomCharge1",&ele_gen_MomCharge1_);
  tree->Branch("ele_gen_MomPt1",&ele_gen_MomPt1_);
  tree->Branch("ele_gen_MomEta1",&ele_gen_MomEta1_);
  tree->Branch("ele_gen_MomPhi1",&ele_gen_MomPhi1_);
  tree->Branch("ele_gen_MomBeta1",&ele_gen_MomBeta1_);
  tree->Branch("ele_gen_MomGamma1",&ele_gen_MomGamma1_);
  tree->Branch("ele_gen_MomLxyz1",&ele_gen_MomLxyz1_);
  tree->Branch("ele_gen_MomLz1",&ele_gen_MomLz1_);
  tree->Branch("ele_gen_MomLxy1",&ele_gen_MomLxy1_);

}
void  BigNtuple::fill_pv_genInfo(const reco::GenParticle prt , const reco::Candidate*  mom){

  float vx = prt.vx(), vy = prt.vy(), vz = prt.vz();
  float mx = mom->vx(), my = mom->vy(), mz = mom->vz();
  float dx = vx - mx, dy = vy - my, dz = vz - mz;

  if(abs(prt.pdgId()) == 13 ){
    mu_gen_PID1_.push_back(prt.pdgId() ); 
    mu_gen_Status1_.push_back(prt.status());
    // genkinematics                                                                                                                          
    mu_gen_Pt1_.push_back(prt.pt());
    mu_gen_Eta1_.push_back(prt.eta());
    mu_gen_Phi1_.push_back(prt.phi());
    mu_gen_Charge1_.push_back(prt.charge());
    // vertexposition                                                                                                                         
    mu_gen_VX1_.push_back(vx);
    mu_gen_VY1_.push_back(vy);
    mu_gen_VZ1_.push_back(vz);
  // vertexflight                                                                                                                           
    mu_gen_Lxy1_.push_back(std::sqrt( vx * vx + vy * vy ));
    mu_gen_Lxyz1_.push_back(std::sqrt( vx * vx + vy * vy + vz * vz));
    // mother beta, gamma,ctau                                                                                                                
    float   beta_mom  = mom->p() / mom->energy();
    float   gamma_mom = mom->energy() / mom->mass();
    // mother quantities related to thedecay                                                                                                  
    // gen id andstatus                                                                                                                       
    mu_gen_MomPID1_.push_back(mom->pdgId());
    mu_gen_MomStatus1_.push_back(mom->status());
    mu_gen_MomCharge1_.push_back(mom->charge());
    mu_gen_MomMass1_.push_back(mom->mass());
    mu_gen_MomPt1_.push_back(mom->pt());
    mu_gen_MomEta1_.push_back(mom->eta());
    mu_gen_MomPhi1_.push_back(mom->phi());
    mu_gen_MomBeta1_.push_back(beta_mom);
    mu_gen_MomGamma1_.push_back(gamma_mom);
    mu_gen_MomLxyz1_.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
    mu_gen_MomLz1_.push_back(dz);
    mu_gen_MomLxy1_.push_back(std::sqrt(dx*dx + dy*dy));
  }

  if(abs(prt.pdgId()) == 11 ){
    ele_gen_PID1_.push_back(prt.pdgId() ); 
    ele_gen_Status1_.push_back(prt.status());
    // genkinematics                                                                                                                          
    ele_gen_Pt1_.push_back(prt.pt());
    ele_gen_Eta1_.push_back(prt.eta());
    ele_gen_Phi1_.push_back(prt.phi());
    ele_gen_Charge1_.push_back(prt.charge());
    // vertexposition                                                                                                                         
    ele_gen_VX1_.push_back(vx);
    ele_gen_VY1_.push_back(vy);
    ele_gen_VZ1_.push_back(vz);
  // vertexflight                                                                                                                           
    ele_gen_Lxy1_.push_back(std::sqrt( vx * vx + vy * vy ));
    ele_gen_Lxyz1_.push_back(std::sqrt( vx * vx + vy * vy + vz * vz));
    // mother beta, gamma,ctau                                                                                                                
    float   beta_mom  = mom->p() / mom->energy();
    float   gamma_mom = mom->energy() / mom->mass();
    // mother quantities related to thedecay                                                                                                  
    // gen id andstatus                                                                                                                       
    ele_gen_MomPID1_.push_back(mom->pdgId());
    ele_gen_MomStatus1_.push_back(mom->status());
    ele_gen_MomCharge1_.push_back(mom->charge());
    ele_gen_MomMass1_.push_back(mom->mass());
    ele_gen_MomPt1_.push_back(mom->pt());
    ele_gen_MomEta1_.push_back(mom->eta());
    ele_gen_MomPhi1_.push_back(mom->phi());
    ele_gen_MomBeta1_.push_back(beta_mom);
    ele_gen_MomGamma1_.push_back(gamma_mom);
    ele_gen_MomLxyz1_.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
    ele_gen_MomLz1_.push_back(dz);
    ele_gen_MomLxy1_.push_back(std::sqrt(dx*dx + dy*dy));
  }
}

void BigNtuple::set_sv_genInfo(TTree* tree) {

  tree->Branch("mu_gen_PID2" , &mu_gen_PID2_);
  tree->Branch("mu_gen_Status2",&mu_gen_Status2_);
  tree->Branch("mu_gen_Charge2",&mu_gen_Charge2_);
  tree->Branch("mu_gen_Pt2",&mu_gen_Pt2_);
  tree->Branch("mu_gen_Eta2",&mu_gen_Eta2_);
  tree->Branch("mu_gen_Phi2",&mu_gen_Phi2_);
  tree->Branch("mu_gen_VX2",&mu_gen_VX2_);
  tree->Branch("mu_gen_VY2",&mu_gen_VY2_);
  tree->Branch("mu_gen_VZ2",&mu_gen_VZ2_);
  tree->Branch("mu_gen_PartVLxy2",&mu_gen_Lxy2_);
  tree->Branch("mu_gen_PartVLxyz2",&mu_gen_Lxyz2_);
  tree->Branch("mu_gen_MomPID2",&mu_gen_MomPID2_);
  tree->Branch("mu_gen_MomStatus2",&mu_gen_MomStatus2_);
  tree->Branch("mu_gen_MomMass2",&mu_gen_MomMass2_);
  tree->Branch("mu_gen_MomCharge2",&mu_gen_MomCharge2_);
  tree->Branch("mu_gen_MomPt2",&mu_gen_MomPt2_);
  tree->Branch("mu_gen_MomEta2",&mu_gen_MomEta2_);
  tree->Branch("mu_gen_MomPhi2",&mu_gen_MomPhi2_);
  tree->Branch("mu_gen_MomBeta2",&mu_gen_MomBeta2_);
  tree->Branch("mu_gen_MomGamma2",&mu_gen_MomGamma2_);
  tree->Branch("mu_gen_MomLxyz2",&mu_gen_MomLxyz2_);
  tree->Branch("mu_gen_MomLz2",&mu_gen_MomLz2_);
  tree->Branch("mu_gen_MomLxy2",&mu_gen_MomLxy2_);
  tree->Branch("mu_gen_MomCTau02" ,&mu_gen_MomCTau02_);

  tree->Branch("ele_gen_PID2" , &ele_gen_PID2_);
  tree->Branch("ele_gen_Status2",&ele_gen_Status2_);
  tree->Branch("ele_gen_Charge2",&ele_gen_Charge2_);
  tree->Branch("ele_gen_Pt2",&ele_gen_Pt2_);
  tree->Branch("ele_gen_Eta2",&ele_gen_Eta2_);
  tree->Branch("ele_gen_Phi2",&ele_gen_Phi2_);
  tree->Branch("ele_gen_VX2",&ele_gen_VX2_);
  tree->Branch("ele_gen_VY2",&ele_gen_VY2_);
  tree->Branch("ele_gen_VZ2",&ele_gen_VZ2_);
  tree->Branch("ele_gen_PartVLxy2",&ele_gen_Lxy2_);
  tree->Branch("ele_gen_PartVLxyz2",&ele_gen_Lxyz2_);
  tree->Branch("ele_gen_MomPID2",&ele_gen_MomPID2_);
  tree->Branch("ele_gen_MomStatus2",&ele_gen_MomStatus2_);
  tree->Branch("ele_gen_MomMass2",&ele_gen_MomMass2_);
  tree->Branch("ele_gen_MomCharge2",&ele_gen_MomCharge2_);
  tree->Branch("ele_gen_MomPt2",&ele_gen_MomPt2_);
  tree->Branch("ele_gen_MomEta2",&ele_gen_MomEta2_);
  tree->Branch("ele_gen_MomPhi2",&ele_gen_MomPhi2_);
  tree->Branch("ele_gen_MomBeta2",&ele_gen_MomBeta2_);
  tree->Branch("ele_gen_MomGamma2",&ele_gen_MomGamma2_);
  tree->Branch("ele_gen_MomLxyz2",&ele_gen_MomLxyz2_);
  tree->Branch("ele_gen_MomLz2",&ele_gen_MomLz2_);
  tree->Branch("ele_gen_MomLxy2",&ele_gen_MomLxy2_);
  tree->Branch("ele_gen_MomCTau02" ,&ele_gen_MomCTau02_);

  tree->Branch("had_gen_PID" , &had_gen_PID_);
  tree->Branch("had_gen_Status",&had_gen_Status_);
  tree->Branch("had_gen_Charge",&had_gen_Charge_);
  tree->Branch("had_gen_Pt",&had_gen_Pt_);
  tree->Branch("had_gen_Eta",&had_gen_Eta_);
  tree->Branch("had_gen_Phi",&had_gen_Phi_);
  tree->Branch("had_gen_Mass",&had_gen_Mass_);

  tree->Branch("quarks_gen_PID" , &quarks_gen_PID_);
  tree->Branch("quarks_gen_Status",&quarks_gen_Status_);
  tree->Branch("quarks_gen_Charge",&quarks_gen_Charge_);
  tree->Branch("quarks_gen_Pt",&quarks_gen_Pt_);
  tree->Branch("quarks_gen_Eta",&quarks_gen_Eta_);
  tree->Branch("quarks_gen_Mass",&quarks_gen_Mass_);

}

void  BigNtuple::fill_sv_genInfo(const reco::GenParticle prt , const reco::Candidate*  mom){

  float vx = prt.vx(), vy = prt.vy(), vz = prt.vz();
  float mx = mom->vx(), my = mom->vy(), mz = mom->vz();
  float dx = vx - mx, dy = vy - my, dz = vz - mz;

  if(abs(prt.pdgId()) == 13 ){
    mu_gen_PID2_.push_back(prt.pdgId() ); 
    mu_gen_Status2_.push_back(prt.status());
    // genkinematics                                                                                                                          
    mu_gen_Pt2_.push_back(prt.pt());
    mu_gen_Eta2_.push_back(prt.eta());
    mu_gen_Phi2_.push_back(prt.phi());
    mu_gen_Charge2_.push_back(prt.charge());
    // vertexposition                                                                                                                         
    mu_gen_VX2_.push_back(vx);
    mu_gen_VY2_.push_back(vy);
    mu_gen_VZ2_.push_back(vz);
  // vertexflight                                                                                                                           
    mu_gen_Lxy2_.push_back(std::sqrt( vx * vx + vy * vy ));
    mu_gen_Lxyz2_.push_back(std::sqrt( vx * vx + vy * vy + vz * vz));
    // mother beta, gamma,ctau                                                                                                                
    float   beta_mom  = mom->p() / mom->energy();
    float   gamma_mom = mom->energy() / mom->mass();
    // mother quantities related to thedecay                                                                                                  
    // gen id andstatus                                                                                                                       
    mu_gen_MomPID2_.push_back(mom->pdgId());
    mu_gen_MomStatus2_.push_back(mom->status());
    mu_gen_MomCharge2_.push_back(mom->charge());
    mu_gen_MomMass2_.push_back(mom->mass());
    mu_gen_MomPt2_.push_back(mom->pt());
    mu_gen_MomEta2_.push_back(mom->eta());
    mu_gen_MomPhi2_.push_back(mom->phi());
    mu_gen_MomBeta2_.push_back(beta_mom);
    mu_gen_MomGamma2_.push_back(gamma_mom);
    mu_gen_MomLxyz2_.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
    mu_gen_MomLz2_.push_back(dz);
    mu_gen_MomLxy2_.push_back(std::sqrt(dx*dx + dy*dy));
    mu_gen_MomCTau02_.push_back(std::sqrt(dx*dx + dy*dy + dz*dz) / (beta_mom * gamma_mom));
  }

  if(abs(prt.pdgId()) == 11 ){
    ele_gen_PID2_.push_back(prt.pdgId() ); 
    ele_gen_Status2_.push_back(prt.status());
    // genkinematics                                                                                                                          
    ele_gen_Pt2_.push_back(prt.pt());
    ele_gen_Eta2_.push_back(prt.eta());
    ele_gen_Phi2_.push_back(prt.phi());
    ele_gen_Charge2_.push_back(prt.charge());
    // vertexposition                                                                                                                         
    ele_gen_VX2_.push_back(vx);
    ele_gen_VY2_.push_back(vy);
    ele_gen_VZ2_.push_back(vz);
  // vertexflight                                                                                                                           
    ele_gen_Lxy2_.push_back(std::sqrt( vx * vx + vy * vy ));
    ele_gen_Lxyz2_.push_back(std::sqrt( vx * vx + vy * vy + vz * vz));
    // mother beta, gamma,ctau                                                                                                                
    float   beta_mom  = mom->p() / mom->energy();
    float   gamma_mom = mom->energy() / mom->mass();
    // mother quantities related to thedecay                                                                                                  
    // gen id andstatus                                                                                                                       
    ele_gen_MomPID2_.push_back(mom->pdgId());
    ele_gen_MomStatus2_.push_back(mom->status());
    ele_gen_MomCharge2_.push_back(mom->charge());
    ele_gen_MomMass2_.push_back(mom->mass());
    ele_gen_MomPt2_.push_back(mom->pt());
    ele_gen_MomEta2_.push_back(mom->eta());
    ele_gen_MomPhi2_.push_back(mom->phi());
    ele_gen_MomBeta2_.push_back(beta_mom);
    ele_gen_MomGamma2_.push_back(gamma_mom);
    ele_gen_MomLxyz2_.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
    ele_gen_MomLz2_.push_back(dz);
    ele_gen_MomLxy2_.push_back(std::sqrt(dx*dx + dy*dy));
    ele_gen_MomCTau02_.push_back(std::sqrt(dx*dx + dy*dy + dz*dz) / (beta_mom * gamma_mom));
  }
  // final state hadrons
  if( prt.status() == 1 && 
      prt.mother()->pdgId() != 9900014 && 
      prt.mother()->pdgId() != 13 &&
      prt.mother()->pdgId() != 11 ){

    had_gen_PID_.push_back(prt.pdgId() );
    had_gen_Status_.push_back(prt.status()); 
    had_gen_Pt_.push_back(prt.pt());
    had_gen_Eta_.push_back(prt.eta());
    had_gen_Phi_.push_back(prt.phi());
    had_gen_Charge_.push_back(prt.charge());
    had_gen_Mass_.push_back(prt.charge());
  }
  // quarks @ sv
  if( prt.status() == 23 &&
      prt.mother()->pdgId() == 9900014 &&
      (abs(prt.pdgId()) != 13 || abs(prt.pdgId()) != 11)  ){
    quarks_gen_PID_.push_back(prt.pdgId() );
    quarks_gen_Status_.push_back(prt.status());
    quarks_gen_Pt_.push_back(prt.pt());
    quarks_gen_Eta_.push_back(prt.eta());
    quarks_gen_Phi_.push_back(prt.phi());
    quarks_gen_Charge_.push_back(prt.charge());
    quarks_gen_Mass_.push_back(prt.charge());
  }
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

}

void BigNtuple::fill_trigInfo(const edm::TriggerResults& triggerResults, const edm::TriggerNames& trigNames){

  for (size_t i = 0; i < trigNames.size(); ++i) {
    const std::string &name = trigNames.triggerName(i);
    bool fired = triggerResults.accept(i);
    if(!fired) continue;

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

  }

  passIsoMu24All_ = passIsoMu24All_   || passIsoMu24_ || passIsoMuTk24_ ;
  passIsoMu27All_ = passIsoMu27All_   || passIsoMu27_ || passIsoMuTk27_ ;

}

void BigNtuple::set_pileupInfo(TTree* tree){

  tree->Branch("npT" , &npT_, "npT/O");
  tree->Branch("npIT" , &npIT_, "npIT/O");
  tree->Branch("PU_Weight" , &pu_Weight_, "PU_Weight/O");
  tree->Branch("PU_WeightUp" , &pu_WeightUp_, "PU_WeightUp/O");
  tree->Branch("PU_WeightDown" , &pu_WeightDown_, "PU_WeightDown/O");  
}

void BigNtuple::fill_pileupInfo( float npT, float npIT, float PU_Weight, float PU_WeightUp, float PU_WeightDown){

 npT_ = npT;
 npIT_ = npIT;
 pu_Weight_ = PU_Weight;
 pu_WeightUp_ = PU_WeightUp;
 pu_WeightDown_ = PU_WeightDown;

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



void BigNtuple::fill_muInfo(const pat::Muon& mu, const reco::Vertex& pv , double Rho ,double match1 , double match2){

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

  mu_rhoIso_.push_back(Rho); //transverse momentum per unit area                                                                                     

  if(mu.globalTrack().isNonnull() ) {
    mu_normalizedChi2_.push_back(mu.globalTrack()->normalizedChi2());
    mu_numberOfValidPixelHits_.push_back(mu.innerTrack()->hitPattern().numberOfValidPixelHits());
    mu_numberOfValidMuonHits_.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
    mu_numberOftrackerLayersWithMeasurement_.push_back(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
    mu_numberOfMatchedStations_.push_back(mu.numberOfMatchedStations());
    mu_numberOfpixelLayersWithMeasurement_.push_back(mu.innerTrack()->hitPattern().pixelLayersWithMeasurement());
    mu_InnerTrackQuality_.push_back(mu.innerTrack()->quality(reco::TrackBase::highPurity));
  }
  else{
    mu_normalizedChi2_.push_back(-999);
    mu_numberOfValidPixelHits_.push_back(-999);
    mu_numberOfValidMuonHits_.push_back(-999);
    mu_numberOftrackerLayersWithMeasurement_.push_back(-999);
    mu_numberOfMatchedStations_.push_back(-999);
    mu_numberOfpixelLayersWithMeasurement_.push_back(-999);
    mu_InnerTrackQuality_.push_back(-999);
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
  
  reco::MuonTime tofAll = mu.time();
  mu_STATofDirection_.push_back(tofAll.direction());
  mu_STATofNDof_.push_back(tofAll.nDof);
  mu_STATofTimeAtIpInOut_.push_back(tofAll.timeAtIpInOut);
  mu_STATofTimeAtIpInOutErr_.push_back(tofAll.timeAtIpInOutErr);
  mu_STATofTimeAtIpOutIn_.push_back(tofAll.timeAtIpOutIn);
  mu_STATofTimeAtIpOutInErr_.push_back(tofAll.timeAtIpOutInErr);
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
  //============= Parameters related to PF isolation =====================                                                                     
  mu_pfSumPUPt_.push_back(mu.pfIsolationR03().sumPhotonEt);
  mu_PFSumPhotonEt_.push_back(mu.pfIsolationR03().sumPhotonEt);
  mu_pfSumChargedHadronPt_.push_back(mu.pfIsolationR03().sumChargedHadronPt);
  mu_pfSumNeutralHadronEt_.push_back(mu.pfIsolationR03().sumNeutralHadronEt);
  
}

void BigNtuple::set_svInfo(TTree* tree){

    tree->Branch("sv_mu_TrackSize" , &sv_mu_TrackSize_);
    tree->Branch("sv_mu_LXYSig" , &sv_mu_LXYSig_);
    tree->Branch("sv_mu_LXYZSig" , &sv_mu_LXYZSig_);
    tree->Branch("sv_mu_LXY" , &sv_mu_LXY_);
    tree->Branch("sv_mu_LXYZ" , &sv_mu_LXYZ_);
    tree->Branch("sv_mu_mass" , &sv_mu_mass_);
    tree->Branch("sv_mu_charge" , &sv_mu_charge_);
    tree->Branch("sv_mu_eta" , &sv_mu_eta_);
    tree->Branch("sv_mu_phi" , &sv_mu_phi_);
    tree->Branch("sv_mu_pt" , &sv_mu_pt_);
    tree->Branch("sv_mu_p" , &sv_mu_p_);
    tree->Branch("sv_mu_Beta" , &sv_mu_Beta_);
    tree->Branch("sv_mu_Gamma" , &sv_mu_Gamma_);
    tree->Branch("sv_mu_CTau0" , &sv_mu_CTau0_);
    tree->Branch("sv_mu_NDof" , &sv_mu_NDof_);
    tree->Branch("sv_mu_Chi2" , &sv_mu_Chi2_);
    tree->Branch("sv_mu_Angle3D" , &sv_mu_Angle3D_);
    tree->Branch("sv_mu_Angle2D" , &sv_mu_Angle2D_);
    tree->Branch("sv_mu_tracks_charge" , &sv_mu_tracks_charge_);
    tree->Branch("sv_mu_tracks_eta" , &sv_mu_tracks_eta_);
    tree->Branch("sv_mu_tracks_phi" , &sv_mu_tracks_phi_);
    tree->Branch("sv_mu_tracks_pt" , &sv_mu_tracks_pt_);
    tree->Branch("sv_mu_tracks_dxySig" , &sv_mu_tracks_dxySig_);
    tree->Branch("sv_mu_tracks_dxy" , &sv_mu_tracks_dxy_);
    tree->Branch("sv_mu_tracks_dxyz" , &sv_mu_tracks_dxyz_);
    tree->Branch("sv_mu_tracks_Sumcharge" , &sv_mu_tracks_Sumcharge_);
    tree->Branch("sv_mu_tracks_Sumpt" , &sv_mu_tracks_Sumpt_);
    tree->Branch("sv_mu_match" , &sv_mu_match_);

    tree->Branch("sv_ele_TrackSize" , &sv_ele_TrackSize_);
    tree->Branch("sv_ele_LXYSig" , &sv_ele_LXYSig_);
    tree->Branch("sv_ele_LXYZSig" , &sv_ele_LXYZSig_);
    tree->Branch("sv_ele_LXY" , &sv_ele_LXY_);
    tree->Branch("sv_ele_LXYZ" , &sv_ele_LXYZ_);
    tree->Branch("sv_ele_mass" , &sv_ele_mass_);
    tree->Branch("sv_ele_charge" , &sv_ele_charge_);
    tree->Branch("sv_ele_eta" , &sv_ele_eta_);
    tree->Branch("sv_ele_phi" , &sv_ele_phi_);
    tree->Branch("sv_ele_pt" , &sv_ele_pt_);
    tree->Branch("sv_ele_p" , &sv_ele_p_);
    tree->Branch("sv_ele_Beta" , &sv_ele_Beta_);
    tree->Branch("sv_ele_Gamma" , &sv_ele_Gamma_);
    tree->Branch("sv_ele_CTau0" , &sv_ele_CTau0_);
    tree->Branch("sv_ele_NDof" , &sv_ele_NDof_);
    tree->Branch("sv_ele_Chi2" , &sv_ele_Chi2_);
    tree->Branch("sv_ele_Angle3D" , &sv_ele_Angle3D_);
    tree->Branch("sv_ele_Angle2D" , &sv_ele_Angle2D_);
    tree->Branch("sv_ele_tracks_charge" , &sv_ele_tracks_charge_);
    tree->Branch("sv_ele_tracks_eta" , &sv_ele_tracks_eta_);
    tree->Branch("sv_ele_tracks_phi" , &sv_ele_tracks_phi_);
    tree->Branch("sv_ele_tracks_pt" , &sv_ele_tracks_pt_);
    tree->Branch("sv_ele_tracks_dxySig" , &sv_ele_tracks_dxySig_);
    tree->Branch("sv_ele_tracks_dxy" , &sv_ele_tracks_dxy_);
    tree->Branch("sv_ele_tracks_dxyz" , &sv_ele_tracks_dxyz_);
    tree->Branch("sv_ele_tracks_Sumcharge" , &sv_ele_tracks_Sumcharge_);
    tree->Branch("sv_ele_tracks_Sumpt" , &sv_ele_tracks_Sumpt_);
    tree->Branch("sv_ele_match" , &sv_ele_match_);
    tree->Branch("sv_ele_score" , &sv_ele_score_);

}


void BigNtuple::fill_sv_mu_Info(const reco::Vertex& bestVertex, const reco::Vertex& pv , double match){

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

  sv_mu_TrackSize_.push_back(bestVertex.nTracks());
  sv_mu_LXY_.push_back(std::sqrt( dx * dx + dy * dy ));
  sv_mu_LXYZ_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz ));
  sv_mu_LXYSig_.push_back(std::sqrt( dx * dx + dy * dy ) / std::sqrt(xE * xE + yE * yE));
  sv_mu_LXYZSig_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz) / std::sqrt(xE * xE + yE * yE + zE * zE));
  sv_mu_mass_.push_back(bestVertex.p4().mass());
  sv_mu_eta_.push_back(bestVertex.p4().eta());
  sv_mu_phi_.push_back(bestVertex.p4().phi());
  sv_mu_pt_.push_back(bestVertex.p4().pt());
  sv_mu_p_.push_back(bestVertex.p4().P());
  sv_mu_Beta_.push_back(beta_mom);
  sv_mu_Gamma_.push_back(gamma_mom);
  sv_mu_CTau0_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz) / (beta_mom * gamma_mom));
  sv_mu_NDof_.push_back(svNDof);
  sv_mu_Chi2_.push_back(svChi2);
  sv_mu_Angle3D_.push_back(svAngle3D);
  sv_mu_Angle2D_.push_back(svAngle2D);

  int ch = 0;
  float pt = 0;

  reco::Vertex::trackRef_iterator tt = bestVertex.tracks_begin();
  for(; tt != bestVertex.tracks_end(); ++tt) {
    
    sv_mu_tracks_charge_.push_back((*tt)->charge());
    sv_mu_tracks_eta_.push_back((*tt)->eta());
    sv_mu_tracks_phi_.push_back((*tt)->phi());
    sv_mu_tracks_pt_.push_back((*tt)->pt());
    sv_mu_tracks_dxySig_.push_back(fabs((*tt)->dxy(pv.position()))/fabs((*tt)->dxyError()));
    sv_mu_tracks_dxy_.push_back((*tt)->dxy(pv.position()));
    
    ROOT::Math::SVector<double, 3> lxyz1((*tt)->vx()-pv.position().x(), (*tt)->vy()-pv.position().y(), (*tt)->vz()-pv.position().z());
    float dxyz = (float)ROOT::Math::Mag(lxyz1); // magntude of the vector

    sv_mu_tracks_dxyz_.push_back(dxyz);
    ch+=(*tt)->charge();
    pt+=(*tt)->pt();
  }

  sv_mu_tracks_Sumcharge_.push_back(ch);
  sv_mu_tracks_Sumpt_.push_back(pt);
  sv_mu_match_.push_back(match);

}

void BigNtuple::fill_sv_ele_Info(const reco::Vertex& bestVertex, const reco::Vertex& pv , double match, double score){

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

  sv_ele_TrackSize_.push_back(bestVertex.nTracks());
  sv_ele_LXY_.push_back(std::sqrt( dx * dx + dy * dy ));
  sv_ele_LXYZ_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz ));
  sv_ele_LXYSig_.push_back(std::sqrt( dx * dx + dy * dy ) / std::sqrt(xE * xE + yE * yE));
  sv_ele_LXYZSig_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz) / std::sqrt(xE * xE + yE * yE + zE * zE));
  sv_ele_mass_.push_back(bestVertex.p4().mass());
  sv_ele_eta_.push_back(bestVertex.p4().eta());
  sv_ele_phi_.push_back(bestVertex.p4().phi());
  sv_ele_pt_.push_back(bestVertex.p4().pt());
  sv_ele_p_.push_back(bestVertex.p4().P());
  sv_ele_Beta_.push_back(beta_mom);
  sv_ele_Gamma_.push_back(gamma_mom);
  sv_ele_CTau0_.push_back(std::sqrt( dx * dx + dy * dy + dz * dz) / (beta_mom * gamma_mom));
  sv_ele_NDof_.push_back(svNDof);
  sv_ele_Chi2_.push_back(svChi2);
  sv_ele_Angle3D_.push_back(svAngle3D);
  sv_ele_Angle2D_.push_back(svAngle2D);
  sv_ele_score_.push_back(score);

  int ch = 0;
  float pt = 0;

  reco::Vertex::trackRef_iterator tt = bestVertex.tracks_begin();
  for(; tt != bestVertex.tracks_end(); ++tt) {
    
    sv_ele_tracks_charge_.push_back((*tt)->charge());
    sv_ele_tracks_eta_.push_back((*tt)->eta());
    sv_ele_tracks_phi_.push_back((*tt)->phi());
    sv_ele_tracks_pt_.push_back((*tt)->pt());
    sv_ele_tracks_dxySig_.push_back(fabs((*tt)->dxy(pv.position()))/fabs((*tt)->dxyError()));
    sv_ele_tracks_dxy_.push_back((*tt)->dxy(pv.position()));
    
    ROOT::Math::SVector<double, 3> lxyz1((*tt)->vx()-pv.position().x(), (*tt)->vy()-pv.position().y(), (*tt)->vz()-pv.position().z());
    float dxyz = (float)ROOT::Math::Mag(lxyz1); // magntude of the vector

    sv_ele_tracks_dxyz_.push_back(dxyz);
    ch+=(*tt)->charge();
    pt+=(*tt)->pt();
  }

  sv_ele_tracks_Sumcharge_.push_back(ch);
  sv_ele_tracks_Sumpt_.push_back(pt);
  sv_ele_match_.push_back(match);

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
  tree->Branch("jet_btag_pt",&jet_btag_pt_);
  tree->Branch("jet_btag_eta",&jet_btag_eta_);
  tree->Branch("jet_btag_phi",&jet_btag_phi_);
  tree->Branch("jet_btag_flavor",&jet_btag_flavor_);
  tree->Branch("jet_btag_pfCSVv2IVF_discriminator",&jet_btag_pfCSVv2IVF_discriminator_);
}

void BigNtuple::fill_bjetInfo(const pat::Jet& jet,  const std::string & bDiscr, int flavor){

  jet_btag_pt_.push_back(jet.pt());
  jet_btag_eta_.push_back(jet.eta());
  jet_btag_phi_.push_back(jet.phi());
  jet_btag_flavor_.push_back(flavor);
  jet_btag_pfCSVv2IVF_discriminator_.push_back(jet.bDiscriminator(bDiscr));

}

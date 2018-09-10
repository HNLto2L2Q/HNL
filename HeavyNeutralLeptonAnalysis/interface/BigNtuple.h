#ifndef HNL_HeavyNeutralLeptonAnalysis_BigNtuple
#define HNL_HeavyNeutralLeptonAnalysis_BigNtuple
/* 
	 Class: BigNtuple
	 Simple interface class to hide all the ROOT I/O from the plugin and make it more readable
*/

#include <vector>

#include "TTree.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/Handle.h"

class BigNtuple {
public:
	BigNtuple(){} //default, empty constructor

	void set_evtInfo(TTree* tree);
	void fill_evtInfo(const edm::EventID& id);

	void set_pv_genInfo(TTree* tree);
        void fill_pv_genInfo(const reco::GenParticle prt , const reco::Candidate*  mom);

	void set_sv_genInfo(TTree* tree);
        void fill_sv_genInfo(const reco::GenParticle prt , const reco::Candidate*  mom);

	void set_pvInfo(TTree* tree);
	void fill_pvInfo(const reco::VertexCollection& pvs);

	void set_trigInfo(TTree* tree);
	void fill_trigInfo(const edm::TriggerResults& triggerResults, const edm::TriggerNames& trigNames);

	void set_pileupInfo(TTree* tree);
        void fill_pileupInfo( float npT, float npIT, float pu_Weight, float pu_WeightUp, float pu_WeightDown);

	void set_svInfo(TTree* tree);
        void fill_sv_mu_Info(const reco::Vertex& bestVertex, const reco::Vertex& pv, double match);
        void fill_sv_ele_Info(const reco::Vertex& bestVertex, const reco::Vertex& pv, double match , double score);

	void set_muInfo(TTree* tree);
        void fill_muInfo(const pat::Muon& mu, const reco::Vertex& pv, double Rho , double match1 , double match2 );
	
	void set_jetInfo(TTree* tree);
	void fill_jetInfo(const pat::Jet& jet);

        void set_metInfo(TTree* tree);
        void fill_metInfo(const pat::MET& met);

        void set_eleInfo(TTree* tree);
        void fill_eleInfo(const pat::Electron& ele_ , const reco::Vertex& pv, double Rho, double match1, double match2 , std::auto_ptr<EcalClusterLazyTools> recHitEcal);


        void set_eleIDInfo(TTree* tree);
        void fill_eleIDInfo(float ele_mva , bool ele_veto , bool ele_loose , bool ele_medium , bool ele_tight);

        void set_bjetInfo(TTree* tree);
	void fill_bjetInfo(const pat::Jet& jet,  const std::string & bDiscr, int flavor);
 
	void reset() {
	  BigNtuple dummy; //create a new one
	  *this = dummy; //use assignment to reset
	}

private:
	unsigned int lumi_ = 0;
	unsigned int run_ = 0;
	unsigned long long evt_ = 0;

	float npT_ = 0; 
	float npIT_ = 0;
	float pu_Weight_ = 0;
	float pu_WeightUp_ = 0;
	float pu_WeightDown_ = 0;
	
	// primary vertex infos  -- they shouldn't be vector 
	float pvX_ = -1000;
	float pvY_ = -1000;
	float pvZ_ = -1000;
	float pvXErr_ = -1000;
	float pvYErr_ = -1000;
	float pvZErr_ = -1000;
	float pvMass_ = -1000;
	float pvLxy_  = -1000;
	float pvLxyz_ = -1000;
	float pvLxySig_  = -1000;
	float pvLxyzSig_ = -1000;
	float pvChi2_ = -1000;
	int pvNTrack_ = -1000;
	float pvSumPtSq_ = -1000;
	int numberPV_    = -1000;

	//gen infos mu @ pv
	std::vector<int>     mu_gen_PID1_;
	std::vector<int>     mu_gen_Status1_;
	std::vector<int>     mu_gen_Charge1_;
	std::vector<float>   mu_gen_Pt1_;
	std::vector<float>   mu_gen_Eta1_;
	std::vector<float>   mu_gen_Phi1_;
	std::vector<float>   mu_gen_VX1_;
	std::vector<float>   mu_gen_VY1_;
	std::vector<float>   mu_gen_VZ1_;
	std::vector<float>   mu_gen_Lxy1_;
	std::vector<float>   mu_gen_Lxyz1_;
	std::vector<int>     mu_gen_MomPID1_;
	std::vector<int>     mu_gen_MomStatus1_;
	std::vector<float>   mu_gen_MomMass1_;
	std::vector<int>     mu_gen_MomCharge1_;
	std::vector<float>   mu_gen_MomPt1_;
	std::vector<float>   mu_gen_MomEta1_;
	std::vector<float>   mu_gen_MomPhi1_;
	std::vector<float>   mu_gen_MomBeta1_;
	std::vector<float>   mu_gen_MomGamma1_;
	std::vector<float>   mu_gen_MomLxyz1_;
	std::vector<float>   mu_gen_MomLz1_;
	std::vector<float>   mu_gen_MomLxy1_;

	//gen Info mu @ sv
	std::vector<int>     mu_gen_PID2_;
	std::vector<int>     mu_gen_Status2_;
	std::vector<int>     mu_gen_Charge2_;
	std::vector<float>   mu_gen_Pt2_;
	std::vector<float>   mu_gen_Eta2_;
	std::vector<float>   mu_gen_Phi2_;
	std::vector<float>   mu_gen_VX2_;
	std::vector<float>   mu_gen_VY2_;
	std::vector<float>   mu_gen_VZ2_;
	std::vector<float>   mu_gen_Lxy2_;
	std::vector<float>   mu_gen_Lxyz2_;
	std::vector<int>     mu_gen_MomPID2_;
	std::vector<float>   mu_gen_MomMass2_;
	std::vector<int>     mu_gen_MomCharge2_;
	std::vector<int>     mu_gen_MomStatus2_;
	std::vector<float>   mu_gen_MomPt2_;
	std::vector<float>   mu_gen_MomEta2_;
	std::vector<float>   mu_gen_MomPhi2_;
	std::vector<float>   mu_gen_MomBeta2_;
	std::vector<float>   mu_gen_MomGamma2_;
	std::vector<float>   mu_gen_MomLxyz2_;
	std::vector<float>   mu_gen_MomLz2_;
	std::vector<float>   mu_gen_MomLxy2_;
	std::vector<float>   mu_gen_MomCTau02_;

	//gen infos ele @ pv
	std::vector<int>     ele_gen_PID1_;
	std::vector<int>     ele_gen_Status1_;
	std::vector<int>     ele_gen_Charge1_;
	std::vector<float>   ele_gen_Pt1_;
	std::vector<float>   ele_gen_Eta1_;
	std::vector<float>   ele_gen_Phi1_;
	std::vector<float>   ele_gen_VX1_;
	std::vector<float>   ele_gen_VY1_;
	std::vector<float>   ele_gen_VZ1_;
	std::vector<float>   ele_gen_Lxy1_;
	std::vector<float>   ele_gen_Lxyz1_;
	std::vector<int>     ele_gen_MomPID1_;
	std::vector<int>     ele_gen_MomStatus1_;
	std::vector<float>   ele_gen_MomMass1_;
	std::vector<int>     ele_gen_MomCharge1_;
	std::vector<float>   ele_gen_MomPt1_;
	std::vector<float>   ele_gen_MomEta1_;
	std::vector<float>   ele_gen_MomPhi1_;
	std::vector<float>   ele_gen_MomBeta1_;
	std::vector<float>   ele_gen_MomGamma1_;
	std::vector<float>   ele_gen_MomLxyz1_;
	std::vector<float>   ele_gen_MomLz1_;
	std::vector<float>   ele_gen_MomLxy1_;

	//gen Info ele @ sv
	std::vector<int>     ele_gen_PID2_;
	std::vector<int>     ele_gen_Status2_;
	std::vector<int>     ele_gen_Charge2_;
	std::vector<float>   ele_gen_Pt2_;
	std::vector<float>   ele_gen_Eta2_;
	std::vector<float>   ele_gen_Phi2_;
	std::vector<float>   ele_gen_VX2_;
	std::vector<float>   ele_gen_VY2_;
	std::vector<float>   ele_gen_VZ2_;
	std::vector<float>   ele_gen_Lxy2_;
	std::vector<float>   ele_gen_Lxyz2_;
	std::vector<int>     ele_gen_MomPID2_;
	std::vector<float>   ele_gen_MomMass2_;
	std::vector<int>     ele_gen_MomCharge2_;
	std::vector<int>     ele_gen_MomStatus2_;
	std::vector<float>   ele_gen_MomPt2_;
	std::vector<float>   ele_gen_MomEta2_;
	std::vector<float>   ele_gen_MomPhi2_;
	std::vector<float>   ele_gen_MomBeta2_;
	std::vector<float>   ele_gen_MomGamma2_;
	std::vector<float>   ele_gen_MomLxyz2_;
	std::vector<float>   ele_gen_MomLz2_;
	std::vector<float>   ele_gen_MomLxy2_;
	std::vector<float>   ele_gen_MomCTau02_;

	// final state hadrons
	std::vector<int>     had_gen_PID_;
	std::vector<int>     had_gen_Status_;
	std::vector<int>     had_gen_Charge_;
	std::vector<float>   had_gen_Pt_;
	std::vector<float>   had_gen_Eta_;
	std::vector<float>   had_gen_Phi_;
	std::vector<float>   had_gen_Mass_;

	//quarks @ gen
	std::vector<int>     quarks_gen_PID_;
	std::vector<int>     quarks_gen_Status_;
	std::vector<int>     quarks_gen_Charge_;
	std::vector<float>   quarks_gen_Pt_;
	std::vector<float>   quarks_gen_Eta_;
	std::vector<float>   quarks_gen_Phi_;
	std::vector<float>   quarks_gen_Mass_;

	//trigger infos

	bool passIsoMuTk18_  = 0;
	bool passIsoMuTk20_  = 0;
	bool passIsoMuTk22_  = 0;
	bool passIsoMuTk24_  = 0;
	bool passIsoMuTk27_  = 0;
	bool passIsoMuTk17e_ = 0;
	bool passIsoMuTk22e_ = 0;
	
	bool passIsoMu18_  = 0;
	bool passIsoMu20_  = 0;
	bool passIsoMu22_  = 0;
	bool passIsoMu24_  = 0;
	bool passIsoMu27_  = 0;
	bool passIsoMu17e_ = 0;
	bool passIsoMu22e_ = 0;
	bool passTkMu17_   = 0;
	bool passTkMu20_   = 0;
	
	bool passIsoMu24All_ = 0;
	bool passIsoMu27All_ = 0;
	
	bool passDoubleMu17TrkIsoMu8_     = 0;
	bool passDoubleMu17TrkIsoTkMu8_   = 0;
	bool passDoubleTkMu17TrkIsoTkMu8_ = 0;

	bool passIsoEle27_              = 0;
	bool passNonIsoEle115_          = 0;
	bool passDoubleEle23andEle12DZ_ = 0;
	bool passDoubleEle23andEle12_   = 0;

	bool passDoubleEle33TrkMW_      = 0;
	bool passDoubleEle33MW_         = 0;
	bool passDoubleEle33_           = 0;

	bool passDoubleMu33Ele33_       = 0;

 
	//secondary verteces info due to mu
	std::vector<int>   sv_mu_TrackSize_;
	std::vector<float> sv_mu_LXYSig_;
	std::vector<float> sv_mu_LXYZSig_;
	std::vector<float> sv_mu_LXY_;
	std::vector<float> sv_mu_LXYZ_;
	std::vector<float> sv_mu_mass_;
	std::vector<int>   sv_mu_charge_;
	std::vector<float> sv_mu_eta_;
	std::vector<float> sv_mu_phi_;
	std::vector<float> sv_mu_pt_;
	std::vector<float> sv_mu_p_;
	std::vector<float> sv_mu_Beta_;
	std::vector<float> sv_mu_Gamma_;
	std::vector<float> sv_mu_CTau0_;
	std::vector<float> sv_mu_NDof_;
	std::vector<float> sv_mu_Chi2_;
	std::vector<float> sv_mu_Angle3D_;
	std::vector<float> sv_mu_Angle2D_;

	std::vector<int  >  sv_mu_tracks_charge_;
	std::vector<float>  sv_mu_tracks_eta_;
	std::vector<float>  sv_mu_tracks_phi_;
	std::vector<float>  sv_mu_tracks_pt_;
	std::vector<float>  sv_mu_tracks_dxySig_;
	std::vector<float>  sv_mu_tracks_dxy_;
	std::vector<float>  sv_mu_tracks_dxyz_;

	std::vector<int  > sv_mu_tracks_Sumcharge_;
	std::vector<float> sv_mu_tracks_Sumpt_;
	std::vector<float> sv_mu_match_;

	//secondary verteces info due to ele
	std::vector<int>   sv_ele_TrackSize_;
	std::vector<float> sv_ele_LXYSig_;
	std::vector<float> sv_ele_LXYZSig_;
	std::vector<float> sv_ele_LXY_;
	std::vector<float> sv_ele_LXYZ_;
	std::vector<float> sv_ele_mass_;
	std::vector<int>   sv_ele_charge_;
	std::vector<float> sv_ele_eta_;
	std::vector<float> sv_ele_phi_;
	std::vector<float> sv_ele_pt_;
	std::vector<float> sv_ele_p_;
	std::vector<float> sv_ele_Beta_;
	std::vector<float> sv_ele_Gamma_;
	std::vector<float> sv_ele_CTau0_;
	std::vector<float> sv_ele_NDof_;
	std::vector<float> sv_ele_Chi2_;
	std::vector<float> sv_ele_Angle3D_;
	std::vector<float> sv_ele_Angle2D_;

	std::vector<int  >  sv_ele_tracks_charge_;
	std::vector<float>  sv_ele_tracks_eta_;
	std::vector<float>  sv_ele_tracks_phi_;
	std::vector<float>  sv_ele_tracks_pt_;
	std::vector<float>  sv_ele_tracks_dxySig_;
	std::vector<float>  sv_ele_tracks_dxy_;
	std::vector<float>  sv_ele_tracks_dxyz_;

	std::vector<int  > sv_ele_tracks_Sumcharge_;
	std::vector<float> sv_ele_tracks_Sumpt_;
	std::vector<float> sv_ele_match_;
	std::vector<float> sv_ele_score_;

	//muon infos
	std::vector<float> mu_en_ ;
	std::vector<float> mu_pt_ ;
	std::vector<float> mu_eta_ ;
	std::vector<float> mu_phi_ ;
	std::vector<float> mu_et_ ;
	std::vector<float> mu_charge_ ;
	std::vector<int>   mu_FirstGenMatch_ ;
	std::vector<int>   mu_SecondGenMatch_ ;
	std::vector<float> mu_trackiso_ ;
	std::vector<float> mu_rhoIso_;
	std::vector<float> mu_pfSumChargedHadronPt_ ;
	std::vector<float> mu_pfSumNeutralHadronEt_ ;
	std::vector<float> mu_PFSumPhotonEt_ ;
	std::vector<float> mu_pfSumPUPt_ ;
	std::vector<int>   mu_numberOfValidMuonHits_ ;
	std::vector<float> mu_emIso_ ;
	std::vector<float> mu_hadIso_ ;
	std::vector<float> mu_normalizedChi2_ ;
	std::vector<int>   mu_numberOfMatchedStations_ ;
	std::vector<int>   mu_numberOfValidPixelHits_ ;
	std::vector<int>   mu_numberOftrackerLayersWithMeasurement_ ;
	std::vector<int>   mu_numberOfpixelLayersWithMeasurement_ ;
	std::vector<int>   mu_TrackQuality_ ;
	std::vector<int>   mu_InnerTrackQuality_ ;
	std::vector<float> mu_pxTunePMuonBestTrack_ ;
	std::vector<float> mu_pyTunePMuonBestTrack_ ;
	std::vector<float> mu_pzTunePMuonBestTrack_ ;
	std::vector<float> mu_pTunePMuonBestTrack_ ;
	std::vector<float> mu_etaTunePMuonBestTrack_ ;
	std::vector<float> mu_LXYZ_ ;
	std::vector<float> mu_LXY_ ;
	std::vector<float> mu_ptTunePMuonBestTrack_ ;
	std::vector<float> mu_phiTunePMuonBestTrack_ ;
	std::vector<float> mu_thetaTunePMuonBestTrack_ ;
	std::vector<float> mu_chargeTunePMuonBestTrack_ ;
	std::vector<float> mu_dPToverPTTunePMuonBestTrack_ ;
	std::vector<float> mu_absdxyTunePMuonBestTrack_ ;
	std::vector<float> mu_absdxyErrorTunePMuonBestTrack_ ;
	std::vector<float> mu_absdxySigTunePMuonBestTrack_ ;
	std::vector<float> mu_absdzTunePMuonBestTrack_ ;
	std::vector<float> mu_absdzErrorTunePMuonBestTrack_ ;
	std::vector<float> mu_absdzSigTunePMuonBestTrack_ ;
	std::vector<float> mu_recoDeltaBeta_ ;
	std::vector<float> mu_recoiso_ ;
	std::vector<float> mu_isGlobalMuon_ ;
	std::vector<float> mu_isStandAloneMuon_ ;
	std::vector<float> mu_isPF_ ;
	std::vector<float> mu_isRPCMuon_ ;
	std::vector<float> mu_isTrackerMuon_ ;
	std::vector<float> mu_isGoodMuon_ ;
	std::vector<float> mu_isSoftMuon_ ;
	std::vector<float> mu_isLoose_ ;
	std::vector<float> mu_isTightMuon_ ;
	std::vector<int>    mu_STAnHits_ ;
	std::vector<int>    mu_STAnLost_ ;
	std::vector<int>    mu_STAnStationsWithAnyHits_ ;
	std::vector<int>    mu_STAnCscChambersWithAnyHits_ ;
	std::vector<int>    mu_STAnDtChambersWithAnyHits_ ;
	std::vector<int>    mu_STAnRpcChambersWithAnyHits_ ;
	std::vector<int>    mu_STAinnermostStationWithAnyHits_ ;
	std::vector<int>    mu_STAoutermostStationWithAnyHits_ ;
	std::vector<int>    mu_STAnStationsWithValidHits_ ;
	std::vector<int>    mu_STAnCscChambersWithValidHits_ ;
	std::vector<int>    mu_STAnDtChambersWithValidHit_ ;
	std::vector<int>    mu_STAnRpcChambersWithValidHits_ ;
	std::vector<int>    mu_STAnValidMuonHits_ ;
	std::vector<int>    mu_STAnValidCscHits_ ;
	std::vector<int>    mu_STAnValidDtHits_ ;
	std::vector<int>    mu_STAnValidRpcHits_ ;
	std::vector<int>    mu_STAinnermostStationWithValidHits_ ;
	std::vector<int>    mu_STAoutermostStationWithValidHits_ ;
	std::vector<float>  mu_STATofDirection_ ;
	std::vector<float>  mu_STATofNDof_ ;
	std::vector<float>  mu_STATofTimeAtIpInOut_ ;
	std::vector<float>  mu_STATofTimeAtIpInOutErr_ ;
	std::vector<float>  mu_STATofTimeAtIpOutIn_ ;
	std::vector<float>  mu_STATofTimeAtIpOutInErr_ ;

	//jet info

	std::vector<float>   jet_charge_ ;
	std::vector<float>   jet_et_ ;
	std::vector<float>   jet_pt_ ;
	std::vector<float>   jet_eta_ ;
	std::vector<float>   jet_phi_ ;
	std::vector<float>   jet_theta_ ;
	std::vector<float>   jet_en_ ;
	std::vector<float>   jet_chargedEmEnergy_ ;
	std::vector<float>   jet_neutralEmEnergyFraction_ ;
	std::vector<float>   jet_chargedHadronEnergy_ ;
	std::vector<float>   jet_neutralHadronEnergyFraction_ ;
	std::vector<float>   jet_chargedMuEnergy_ ;
	std::vector<float>   jet_chargedMuEnergyFraction_ ;
	std::vector<float>   jet_chargedMultiplicity_ ;
	std::vector<float>   jet_numberOfDaughters_ ;
	std::vector<float>   jet_muonEnergy_ ;
	std::vector<float>   jet_muonEnergyFraction_ ;
	std::vector<float>   jet_muonMultiplicity_ ;
	std::vector<float>   jet_neutralEmEnergy_ ;
	std::vector<float>   jet_neutralHadronEnergy_ ;
	std::vector<float>   jet_neutralHadronMultiplicity_ ;
	std::vector<float>   jet_neutralMultiplicity_ ;

	//electron info

	std::vector<float>   ele_Et_;
	std::vector<float>   ele_EtFromCaloEn_;    
	std::vector<float>   ele_pt_; 
	std::vector<float>   ele_etaSC_;
	std::vector<float>   ele_phiSC_;
	std::vector<float>   ele_phiWidth_; 
	std::vector<float>   ele_etaWidth_; 
	std::vector<float>   ele_energySC_;
	std::vector<float>   ele_thetaSC_;
	std::vector<float>   ele_preshowerEnergySC_;
	std::vector<float>   ele_etaTrack_; 
	std::vector<float>   ele_phiTrack_;
	std::vector<float>   ele_thetaTrack_;   
	std::vector<float>   ele_x_;
	std::vector<float>   ele_y_;
	std::vector<float>   ele_z_;  
	std::vector<float>   ele_e2x5Max_;
	std::vector<float>   ele_e1x5_;
	std::vector<float>   ele_e5x5_;
	std::vector<float>   ele_e2x5MaxOver5x5_;
	std::vector<float>   ele_e1x5Over5x5_;
	std::vector<float>   ele_sigmaIetaIetaFull5x5_;
	std::vector<float>   ele_e2x5MaxFull5x5_;
	std::vector<float>   ele_e1x5Full5x5_;
	std::vector<float>   ele_e5x5Full5x5_;
	std::vector<float>   ele_e2x5MaxOver5x5Full5x5_;
	std::vector<float>   ele_e1x5Over5x5Full5x5_;  
	std::vector<float>   ele_zTrackPositionAtVtx_;
	std::vector<float>   ele_hadronicOverEm_;
	std::vector<float>   ele_deltaEtaInSC_;
	std::vector<float>   ele_deltaPhiInSC_;
	std::vector<float>   ele_deltaEtaInSeedCluster_;
	std::vector<float>   ele_deltaPhiInSeedCluster_;
	std::vector<float>   ele_sigmaIetaIeta_;  
	std::vector<float>   ele_e2x5Right_;
	std::vector<float>   ele_e2x5Left_;
	std::vector<float>   ele_e2x5Top_;
	std::vector<float>   ele_e2x5Bottom_;
	std::vector<float>   ele_eMax_;
	std::vector<float>   ele_eRight_;
	std::vector<float>   ele_eLeft_;
	std::vector<float>   ele_eTop_;
	std::vector<float>   ele_eBottom_;
	std::vector<float>   ele_e3x3_;
	std::vector<float>   ele_frac51_;
	std::vector<float>   ele_frac15_;
	
	std::vector<int>   ele_rawId_;
	std::vector<int>   ele_ieta_;
	std::vector<int>   ele_nbOfMissingHits_;
	std::vector<int>   ele_charge_;
	std::vector<bool>  ele_isEcalDrivenSeed_;
	std::vector<bool>  ele_isPassConversionVeto_;
  
	std::vector<float>   ele_dxy_;
	std::vector<float>   ele_dz_; 
	std::vector<float>   ele_rhoIso_;
	std::vector<float>   ele_fbrem_;
	std::vector<float>   ele_EoverP_;
	std::vector<float>   ele_Xposition_;   
	std::vector<float>   ele_Yposition_;   
	std::vector<float>   ele_dr03TkSumPt_;
	std::vector<float>   ele_hcalDepth1OverEcal_;
	std::vector<float>   ele_hcalDepth2OverEcal_;
	std::vector<float>   ele_dr03HcalDepth2TowerSumEt_;
	std::vector<float>   ele_hcalDepth2TowerSumEtNoVeto_; 
	std::vector<float>   ele_hcalDepth1TowerSumEtNoVeto_; 
	std::vector<float>   ele_EcalPlusHcald1iso_;
	std::vector<float>   ele_dr03EcalRecHitSumEt_;
	std::vector<float>   ele_dr03HcalDepth1TowerSumEt_;
	std::vector<float>   ele_dr03HcalDepth1TowerSumEtBc_;
	std::vector<float>   ele_pfSumPhotonEt_;
	std::vector<float>   ele_pfSumChargedHadronPt_; 
	std::vector<float>   ele_pfSumNeutralHadronEt_;
	std::vector<float>   ele_pfSumPUPt_;  
	std::vector<float>   ele_pfDeltaBeta_;
	std::vector<float>   ele_FirstGenMatch_;
	std::vector<float>   ele_SecondGenMatch_;

	std::vector<float>   ele_Mva2016_;
	std::vector<float>   ele_CutVeto_; 
	std::vector<float>   ele_CutLoose_;
	std::vector<float>   ele_CutMedium_;
	std::vector<float>   ele_CutTight_;

	/*
	std::vector<float>   ele_Mva_;
	std::vector<float>   ele_MvaFall17Iso_;
	std::vector<float>   ele_MvaFall17NoIso_;
	std::vector<float>   ele_CutBasedVeto_;
	std::vector<float>   ele_CutBasedLoose_;
	std::vector<float>   ele_CutBasedMedium_;
	std::vector<float>   ele_CutBasedTight_;
	*/

	//MET info
	
	float  pfMet_et_ = -1000;
	float  pfMet_pt_ = -1000;
	float  pfMet_phi_ = -1000;
	float  pfMet_en_ = -1000;
	float  pfMet_px_ = -1000;
	float  pfMet_py_ = -1000;
	float  pfMet_pz_ = -1000;
	float  pfMet_sumEt_ = -1000;
	float  caloMet_pt_ = -1000;
	float  caloMet_phi_ = -1000;

	//bJet info
	std::vector<int>   jet_btag_flavor_;
	std::vector<float> jet_btag_pfCSVv2IVF_discriminator_;
	std::vector<float> jet_btag_pt_;
	std::vector<float> jet_btag_eta_;
	std::vector<float> jet_btag_phi_;

}; 


#endif

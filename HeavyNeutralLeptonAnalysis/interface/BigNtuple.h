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

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


class BigNtuple {
public:
	BigNtuple(){} //default, empty constructor

	//getter

	//vertex info
        float get_pv_x()const{return lep1_gen_vx_;}
        float get_pv_y()const{return lep1_gen_vy_;}
	float get_pv_z()const{return lep1_gen_vz_;}

        float get_sv_x()const{return lep2_gen_vx_;}
        float get_sv_y()const{return lep2_gen_vy_;}
        float get_sv_z()const{return lep2_gen_vz_;}

	float get_sv_pt()const{return best_sv_pt_;}
	float get_sv_px()const{return best_sv_px_;}
	float get_sv_py()const{return best_sv_py_;}
	float get_sv_pz()const{return best_sv_pz_;}
	float get_sv_en()const{return best_sv_energy_;}

        float get_sv_recox()const{return best_sv_recox_;}
        float get_sv_recoy()const{return best_sv_recoy_;}
        float get_sv_recoz()const{return best_sv_recoz_;}

        float get_met_pt()const{return pfMet_pt_;}
        float get_met_px()const{return pfMet_px_;}
        float get_met_py()const{return pfMet_py_;}

	float get_pvTosv_rho()const{return pvTosv_rho_;}
	float get_pvTosv_phi()const{return pvTosv_phi_;}
	float get_pvTosv_theta()const{return pvTosv_theta_;}

	//phi info
	float get_lep1_phi()const{return lep1_gen_Phi_;}
	float get_lep1_pt()const{return lep1_gen_Pt_;}
	float get_lep1_eta()const{return lep1_gen_Eta_;}
	float get_lep1_en()const{return lep1_gen_En_;}

	float get_lep2_phi()const{return lep2_gen_Phi_;}
	float get_lep2_eta()const{return lep2_gen_Eta_;}


	//trigger info
	bool get_passIsoMu24All()const{return passIsoMu24All_;}
        bool get_passIsoMu27All()const{return passIsoMu27All_;}

	//setter

	//void set_sv_x(float sv_x){lep2_gen_vx_ = sv_x;}
	//void set_sv_y(float sv_y){lep2_gen_vy_ = sv_y;}
	//void set_sv_z(float sv_z){lep2_gen_vz_ = sv_z;}

	void set_evtInfo(TTree* tree);
	void fill_evtInfo(const edm::EventID& id, int& pile_up_info);

	void set_weightsInfo(TTree* tree);
	void fill_weightsInfo(const edm::Handle<GenEventInfoProduct> genEventInfoHandle, edm::Handle<LHEEventProduct> lheinfo);

        void set_prefiring(TTree* tree);
        void fill_prefiring(double weight, double weightup, double weightdown);

	void set_pv_genInfo(TTree* tree);
        void fill_pv_genInfo(const reco::GenParticle prt,const std::vector<reco::GenParticle>);

	void set_sv_genInfo(TTree* tree);
        void fill_sv_genInfo(const reco::GenParticle prt,const std::vector<reco::GenParticle>);

	void set_pvInfo(TTree* tree);
	void fill_pvInfo(const reco::VertexCollection& pvs);

	void set_trigInfo(TTree* tree);
	void fill_trigInfo(const edm::TriggerResults& triggerResults, const edm::TriggerNames& trigNames);


	void set_sv_Info(TTree* tree);
        void fill_sv_Info(const reco::Vertex& bestVertex, const reco::Vertex& pv, std::pair<float, float> match, bool lept);

	void set_muInfo(TTree* tree);
        void fill_muInfo(const pat::Muon& mu, const reco::Vertex& pv, double Rho , std::pair<double, double> match1 , std::pair<double,double> match2 );

	void set_jetInfo(TTree* tree);
	//void fill_jetInfo(const pat::Jet& jet);
	void fill_jetInfo(const pat::Jet& jet, float smeared,float smearedUp ,float smearedDown ,double un , double unSmeared );

        void set_metInfo(TTree* tree);
        void fill_metInfo(const pat::MET& met);

	void set_transverseMassInfo(TTree* tree);
	void fill_transverseMassInfo(float tr_mass_ivf, float tr_mass_lep1, float tr_mass_ivfPluslep1);

	void set_massCorrection(TTree* tree);
	void fill_massCorrection(double mass_corr);

        void set_eleInfo(TTree* tree);
        void fill_eleInfo(const pat::Electron& ele_ , const reco::Vertex& pv, double Rho, std::pair<double,double> match1, std::pair<double,double> match2 , std::auto_ptr<EcalClusterLazyTools> recHitEcal);


        void set_eleIDInfo(TTree* tree);
        void fill_eleIDInfo(float ele_mva , bool ele_veto , bool ele_loose , bool ele_medium , bool ele_tight);

        //void set_bjetInfo(TTree* tree);
	//void fill_bjetInfo(const pat::Jet& jet, int flavor);

	void reset() {
	  BigNtuple dummy; //create a new one
	  *this = dummy; //use assignment to reset
	}

	//float lep2_gen_vx_ = -1000;
	//float lep2_gen_vy_ = -1000;
	//float lep2_gen_vz_ = -1000;
private:

	unsigned int lumi_ = 0;
	unsigned int run_ = 0;
	int tnpv = -1;
	unsigned long long evt_ = 0;

	float gen_weight_ = 0;
	float lhe_weight_ = 0;
	float lhe_ctau_   = 0;

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
	float pvLxySigma_  = -1000;
	float pvLxyzSigma_ = -1000;
	float pvChi2_ = -1000;
	int pvNTrack_ = -1000;
	float pvSumPtSq_ = -1000;
	int numberPV_    = -1000;

	//gen infos mu @ pv
	int     lep1_gen_PID_     = -1000;
	int     lep1_gen_Charge_  = -1000;
	float   lep1_gen_Pt_      = -1000;
	float   lep1_gen_Eta_      = -1000;
	float   lep1_gen_Phi_      = -1000;
	float   lep1_gen_En_       = -1000;
	float   lep1_gen_vx_       = -1000;
	float   lep1_gen_vy_       = -1000;
	float   lep1_gen_vz_       = -1000;
	float   lep1_gen_Lxy_      = -1000;
	float   lep1_gen_Lxyz_     = -1000;
	int     HNL_gen_PID_       = -1000;
	float   HNL_gen_Mass_      = -1000;
	int     HNL_gen_Charge_    = -1000;
	float   HNL_gen_Pt_        = -1000;
	float   HNL_gen_Eta_       = -1000;
	float   HNL_gen_Phi_       = -1000;

	//gen Info mu @ sv
	int     lep2_gen_PID_      = -1000;
	int     lep2_gen_Charge_   = -1000;
	float   lep2_gen_Pt_       = -1000;
	float   lep2_gen_Eta_      = -1000;
	float   lep2_gen_Phi_      = -1000;
	float   lep2_gen_vx_       = -1000;
	float   lep2_gen_vy_       = -1000;
	float   lep2_gen_vz_       = -1000;
	float   lep2_gen_MomLxyz_  = -1000;
	float   lep2_gen_MomLz_    = -1000;
	float   lep2_gen_MomLxy_   = -1000;
	float   lep2_gen_MomCTau0_ = -1000;


        //float get_sv_x(){return lep2_gen_vx_;}
        //float get_sv_y(){return lep2_gen_vy_;}
        //float get_sv_z(){return lep2_gen_vz_;}

	// final state hadrons
	std::vector<int>     daugh_gen_PID_;
	std::vector<int>     daugh_gen_Charge_;
	std::vector<float>   daugh_gen_Pt_;
	std::vector<float>   daugh_gen_Eta_;
	std::vector<float>   daugh_gen_Phi_;
	std::vector<float>   daugh_gen_Mass_;

	//trigger infos
	bool passEle32_WPTight_Gsf_ = 0;
	bool passMu3_PFJet40_    = 0;
	bool passMu8_TrkIsoVVL_  = 0;
	bool passMu17_TrkIsoVVL_ = 0;

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

	float best_sv_px_ = 999;
	float best_sv_py_ = 999;
	float best_sv_pz_ = 999;
	float best_sv_pt_ = 999;
	float best_sv_energy_ = 999;

        float best_sv_recox_ = 999;
        float best_sv_recoy_ = 999;
	float best_sv_recoz_ = 999;

	//secondary verteces info due to mu
	std::vector<bool> sv_hasMuon_;
	std::vector<int>   sv_numTracks_;
	std::vector<float> sv_x_;
	std::vector<float> sv_y_;
	std::vector<float> sv_z_;
	std::vector<float> sv_xErr_;
	std::vector<float> sv_yErr_;
	std::vector<float> sv_zErr_;
	std::vector<float> sv_LxySig_;
	std::vector<float> sv_LxyzSig_;
	std::vector<float> sv_Lxy_;
	std::vector<float> sv_Lxyz_;
	std::vector<float> sv_mass_;
	std::vector<int>   sv_charge_;
	std::vector<float> sv_eta_;
	std::vector<float> sv_phi_;
	std::vector<float> sv_pt_;
	std::vector<float> sv_p_;
	std::vector<float> sv_px_;
	std::vector<float> sv_py_;
	std::vector<float> sv_pz_;
	std::vector<float> sv_energy_;
	std::vector<float> sv_Beta_;
	std::vector<float> sv_Gamma_;
	std::vector<float> sv_CTau0_;
	std::vector<float> sv_NDof_;
	std::vector<float> sv_Chi2_;
	std::vector<float> sv_Angle3D_;
	std::vector<float> sv_Angle2D_;

	std::vector<std::vector<int  > > sv_tracks_charge_;
	std::vector<std::vector<float> > sv_tracks_eta_;
	std::vector<std::vector<float> > sv_tracks_phi_;
	std::vector<std::vector<float> > sv_tracks_pt_;
	std::vector<std::vector<float> > sv_tracks_p_;
	std::vector<std::vector<float> > sv_tracks_dxySig_;
	std::vector<std::vector<float> > sv_tracks_dxy_;
	std::vector<std::vector<float> > sv_tracks_dxyz_;

	std::vector<float> sv_lx_;
	std::vector<float> sv_ly_;
	std::vector<float> sv_lz_;

	std::vector<int  > sv_tracks_Sumcharge_;
	std::vector<float> sv_tracks_Sumpt_;
	std::vector<float> sv_match_dxy_;
	std::vector<float> sv_match_dxyz_;


	//muon infos
	std::vector<float> mu_en_ ;
	std::vector<float> mu_pt_ ;
	std::vector<float> mu_eta_ ;
	std::vector<float> mu_phi_ ;
	std::vector<float> mu_et_ ;
	std::vector<float> mu_charge_ ;
	std::vector<std::pair<double,double>>   mu_FirstGenMatch_ ;
	std::vector<std::pair<double,double>>   mu_SecondGenMatch_ ;
	std::vector<float> mu_trackiso_ ;
	std::vector<float> mu_rhoIso_;
	std::vector<float> mu_pfSumChargedHadronPt_ ;
	std::vector<float> mu_pfSumNeutralHadronEt_ ;
	std::vector<float> mu_PFSumPhotonEt_ ;
	std::vector<float> mu_pfSumPUPt_ ;
	std::vector<int>   mu_numberOfValidMuonHits_ ;
	std::vector<float> mu_emIso_ ;
	std::vector<float> mu_hadIso_ ;
	std::vector<float> mu_segmentCompatibilityMuonBestTrack_;
	std::vector<float> mu_trkKinkMuonBestTrack_;
	std::vector<float> mu_chi2LocalPositionMuonBestTrack_;
	std::vector<float> mu_normalizedChi2_ ;
	std::vector<int>   mu_numberOfMatchedStations_ ;
	std::vector<int>   mu_numberOfValidPixelHits_ ;
	std::vector<int>   mu_numberOftrackerLayersWithMeasurement_ ;
	std::vector<int>   mu_numberOfpixelLayersWithMeasurement_ ;
	std::vector<int>   mu_TrackQuality_ ;
	std::vector<int>   mu_InnerTrackQuality_ ;
	std::vector<float> mu_InnerTrackValidFraction_;
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
	std::vector<float> mu_isPFMuon_ ;
	std::vector<float> mu_isRPCMuon_ ;
	std::vector<float> mu_isTrackerMuon_ ;
	std::vector<float> mu_isGoodMuon_ ;
	std::vector<float> mu_isSoftMuon_ ;
	std::vector<float> mu_isLooseMuon_ ;
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

	int mu_cmb_nDof_ = -1;
	int mu_rpc_nDof_ = -1;
	double mu_cmb_time_ = -999.;
	double mu_rpc_time_ = -999.;
	double mu_cmb_timeErr_ = -999.;
	double mu_rpc_timeErr_ = -999.;

	//jet info
	std::vector<float>   jetPt_JECUp_;
	std::vector<float>   jetPt_JECDown_;
	std::vector<float>   jetSmearedPt_;
	std::vector<float>   jetSmearedPt_JERUp_;
	std::vector<float>   jetSmearedPt_JERDown_;
	std::vector<float>   jetSmearedPt_unUp_;
	std::vector<float>   jetSmearedPt_unDown_;

	std::vector<float>   jet_charge_ ;
	std::vector<float>   jet_et_ ;
	std::vector<float>   jet_pt_ ;
	std::vector<float>   jet_eta_ ;
	std::vector<float>   jet_phi_ ;
	std::vector<float>   jet_theta_ ;
	std::vector<float>   jet_en_ ;
	std::vector<float>   jet_chargedEmEnergy_ ;
	std::vector<float>   jet_chargedEmEnergyFraction_ ;
	std::vector<float>   jet_neutralEmEnergyFraction_ ;
	std::vector<float>   jet_chargedHadronEnergy_ ;
	std::vector<float>   jet_chargedHadronEnergyFraction_ ;
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
	std::vector<float>   jet_pileUpid_ ;
	std::vector<float>   jet_L1ptcorrection_ ;
	std::vector<float>   jet_L2ptcorrection_ ;
	std::vector<float>   jet_L3ptcorrection_ ;
	std::vector<float>   jet_ptuncorrected_ ;
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
	std::vector<std::pair<double,double>>   ele_FirstGenMatch_;
	std::vector<std::pair<double,double>>   ele_SecondGenMatch_;

	std::vector<float>   ele_Mva2016_;
	std::vector<float>   ele_CutVeto_;
	std::vector<float>   ele_CutLoose_;
	std::vector<float>   ele_CutMedium_;
	std::vector<float>   ele_CutTight_;
	std::vector<float>   ele_isEB_;
	std::vector<float>   ele_isEE_;
	std::vector<float>   ele_eSuperClusterOverP_;
	std::vector<float>   ele_ecalEnergy_;
	std::vector<float>   ele_dEtaInSeed_;
	std::vector<float>   ele_InvMinusPInv_;

	std::vector<float>   ele_PtCorr_;
	std::vector<float>   ele_PtScaleUp_;
	std::vector<float>   ele_PtScaleDown_;
	std::vector<float>   ele_PtResUp_;
	std::vector<float>   ele_PtResDown_;
	std::vector<float>   ele_ECorr_;
	std::vector<float>   ele_EScaleUp_;
	std::vector<float>   ele_EScaleDown_;
	std::vector<float>   ele_EResUp_;
	std::vector<float>   ele_EResDown_;

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
	float metJECDown_ = -1000;
	float metJECUp_ = -1000;
	float metUnclDown_ = -1000;
	float metUnclUp_ = -1000;
	float metPhiJECDown_ = -1000;
	float metPhiJECUp_ = -1000;
	float metPhiUnclDown_ = -1000;
	float metPhiUnclUp_ = -1000;

	//transverse mass info
	float tranvsverseMass_ivf_ = -999;
	float tranvsverseMass_lep1_ = -999;
	float tranvsverseMass_ivfPluslep1_ = -999;

	//correction mass info
	float pvTosv_rho_ = -999; //respect to primary vertex
	float pvTosv_phi_ = -999;
	float pvTosv_theta_ = -999;
	float sv_mass_corr_ = -999;

	//bJet info
	std::vector<int>   jet_btag_flavor_;
	std::vector<float> jet_btag_pfDeepCSV_bb_discriminator_;
	std::vector<float> jet_btag_pfDeepCSV_bbb_discriminator_;
	std::vector<float> jet_btag_pfDeepCSV_bc_discriminator_;
	std::vector<float> jet_btag_pt_;
	std::vector<float> jet_btag_eta_;
	std::vector<float> jet_btag_phi_;

	double prefiring_weight_ = -99999;
	double prefiring_weightup_ = -99999;
	double prefiring_weightdown_ = -99999;



};


#endif

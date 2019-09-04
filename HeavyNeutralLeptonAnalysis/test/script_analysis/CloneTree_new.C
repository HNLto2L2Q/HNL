//**********************************************************************************************************************************
// Remove some branches + selects the events + add variables -- for muonic channel
//***************************************** To Compile******************************************************************************
// g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit CloneTree.C -o CloneTree.exe
//**********************************************************************************************************************************

#ifndef __CINT__
//#include "RooGlobalFunc.h"
//------------------------------------------------
     
#endif
#include <iostream>
#include <fstream>
#include "TTree.h"
#include <iostream>


//#include "RooMCStudy.h"
//#include "RooFitResult.h"
//#include "RooStats/SPlot.h"
#include <vector>
#include <string>
#include <iostream>
//#include "RooRandom.h"
//#include "RooMinuit.h"
#include "TRandom3.h"
#include <time.h>
#include <TROOT.h>
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <iostream>
#include <TMath.h>
#include "TH1D.h"
#include "TH2.h"
#include "TString.h"
#include <TLatex.h>
#include "TText.h"
#include "TPaveText.h"
#include "TLorentzVector.h"


using namespace std;



int main(int argc, char **argv){
 
  std::ifstream infile_data;
  infile_data.open(static_cast<std::string>(argv[1]));
  std::string inputfilename_data;

 
  Float_t isoCut = 0.15;
  bool isMC = true;

  TChain * oldtree = new TChain("HeavyNeutralLepton/tree_","");

  std::string name_tree;
  std::string str;
  while (std::getline(infile_data, str)){
  cout << str << "\n"; 

  oldtree->Add(static_cast<TString>(str));
  int cont = 0;
  std::vector<std::string> vs;
  char* token;
  const char* delim = "/";
  token=strtok((char*)str.c_str(), delim);

  while (token != NULL)
    {
  	vs.push_back(std::string(token));
  	token = strtok(NULL, delim);
  	cont++;
    }


    name_tree="skimmed_" + vs[cont-1];
    
    TFile *newfile = new TFile(static_cast<TString>(name_tree),"recreate");
    TTree *newtree  = new TTree("newtree","Analysis Tree");
    Long64_t nentries = oldtree->GetEntriesFast();
    cout  <<nentries<<endl;
    cout << str <<endl;    
    //======================= Old Tree Variables ==========================================// 
    // Trigger variables 
    Bool_t passIsoMu24=0;
    Bool_t  passEle32_WPTight_Gsf=0;
    
    oldtree->SetBranchAddress("passIsoMu24",&passIsoMu24);
    oldtree->SetBranchAddress("passEle32_WPTight_Gsf",&passEle32_WPTight_Gsf);
    
    //vector<Float_t>   *PU_Weight = 0;
    //oldtree->SetBranchAddress("PU_Weight",&PU_Weight);
    
    vector<Float_t>   *mu_isTightMuon_ = 0;
    vector<Float_t>   *mu_isLooseMuon_ = 0;
    vector<Float_t>   *mu_pt_ = 0;
    vector<Float_t>   *mu_eta_ =0;
    vector<Float_t>   *mu_phi_ =0; 
    vector<Float_t>   *mu_charge_ =0; 
    vector<Float_t>   *mu_en_ =0;
    vector<Float_t>   *mu_et_ =0;
    vector<Float_t>   *deltaBeta_ =0;
    vector<Float_t>   *mu_rhoIso_ = 0 ;
    vector<Float_t>   *mu_trackiso_ = 0 ;
    vector<Float_t>   *mu_pfSumChargedHadronPt_ = 0 ;
    vector<Float_t>   *mu_pfSumNeutralHadronEt_ = 0 ;
    vector<Float_t>   *mu_PFSumPhotonEt_ = 0 ;
    vector<Float_t>   *mu_pfSumPUPt_ = 0 ;
    vector<Float_t>   *mu_emIso_ = 0 ;
    vector<Float_t>   *mu_hadIso_ = 0 ;
    //vector<Float_t>   *mu_normalizedChi2_ = 0 ;
    vector<Float_t>   *mu_dPToverPTTunePMuonBestTrack_ = 0 ;
    vector<Float_t>   *mu_absdxyTunePMuonBestTrack_ = 0 ;
    vector<Float_t>   *mu_absdxyErrorTunePMuonBestTrack_ = 0 ;
    vector<Float_t>   *mu_absdxySigTunePMuonBestTrack_ = 0 ;
    vector<Float_t>   *mu_absdzTunePMuonBestTrack_ = 0 ;
    vector<Float_t>   *mu_absdzErrorTunePMuonBestTrack_ = 0 ;
    vector<Float_t>   *mu_absdzSigTunePMuonBestTrack_ = 0 ;
    vector<Float_t>   *mu_recoiso_ = 0 ;
    // vector<Float_t>   *mu_STATofDirection_ = 0 ;
    vector<Float_t>   *mu_STATofNDof_ = 0 ;
    vector<Float_t>   *mu_STATofTimeAtIpInOut_ = 0 ;
    vector<Float_t>   *mu_STATofTimeAtIpInOutErr_ = 0 ;
    vector<Float_t>   *mu_STATofTimeAtIpOutIn_ = 0 ;
    vector<Float_t>   *mu_STATofTimeAtIpOutInErr_ = 0 ;
    //vector<Int_t>     *mu_numberOfValidMuonHits_ = 0 ;
    //vector<Int_t>     *mu_numberOfMatchedStations_ = 0 ;
    //vector<Int_t>     *mu_numberOfValidPixelHits_ = 0 ;
    //vector<Int_t>     *mu_TrackQuality_ = 0 ;
    // vector<Int_t>     *mu_InnerTrackQuality_ = 0 ;
    vector<Int_t>     *mu_FirstGenMatch_ = 0;
    vector<Int_t>     *mu_SecondGenMatch_ = 0;
    
    vector<Int_t>   *sv_mu_TrackSize_ = 0 ;
    vector<Float_t> *sv_mu_LXYSig_ = 0 ;
    vector<Float_t> *sv_mu_LXYZSig_ = 0 ;
    vector<Float_t> *sv_mu_LX_ = 0 ;
    vector<Float_t> *sv_mu_LY_ = 0 ;
    vector<Float_t> *sv_mu_LZ_ = 0 ;
    vector<Float_t> *sv_mu_LXY_ = 0 ;
    vector<Float_t> *sv_mu_LXYZ_ = 0 ;
    vector<Float_t> *sv_mu_mass_ = 0 ;
    vector<Float_t> *sv_mu_eta_ = 0 ;
    vector<Float_t> *sv_mu_phi_ = 0 ;
    vector<Float_t> *sv_mu_pt_ = 0 ;
    vector<Float_t> *sv_mu_p_ = 0 ;
    vector<Float_t> *sv_mu_en_ = 0 ;
    vector<Float_t> *sv_mu_Beta_ = 0 ;
    vector<Float_t> *sv_mu_Gamma_ = 0 ;
    vector<Float_t> *sv_mu_CTau0_ = 0 ;
    vector<Float_t> *sv_mu_NDof_ = 0 ;
    vector<Float_t> *sv_mu_Chi2_ = 0 ;
    vector<Float_t> *sv_mu_Angle3D_ = 0 ;
    vector<Float_t> *sv_mu_Angle2D_ = 0 ;
    vector<Int_t>   *sv_mu_tracks_Sumcharge_ = 0 ;
    vector<Float_t> *sv_mu_tracks_Sumpt_ = 0 ;
    vector<Float_t> *sv_mu_match_dxyz_ = 0 ;
    
    vector<Float_t>   *jet_charge = 0;
    vector<Float_t>   *jet_et = 0;
    vector<Float_t>   *jet_pt = 0;
    vector<Float_t>   *jet_eta = 0;
    vector<Float_t>   *jet_phi = 0;
    vector<Float_t>   *jet_theta = 0;
    vector<Float_t>   *jet_en = 0;
    vector<Float_t>   *jet_chargedEmEnergy = 0;
    vector<Float_t>   *jet_neutralEmEnergyFraction = 0;
    vector<Float_t>   *jet_chargedHadronEnergy = 0;
    vector<Float_t>   *jet_neutralHadronEnergyFraction = 0;
    vector<Float_t>   *jet_chargedMuEnergy = 0;
    vector<Float_t>   *jet_chargedMuEnergyFraction = 0;
    vector<Float_t>   *jet_numberOfDaughters = 0;
    vector<Float_t>   *jet_muonEnergy = 0;
    vector<Float_t>   *jet_muonEnergyFraction = 0;
    vector<Float_t>   *jet_muonMultiplicity = 0;
    vector<Float_t>   *jet_neutralEmEnergy = 0;
    vector<Float_t>   *jet_neutralHadronEnergy = 0;
    vector<Float_t>   *jet_neutralHadronMultiplicity = 0;
    vector<Float_t>   *jet_neutralMultiplicity = 0;
    
    Float_t   pfMet_et;
    Float_t   pfMet_pt;
    Float_t   pfMet_px;
    Float_t   pfMet_py;
    Float_t   pfMet_pz;
    Float_t   pfMet_phi;
    Float_t   pfMet_en;
    Float_t   pfMet_sumEt;
    Float_t   caloMet_pt;
    Float_t   caloMet_phi;
    /*
      vector<Float_t>   *jet_btag_pt = 0;
      vector<Float_t>   *jet_btag_eta = 0;
      vector<Float_t>   *jet_btag_phi = 0;
      vector<Int_t>     *jet_btag_flavor = 0;
      vector<Float_t>   *jet_btag_pfDeepCSV_bb_discriminator = 0;
      vector<Float_t>   *jet_btag_pfDeepCSV_bbb_discriminator = 0;
      vector<Float_t>   *jet_btag_pfDeepCSV_bc_discriminator = 0;
    */
    Float_t transverseMass_ivf_;
    Float_t transverseMass_lep1_;
    Float_t transverseMass_ivfPluslep1_;
    
    Float_t sv_mass_corr_;
    
    Int_t mu_cmb_nDof_, mu_rpc_nDof_;
    Double_t mu_cmb_time_, mu_rpc_time_;
    Double_t mu_cmb_timeErr_, mu_rpc_timeErr_;
    
    oldtree->SetBranchAddress("mu_isTightMuon",&mu_isTightMuon_);
    oldtree->SetBranchAddress("mu_isLooseMuon",&mu_isLooseMuon_);
    oldtree->SetBranchAddress("mu_pt",&mu_pt_);
    oldtree->SetBranchAddress("mu_eta",&mu_eta_);
    oldtree->SetBranchAddress("mu_phi",&mu_phi_);
    oldtree->SetBranchAddress("mu_charge",&mu_charge_);
    oldtree->SetBranchAddress("mu_en",&mu_en_);
    oldtree->SetBranchAddress("mu_et",&mu_et_);
    oldtree->SetBranchAddress("mu_recoDeltaBeta",&deltaBeta_);
    oldtree->SetBranchAddress("mu_rhoIso", &mu_rhoIso_);
    oldtree->SetBranchAddress("mu_trackiso", &mu_trackiso_);
    oldtree->SetBranchAddress("mu_pfSumChargedHadronPt", &mu_pfSumChargedHadronPt_);
    oldtree->SetBranchAddress("mu_pfSumNeutralHadronEt", &mu_pfSumNeutralHadronEt_);
    oldtree->SetBranchAddress("mu_PFSumPhotonEt", &mu_PFSumPhotonEt_);
    oldtree->SetBranchAddress("mu_pfSumPUPt", &mu_pfSumPUPt_);
    oldtree->SetBranchAddress("mu_emIso", &mu_emIso_);
    oldtree->SetBranchAddress("mu_hadIso", &mu_hadIso_);
    //oldtree->SetBranchAddress("mu_normalizedChi2", &mu_normalizedChi2_);
    oldtree->SetBranchAddress("mu_dPToverPTTunePMuonBestTrack", &mu_dPToverPTTunePMuonBestTrack_);
    oldtree->SetBranchAddress("mu_absdxyTunePMuonBestTrack", &mu_absdxyTunePMuonBestTrack_);
    oldtree->SetBranchAddress("mu_absdxyErrorTunePMuonBestTrack", &mu_absdxyErrorTunePMuonBestTrack_);
    oldtree->SetBranchAddress("mu_absdxySigTunePMuonBestTrack", &mu_absdxySigTunePMuonBestTrack_);
    oldtree->SetBranchAddress("mu_absdzTunePMuonBestTrack", &mu_absdzTunePMuonBestTrack_);
    oldtree->SetBranchAddress("mu_absdzErrorTunePMuonBestTrack", &mu_absdzErrorTunePMuonBestTrack_);
    oldtree->SetBranchAddress("mu_absdzSigTunePMuonBestTrack", &mu_absdzSigTunePMuonBestTrack_);
    oldtree->SetBranchAddress("mu_recoiso", &mu_recoiso_);
    //oldtree->SetBranchAddress("mu_STATofDirection", &mu_STATofDirection_);
    oldtree->SetBranchAddress("mu_STATofNDof", &mu_STATofNDof_);
    oldtree->SetBranchAddress("mu_STATofTimeAtIpInOut", &mu_STATofTimeAtIpInOut_);
    oldtree->SetBranchAddress("mu_STATofTimeAtIpInOutErr", &mu_STATofTimeAtIpInOutErr_);
    oldtree->SetBranchAddress("mu_STATofTimeAtIpOutIn", &mu_STATofTimeAtIpOutIn_);
    oldtree->SetBranchAddress("mu_STATofTimeAtIpOutInErr", &mu_STATofTimeAtIpOutInErr_);
    //oldtree->SetBranchAddress("mu_numberOfValidMuonHits", &mu_numberOfValidMuonHits_);
    //oldtree->SetBranchAddress("mu_numberOfMatchedStations", &mu_numberOfMatchedStations_);
    //oldtree->SetBranchAddress("mu_numberOfValidPixelHits", &mu_numberOfValidPixelHits_);
    //oldtree->SetBranchAddress("mu_TrackQuality", &mu_TrackQuality_);
    //oldtree->SetBranchAddress("mu_InnerTrackQuality", &mu_InnerTrackQuality_);
    //    oldtree->SetBranchAddress("mu_FirstGenMatch", &mu_FirstGenMatch_);
    //oldtree->SetBranchAddress("mu_SecondGenMatch", &mu_SecondGenMatch_);
    
    //oldtree->SetBranchAddress("sv_tracks_pt", &sv_mu_TrackSize_);
    oldtree->SetBranchAddress("sv_LxySig", &sv_mu_LXYSig_);
    oldtree->SetBranchAddress("sv_LxyzSig", &sv_mu_LXYZSig_);
    oldtree->SetBranchAddress("sv_Lx", &sv_mu_LX_);
    oldtree->SetBranchAddress("sv_Ly", &sv_mu_LY_);
    oldtree->SetBranchAddress("sv_Lz", &sv_mu_LZ_);
    oldtree->SetBranchAddress("sv_Lxy", &sv_mu_LXY_);
    oldtree->SetBranchAddress("sv_Lxyz", &sv_mu_LXYZ_);
    oldtree->SetBranchAddress("sv_mass", &sv_mu_mass_);
    oldtree->SetBranchAddress("sv_eta", &sv_mu_eta_);   
    oldtree->SetBranchAddress("sv_phi", &sv_mu_phi_);  
    oldtree->SetBranchAddress("sv_pt", &sv_mu_pt_);  
    oldtree->SetBranchAddress("sv_p", &sv_mu_p_);  
    oldtree->SetBranchAddress("sv_energy", &sv_mu_en_);
    oldtree->SetBranchAddress("sv_Beta", &sv_mu_Beta_);
    oldtree->SetBranchAddress("sv_Gamma", &sv_mu_Gamma_); 
    oldtree->SetBranchAddress("sv_CTau0", &sv_mu_CTau0_);  
    oldtree->SetBranchAddress("sv_NDof", &sv_mu_NDof_);   
    oldtree->SetBranchAddress("sv_Chi2", &sv_mu_Chi2_);   
    oldtree->SetBranchAddress("sv_Angle3D", &sv_mu_Angle3D_);
    oldtree->SetBranchAddress("sv_Angle2D", &sv_mu_Angle2D_);  
    oldtree->SetBranchAddress("sv_tracks_Sumcharge", &sv_mu_tracks_Sumcharge_); 
    oldtree->SetBranchAddress("sv_tracks_Sumpt", &sv_mu_tracks_Sumpt_);  
    //oldtree->SetBranchAddress("sv_match_dxyz", &sv_mu_match_dxyz_);
    /*
      oldtree->SetBranchAddress("jet_btag_pt", &jet_btag_pt);
      oldtree->SetBranchAddress("jet_btag_eta", &jet_btag_eta);
      oldtree->SetBranchAddress("jet_btag_phi", &jet_btag_phi);
      oldtree->SetBranchAddress("jet_btag_flavor", &jet_btag_flavor);
      oldtree->SetBranchAddress("pfDeepCSV_bbb", &jet_btag_pfDeepCSV_bbb_discriminator);
      oldtree->SetBranchAddress("pfDeepCSV_bb", &jet_btag_pfDeepCSV_bb_discriminator);
      oldtree->SetBranchAddress("pfDeepCSV_bc", &jet_btag_pfDeepCSV_bc_discriminator);
    */
    oldtree->SetBranchAddress("pfMet_et", &pfMet_et);
    oldtree->SetBranchAddress("pfMet_pt", &pfMet_pt);
    oldtree->SetBranchAddress("pfMet_phi", &pfMet_phi);
    oldtree->SetBranchAddress("pfMet_en", &pfMet_en);
    oldtree->SetBranchAddress("pfMet_sumEt", &pfMet_sumEt);
    oldtree->SetBranchAddress("caloMet_pt", &caloMet_pt);
    oldtree->SetBranchAddress("caloMet_phi", &caloMet_phi);
    oldtree->SetBranchAddress("pfMet_px", &pfMet_px);
    oldtree->SetBranchAddress("pfMet_py", &pfMet_py);
    oldtree->SetBranchAddress("pfMet_pz", &pfMet_pz);    

    oldtree->SetBranchAddress("jet_charge", &jet_charge);
    oldtree->SetBranchAddress("jet_et", &jet_et);
    oldtree->SetBranchAddress("jet_pt", &jet_pt);
    oldtree->SetBranchAddress("jet_eta", &jet_eta);
    oldtree->SetBranchAddress("jet_phi", &jet_phi);
    oldtree->SetBranchAddress("jet_theta", &jet_theta);
    oldtree->SetBranchAddress("jet_en", &jet_en);
    oldtree->SetBranchAddress("jet_chargedEmEnergy", &jet_chargedEmEnergy);
    oldtree->SetBranchAddress("jet_neutralEmEnergyFraction", &jet_neutralEmEnergyFraction);
    oldtree->SetBranchAddress("jet_chargedHadronEnergy", &jet_chargedHadronEnergy);
    oldtree->SetBranchAddress("jet_neutralHadronEnergyFraction", &jet_neutralHadronEnergyFraction);
    oldtree->SetBranchAddress("jet_chargedMuEnergy", &jet_chargedMuEnergy);
    oldtree->SetBranchAddress("jet_chargedMuEnergyFraction", &jet_chargedMuEnergyFraction);
    oldtree->SetBranchAddress("jet_numberOfDaughters", &jet_numberOfDaughters);
    oldtree->SetBranchAddress("jet_muonEnergy", &jet_muonEnergy);
    oldtree->SetBranchAddress("jet_muonEnergyFraction", &jet_muonEnergyFraction);
    oldtree->SetBranchAddress("jet_muonMultiplicity", &jet_muonMultiplicity);
    oldtree->SetBranchAddress("jet_neutralEmEnergy", &jet_neutralEmEnergy);
    oldtree->SetBranchAddress("jet_neutralHadronEnergy", &jet_neutralHadronEnergy);
    oldtree->SetBranchAddress("jet_neutralHadronMultiplicity", &jet_neutralHadronMultiplicity);
    oldtree->SetBranchAddress("jet_neutralMultiplicity", &jet_neutralMultiplicity);
    
    
    //transverse mass
    //oldtree->SetBranchAddress("tranvsverseMass_ivf", &transverseMass_ivf_);
    //oldtree->SetBranchAddress("tranvsverseMass_lep1", &transverseMass_lep1_);
    //oldtree->SetBranchAddress("tranvsverseMass_ivfPluslep1", &transverseMass_ivfPluslep1_);
    
    //sv mass correction
    //oldtree->SetBranchAddress("sv_mass_corr", &sv_mass_corr_);
    
    
    //RPC info
    
    oldtree->SetBranchAddress("mu_cmb_nDof", &mu_cmb_nDof_);
    oldtree->SetBranchAddress("mu_rpc_nDof", &mu_rpc_nDof_);
    oldtree->SetBranchAddress("mu_cmb_time", &mu_cmb_time_);
    oldtree->SetBranchAddress("mu_rpc_time", &mu_rpc_time_);
    oldtree->SetBranchAddress("mu_cmb_timeErr", &mu_cmb_timeErr_);
    oldtree->SetBranchAddress("mu_rpc_timeErr", &mu_rpc_timeErr_);
    
    
    //Throw away some branches -- thos are all the  one I don't want to keep in the ntuples
    //IMPORTANT: we cannot throw away the branches with the variable we are using in the loop!
    oldtree->SetBranchStatus("ele_*",    0);
    
    //================================ PU Weight ============================================// 
    //Float_t pu_weight;
    
    //TBranch* branch_pu_weight = newtree->Branch("pu_weight",&pu_weight,"pu_weight/F");
    //======================= First Muon Variables ==========================================//   
    Float_t mu_promptPt_, mu_promptEta_, mu_promptPhi_, mu_promptCharge_, mu_promptEt_, mu_promptE_;
    
    Float_t mu_promptRhoIso_,        mu_promptTrackiso_,       mu_promptPfSumChHadPt_,
      mu_promptPFSumPhotonEt_,  mu_promptPfSumPUPt_,      mu_promptEmIso_,
      mu_promptHadIso_,        mu_promptDPToverPT_,
      mu_promptAbsdxy_,        mu_promptAbsdxyError_,    mu_promptAbsdxySig_,
      mu_promptAbsdz_,         mu_promptAbsdzError_,     mu_promptAbsdzSig_,
      mu_promptRecoDeltaBeta_, mu_promptRecoiso_,  
      mu_promptNDof_,          mu_promptTimeAtIpInOut_,  mu_promptTimeAtIpInOutErr_,
      mu_promptTimeAtIpOutIn_, mu_promptTimeAtIpOutInErr_,  mu_promptPfSumNHadEt_;
    
    Int_t   mu_promptGenMatch_;
    
    
    TBranch* branch_mu_promptPt = newtree->Branch("mu_promptPt",&mu_promptPt_,"mu_promptPt/F");
    TBranch* branch_mu_promptEta = newtree->Branch("mu_promptEta",&mu_promptEta_,"mu_promptEta/F");
    TBranch* branch_mu_promptPhi = newtree->Branch("mu_promptPhi",&mu_promptPhi_,"mu_promptPhi/F");
    TBranch* branch_mu_promptCharge = newtree->Branch("mu_promptCharge",&mu_promptCharge_,"mu_promptCharge/F");
    TBranch* branch_mu_promptE = newtree->Branch("mu_promptE",&mu_promptE_,"mu_promptE/F");
    TBranch* branch_mu_promptEt = newtree->Branch("mu_promptEt",&mu_promptEt_,"mu_promptEt/F");
    TBranch* branch_mu_promptRhoIso = newtree->Branch("mu_promptRhoIso", &mu_promptRhoIso_, "mu_promptRhoIso/F");
    TBranch* branch_mu_promptTrackiso = newtree->Branch("mu_promptTrackiso", &mu_promptTrackiso_, "mu_promptTrackiso/F");
    TBranch* branch_mu_promptPfSumChHadPt = newtree->Branch("mu_promptPfSumChHadPt", &mu_promptPfSumChHadPt_, "mu_promptPfSumChHadPt/F");
    TBranch* branch_mu_promptPfSumNHadEt = newtree->Branch("mu_promptPfSumNHadEt", &mu_promptPfSumNHadEt_, "mu_promptPfSumNHadEt/F");
    TBranch* branch_mu_promptPFSumPhotonEt = newtree->Branch("mu_promptPFSumPhotonEt", &mu_promptPFSumPhotonEt_, "mu_promptPFSumPhotonEt/F");
    TBranch* branch_mu_promptPfSumPUPt = newtree->Branch("mu_promptPfSumPUPt", &mu_promptPfSumPUPt_, "mu_promptPfSumPUPt/F");
    TBranch* branch_mu_promptEmIso = newtree->Branch("mu_promptEmIso", &mu_promptEmIso_, "mu_promptEmIso/F");
    TBranch* branch_mu_promptHadIso = newtree->Branch("mu_promptHadIso", &mu_promptHadIso_, "mu_promptHadIso/F");
    //TBranch* branch_mu_promptNormalizedChi2 = newtree->Branch("mu_promptNormalizedChi2", &mu_promptNormalizedChi2_, "mu_promptNormalizedChi2/F");
    TBranch* branch_mu_promptDPToverPT = newtree->Branch("mu_promptDPToverPT", &mu_promptDPToverPT_, "mu_promptDPToverPT/F");
    TBranch* branch_mu_promptAbsdxy = newtree->Branch("mu_promptAbsdxy", &mu_promptAbsdxy_, "mu_promptAbsdxy/F");
    TBranch* branch_mu_promptAbsdxyError = newtree->Branch("mu_promptAbsdxyError", &mu_promptAbsdxyError_, "mu_promptAbsdxyError/F");
    TBranch* branch_mu_promptAbsdxySig = newtree->Branch("mu_promptAbsdxySig", &mu_promptAbsdxySig_, "mu_promptAbsdxySig/F");
    TBranch* branch_mu_promptAbsdz = newtree->Branch("mu_promptAbsdz", &mu_promptAbsdz_, "mu_promptAbsdz/F");
    TBranch* branch_mu_promptAbsdzError = newtree->Branch("mu_promptAbsdzError", &mu_promptAbsdzError_, "mu_promptAbsdzError/F");
    TBranch* branch_mu_promptAbsdzSig = newtree->Branch("mu_promptAbsdzSig", &mu_promptAbsdz_, "mu_promptAbsdzSig/F");
    TBranch* branch_mu_promptRecoDeltaBeta = newtree->Branch("mu_promptRecoDeltaBeta", &mu_promptRecoDeltaBeta_, "mu_promptRecoDeltaBeta/F");
    TBranch* branch_mu_promptRecoiso = newtree->Branch("mu_promptRecoiso", &mu_promptRecoiso_, "mu_promptRecoiso/F");
    //TBranch* branch_mu_promptDirection = newtree->Branch("mu_promptDirection", &mu_promptDirection_, "mu_promptDirection/F");
    TBranch* branch_mu_promptNDof = newtree->Branch("mu_promptNDof", &mu_promptNDof_, "mu_promptNDof/F");
    TBranch* branch_mu_promptTimeAtIpInOut = newtree->Branch("mu_promptTimeAtIpInOut", &mu_promptTimeAtIpInOut_, "mu_promptTimeAtIpInOut/F");
    TBranch* branch_mu_promptTimeAtIpInOutErr = newtree->Branch("mu_promptTimeAtIpInOutErr", &mu_promptTimeAtIpInOutErr_, "mu_promptTimeAtIpInOutErr/F");
    TBranch* branch_mu_promptTimeAtIpOutIn = newtree->Branch("mu_promptTimeAtIpOutIn", &mu_promptTimeAtIpOutIn_, "mu_promptTimeAtIpOutIn/F");
    TBranch* branch_mu_promptTimeAtIpOutInErr = newtree->Branch("mu_promptTimeAtIpOutInErr", &mu_promptTimeAtIpOutInErr_, "mu_promptTimeAtIpOutInErr/F");
    //TBranch* branch_mu_promptMatchedStations = newtree->Branch("mu_promptMatchedStations", &mu_promptMatchedStations_, "mu_promptMatchedStations/I");
    //TBranch* branch_mu_promptValidPixelHits = newtree->Branch("mu_promptValidPixelHits", &mu_promptValidPixelHits_, "mu_promptValidPixelHits/I");
    //TBranch* branch_mu_promptTrackQuality = newtree->Branch("mu_promptTrackQuality", &mu_promptTrackQuality_, "mu_promptTrackQuality/I");
    //TBranch* branch_mu_promptInrTrackQuality = newtree->Branch("mu_promptInrTrackQuality", &mu_promptInrTrackQuality_, "mu_promptInrTrackQuality/I");
    //TBranch* branch_mu_promptValidMuonHits = newtree->Branch("mu_promptValidMuonHits", &mu_promptValidMuonHits_, "mu_promptValidMuonHits/I");
    TBranch* branch_mu_promptGenMatch = newtree->Branch("mu_promptGenMatch", &mu_promptGenMatch_, "mu_promptGenMatch/I");
    
    
    //======================= Second Muon Variables ==========================================//
    Float_t mu_secondPt_, mu_secondEta_, mu_secondPhi_, mu_secondCharge_, mu_secondEt_, mu_secondE_;
    
    Float_t mu_secondRhoIso_,        mu_secondTrackiso_,  mu_secondPfSumChHadPt_,
      mu_secondPFSumPhotonEt_, mu_secondPfSumPUPt_,      mu_secondEmIso_,        
      mu_secondHadIso_,         mu_secondDPToverPT_,   
      mu_secondAbsdxy_,        mu_secondAbsdxyError_,    mu_secondAbsdxySig_, 
      mu_secondAbsdz_,         mu_secondAbsdzError_,     mu_secondAbsdzSig_,
      mu_secondRecoDeltaBeta_, mu_secondRecoiso_,  
      mu_secondNDof_,          mu_secondTimeAtIpInOut_,  mu_secondTimeAtIpInOutErr_, 
      mu_secondTimeAtIpOutIn_, mu_secondTimeAtIpOutInErr_,  mu_secondPfSumNHadEt_;
    
    Float_t mu_DeltaBetaR3_, mu_DiMuMass_, mu_Size_, mu_DeltaR_;
    
    Int_t    mu_secondGenMatch_;
    
    TBranch* branch_mu_secondPt  = newtree->Branch("mu_secondPt",&mu_secondPt_,"mu_secondPt/F");
    TBranch* branch_mu_secondEta = newtree->Branch("mu_secondEta",&mu_secondEta_,"mu_secondEta/F");
    TBranch* branch_mu_secondPhi = newtree->Branch("mu_secondPhi",&mu_secondPhi_,"mu_secondPhi/F");
    TBranch* branch_mu_secondCharge = newtree->Branch("mu_secondCharge",&mu_secondCharge_,"mu_secondCharge/F");
    TBranch* branch_mu_secondE  = newtree->Branch("mu_secondE",&mu_secondE_,"mu_secondE/F");
    TBranch* branch_mu_secondEt = newtree->Branch("mu_secondEt",&mu_secondEt_,"mu_secondEt/F");
    TBranch* branch_mu_secondRhoIso   = newtree->Branch("mu_secondRhoIso", &mu_secondRhoIso_, "mu_secondRhoIso/F");
    TBranch* branch_mu_secondTrackiso = newtree->Branch("mu_secondTrackiso", &mu_secondTrackiso_, "mu_secondTrackiso/F");
    TBranch* branch_mu_secondPfSumChHadPt  = newtree->Branch("mu_secondPfSumChHadPt", &mu_secondPfSumChHadPt_, "mu_secondPfSumChHadPt/F");
    TBranch* branch_mu_secondPfSumNHadEt   = newtree->Branch("mu_secondPfSumNHadEt", &mu_secondPfSumNHadEt_, "mu_secondPfSumNHadEt/F");
    TBranch* branch_mu_secondPFSumPhotonEt = newtree->Branch("mu_secondPFSumPhotonEt", &mu_secondPFSumPhotonEt_, "mu_secondPFSumPhotonEt/F");
    TBranch* branch_mu_secondPfSumPUPt = newtree->Branch("mu_secondPfSumPUPt", &mu_secondPfSumPUPt_, "mu_secondPfSumPUPt/F");
    TBranch* branch_mu_secondEmIso  = newtree->Branch("mu_secondEmIso", &mu_secondEmIso_, "mu_secondEmIso/F");
    TBranch* branch_mu_secondHadIso = newtree->Branch("mu_secondHadIso", &mu_secondHadIso_, "mu_secondHadIso/F");
    //TBranch* branch_mu_secondNormalizedChi2 = newtree->Branch("mu_secondNormalizedChi2", &mu_secondNormalizedChi2_, "mu_secondNormalizedChi2/F");
    TBranch* branch_mu_secondDPToverPT = newtree->Branch("mu_secondDPToverPT", &mu_secondDPToverPT_, "mu_secondDPToverPT/F");
    TBranch* branch_mu_secondAbsdxy = newtree->Branch("mu_secondAbsdxy", &mu_secondAbsdxy_, "mu_secondAbsdxy/F");
    TBranch* branch_mu_secondAbsdxyError = newtree->Branch("mu_secondAbsdxyError", &mu_secondAbsdxyError_, "mu_secondAbsdxyError/F");
    TBranch* branch_mu_secondAbsdxySig   = newtree->Branch("mu_secondAbsdxySig", &mu_secondAbsdxySig_, "mu_secondAbsdxySig/F");
    TBranch* branch_mu_secondAbsdz       = newtree->Branch("mu_secondAbsdz", &mu_secondAbsdz_, "mu_secondAbsdz/F");
    TBranch* branch_mu_secondAbsdzError  = newtree->Branch("mu_secondAbsdzError", &mu_secondAbsdzError_, "mu_secondAbsdzError/F");
    TBranch* branch_mu_secondAbsdzSig = newtree->Branch("mu_secondAbsdzSig", &mu_secondAbsdz_, "mu_secondAbsdzSig/F");
    TBranch* branch_mu_secondRecoDeltaBeta = newtree->Branch("mu_secondRecoDeltaBeta", &mu_secondRecoDeltaBeta_, "mu_secondRecoDeltaBeta/F");
    TBranch* branch_mu_secondRecoiso   = newtree->Branch("mu_secondRecoiso", &mu_secondRecoiso_, "mu_secondRecoiso/F");
    //TBranch* branch_mu_secondDirection = newtree->Branch("mu_secondDirection", &mu_secondDirection_, "mu_secondDirection/F");
    TBranch* branch_mu_secondNDof = newtree->Branch("mu_secondNDof", &mu_secondNDof_, "mu_secondNDof/F");
    TBranch* branch_mu_secondTimeAtIpInOut    = newtree->Branch("mu_secondTimeAtIpInOut", &mu_secondTimeAtIpInOut_, "mu_secondTimeAtIpInOut/F");
    TBranch* branch_mu_secondTimeAtIpInOutErr = newtree->Branch("mu_secondTimeAtIpInOutErr", &mu_secondTimeAtIpInOutErr_, "mu_secondTimeAtIpInOutErr/F");
    TBranch* branch_mu_secondTimeAtIpOutIn = newtree->Branch("mu_secondTimeAtIpOutIn", &mu_secondTimeAtIpOutIn_, "mu_secondTimeAtIpOutIn/F");
    TBranch* branch_mu_secondTimeAtIpOutInErr = newtree->Branch("mu_secondTimeAtIpOutInErr", &mu_secondTimeAtIpOutInErr_, "mu_secondTimeAtIpOutInErr/F");
    //TBranch* branch_mu_secondMatchedStations  = newtree->Branch("mu_secondMatchedStations", &mu_secondMatchedStations_, "mu_secondMatchedStations/I");
    //TBranch* branch_mu_secondValidPixelHits   = newtree->Branch("mu_secondValidPixelHits", &mu_secondValidPixelHits_, "mu_secondValidPixelHits/I");
    //TBranch* branch_mu_secondTrackQuality    = newtree->Branch("mu_secondTrackQuality", &mu_secondTrackQuality_, "mu_secondTrackQuality/I");
    //TBranch* branch_mu_secondInrTrackQuality = newtree->Branch("mu_secondInrTrackQuality", &mu_secondInrTrackQuality_, "mu_secondInrTrackQuality/I");
    //TBranch* branch_mu_secondValidMuonHits   = newtree->Branch("mu_secondValidMuonHits", &mu_secondValidMuonHits_, "mu_secondValidMuonHits/I");
    TBranch* branch_mu_secondGenMatch        = newtree->Branch("mu_secondGenMatch", &mu_secondGenMatch_, "mu_secondGenMatch/I");
    
    TBranch* branch_mu_DeltaBetaR3 = newtree->Branch("mu_DeltaBetaR3", &mu_DeltaBetaR3_, "mu_DeltaBetaR3/F");
    TBranch* branch_mu_DiMuMass    = newtree->Branch("mu_DiMuMass", &mu_DiMuMass_, "mu_DiMuMass/F");
    TBranch* branch_mu_DeltaR      = newtree->Branch("mu_DeltaR", &mu_DeltaR_, "mu_DeltaR/F");
    TBranch* branch_mu_Size        = newtree->Branch("mu_Size", &mu_Size_, "mu_Size/F");
    
    //======================= Second Vertex Variables ==========================================//
    Float_t  sv_LXYSig_new_,  sv_LXYZSig_new_, sv_LXY_new_, sv_LX_new_,sv_LY_new_,sv_LZ_new_, sv_LXYZ_new_, sv_mass_new_, 
      sv_eta_new_,     sv_phi_new_,     sv_pt_new_,   sv_p_new_,  sv_en_new_,  sv_Beta_new_, 
      sv_Gamma_new_,   sv_CTau0_new_,   sv_NDof_new_, sv_Chi2_new_, sv_Angle3D_new_, 
      sv_Angle2D_new_, sv_tracks_Sumpt_new_, sv_match_dxyz_new_; 
    
    
    Int_t sv_TrackSize_new_, sv_tracks_Sumcharge_new_;
    
    TBranch* branch_sv_TrackSize  = newtree->Branch("sv_TrackSize", &sv_TrackSize_new_, "sv_TrackSize/I");
    TBranch* branch_sv_LXYSig     = newtree->Branch("sv_LXYSig", &sv_LXYSig_new_, "sv_LXYSig/F");
    TBranch* branch_sv_LXYZSig    = newtree->Branch("sv_LXYZSig", &sv_LXYZSig_new_, "sv_LXYZSig/F");
    TBranch* branch_sv_LX        = newtree->Branch("sv_LX", &sv_LX_new_, "sv_LX/F");
    TBranch* branch_sv_LY        = newtree->Branch("sv_LY", &sv_LY_new_, "sv_LY/F");
    TBranch* branch_sv_LZ        = newtree->Branch("sv_LZ", &sv_LZ_new_, "sv_LZ/F");
    TBranch* branch_sv_LXY        = newtree->Branch("sv_LXY", &sv_LXY_new_, "sv_LXY/F");
    TBranch* branch_sv_LXYZ       = newtree->Branch("sv_LXYZ", &sv_LXYZ_new_, "sv_LXYZ/F");
    TBranch* branch_sv_mass       = newtree->Branch("sv_mass", &sv_mass_new_, "sv_mass/F");
    TBranch* branch_sv_eta        = newtree->Branch("sv_eta", &sv_eta_new_, "sv_eta/F");
    TBranch* branch_sv_phi        = newtree->Branch("sv_phi", &sv_phi_new_, "sv_phi/F");
    TBranch* branch_sv_pt         = newtree->Branch("sv_pt", &sv_pt_new_, "sv_pt/F");
    TBranch* branch_sv_p          = newtree->Branch("sv_p", &sv_p_new_, "sv_p/F");
    TBranch* branch_sv_en         = newtree->Branch("sv_en", &sv_en_new_, "sv_en/F");
    TBranch* branch_sv_Beta       = newtree->Branch("sv_Beta", &sv_Beta_new_, "sv_Beta/F");
    TBranch* branch_sv_Gamma      = newtree->Branch("sv_Gamma", &sv_Gamma_new_, "sv_Gamma/F");
    TBranch* branch_sv_CTau0      = newtree->Branch("sv_CTau0", &sv_CTau0_new_, "sv_CTau0/F");
    TBranch* branch_sv_NDof       = newtree->Branch("sv_NDof", &sv_NDof_new_, "sv_NDof/F");
    TBranch* branch_sv_Chi2       = newtree->Branch("sv_Chi2_mew", &sv_Chi2_new_, "sv_Chi2/F");
    TBranch* branch_sv_Angle3D    = newtree->Branch("sv_Angle3D", &sv_Angle3D_new_, "sv_Angle3D/F");
    TBranch* branch_sv_Angle2D    = newtree->Branch("sv_Angle2D", &sv_Angle2D_new_, "sv_Angle2D/F");
    TBranch* branch_sv_tracks_Sumcharge = newtree->Branch("sv_tracks_Sumcharge", &sv_tracks_Sumcharge_new_, "sv_tracks_Sumcharge/I");
    TBranch* branch_sv_tracks_Sumpt     = newtree->Branch("sv_tracks_Sumpt", &sv_tracks_Sumpt_new_, "sv_tracks_Sumpt/F");
    TBranch* branch_sv_match_dxyz            = newtree->Branch("sv_match_dxyz_", &sv_match_dxyz_new_, "sv_match_dxyz/F");
    //================================= Jets Variables ==============================================//
    Float_t jet_charge_, jet_et_,    jet_pt_, jet_eta_, jet_NEmEnFraction_,
      jet_phi_,    jet_theta_, jet_en_, jet_chEmEn_, jet_NHadEnFraction_,
      jet_chMuEn_,         jet_chMuEnFraction_, jet_numberOfDaughters_, 
      jet_muonEnergyFraction_, jet_muonMultiplicity_, 
      jet_neutralEmEnergy_,  jet_neutralHadronEnergy_, jet_NHadMultiplicity_, 
      jet_NMultiplicity_, jet_chHadEn_, jet_muonEnergy_;
    
    TBranch* branch_jet_charge_  = newtree->Branch("jet_charge_", &jet_charge_, "jet_charge_/F");
    TBranch* branch_jet_et_  = newtree->Branch("jet_et_", &jet_et_, "jet_et_/F");
    TBranch* branch_jet_pt_  = newtree->Branch("jet_pt_", &jet_pt_, "jet_pt_/F");
    TBranch* branch_jet_eta_  = newtree->Branch("jet_eta_", &jet_eta_, "jet_eta_/F");
    TBranch* branch_jet_phi_  = newtree->Branch("jet_phi_", &jet_phi_, "jet_phi_/F");
    TBranch* branch_jet_theta_  = newtree->Branch("jet_theta_", &jet_theta_, "jet_theta_/F");
    TBranch* branch_jet_en_  = newtree->Branch("jet_en_", &jet_en_, "jet_en_/F");
    TBranch* branch_jet_chargedEmEnergy_  = newtree->Branch("jet_chEmEn_", &jet_chEmEn_, "jet_chEmEn_/F");
    TBranch* branch_jet_NEmEnFraction_  = newtree->Branch("jet_NEmEnFraction_", &jet_NEmEnFraction_, "jet_NEmEnFraction_/F");
    TBranch* branch_jet_chHadEn_  = newtree->Branch("jet_chHadEn_", &jet_chHadEn_, "jet_chHadEn_/F");
    TBranch* branch_jet_NHadEnFraction_  = newtree->Branch("jet_NHadEnFraction_", &jet_NHadEnFraction_, "jet_NHadEnFraction_/F");
    TBranch* branch_jet_chMuEn_  = newtree->Branch("jet_chMuEn_", &jet_chMuEn_, "jet_chMuEn_/F");
    TBranch* branch_jet_chMuEnFraction_  = newtree->Branch("jet_chMuEnFraction_", &jet_chMuEnFraction_, "jet_chMuEnFraction_/F");
    TBranch* branch_jet_numberOfDaughters_  = newtree->Branch("jet_numberOfDaughters_", &jet_numberOfDaughters_, "jet_numberOfDaughters_/F");
    TBranch* branch_jet_muonEnergy_  = newtree->Branch("jet_muonEnergy_", &jet_muonEnergy_, "jet_muonEnergy_/F");
    TBranch* branch_jet_muonEnergyFraction_  = newtree->Branch("jet_muonEnergyFraction_", &jet_muonEnergyFraction_, "jet_muonEnergyFraction_/F");
    TBranch* branch_jet_muonMultiplicity_  = newtree->Branch("jet_muonMultiplicity_", &jet_muonMultiplicity_, "jet_muonMultiplicity_/F");
    TBranch* branch_jet_neutralEmEnergy_  = newtree->Branch("jet_neutralEmEnergy_", &jet_neutralEmEnergy_, "jet_neutralEmEnergy_/F");
    TBranch* branch_jet_neutralHadronEnergy_  = newtree->Branch("jet_neutralHadronEnergy_", &jet_neutralHadronEnergy_, "jet_neutralHadronEnergy_/F");
    TBranch* branch_jet_NHadMultiplicity_  = newtree->Branch("jet_NHadMultiplicity_", &jet_NHadMultiplicity_, "jet_NHadMultiplicity_/F");
    TBranch* branch_jet_NMultiplicity_  = newtree->Branch("jet_NMultiplicity_", &jet_NMultiplicity_, "jet_NMultiplicity_/F");
    
    //============================ Missing Energy Variables =========================================//
    Float_t pfMet_et_, pfMet_pt_, pfMet_phi_, pfMet_en_, pfMet_sumEt_, caloMet_phi_,caloMet_pt_;
    
    TBranch* branch_pfMet_et_    = newtree->Branch("pfMet_et_", &pfMet_et_, "pfMet_et_/F");
    TBranch* branch_pfMet_pt_    = newtree->Branch("pfMet_pt_", &pfMet_pt_, "pfMet_pt_/F");
    TBranch* branch_pfMet_phi_   = newtree->Branch("pfMet_phi_", &pfMet_phi_, "pfMet_phi_/F");
    TBranch* branch_pfMet_en_    = newtree->Branch("pfMet_en_", &pfMet_en_, "pfMet_en_/F");
    TBranch* branch_pfMet_sumEt_ = newtree->Branch("pfMet_sumEt_", &pfMet_sumEt_, "pfMet_sumEt_/F");
    TBranch* branch_caloMet_pt_  = newtree->Branch("caloMet_pt_", &caloMet_pt_, "caloMet_pt_/F");
    TBranch* branch_caloMet_phi_ = newtree->Branch("caloMet_phi_", &caloMet_phi_, "caloMet_phi_/F");
    
    //=============================== b-tagging Variables ===========================================//
    Float_t jet_btag_pt_, jet_btag_eta_,  jet_btag_phi_, jet_btag_bbb_discriminator_, jet_btag_bb_discriminator_, jet_btag_bc_discriminator_;  
    Int_t   jet_btag_flavor_;
    
    TBranch* branch_jet_btag_pt_  = newtree->Branch("jet_btag_pt_", &jet_btag_pt_, "jet_btag_pt_/F");
    TBranch* branch_jet_btag_eta_  = newtree->Branch("jet_btag_eta_", &jet_btag_eta_, "jet_btag_eta_/F");
    TBranch* branch_jet_btag_phi_  = newtree->Branch("jet_btag_phi_", &jet_btag_phi_, "jet_btag_phi_/F");
    TBranch* branch_jet_btag_flavor_  = newtree->Branch("jet_btag_flavor_", &jet_btag_flavor_, "jet_btag_flavor_/F");
    TBranch* branch_jet_btag_bbb_discriminator_  = newtree->Branch("jet_btag_bbb_discriminator_", &jet_btag_bbb_discriminator_, "jet_btag_bbb_discriminator_/F");
    TBranch* branch_jet_btag_bb_discriminator_  = newtree->Branch("jet_btag_bb_discriminator_", &jet_btag_bb_discriminator_, "jet_btag_bb_discriminator_/F");
    TBranch* branch_jet_btag_bc_discriminator_  = newtree->Branch("jet_btag_bc_discriminator_", &jet_btag_bc_discriminator_, "jet_btag_bc_discriminator_/F");
    
    
    //=============================== transverse mass Variables ===========================================//
    Float_t trMass_ivf_, trMass_lep1_, trMass_ivfPluslep1_;
    
    TBranch* branch_transverseMass_ivf_  = newtree->Branch("transverseMass_ivf", &trMass_ivf_, "trMass_ivf_/F");
    TBranch* branch_transverseMass_lep1_  = newtree->Branch("transverseMass_lep1", &trMass_lep1_, "trMass_lep1_/F");
    TBranch* branch_transverseMass_ivfPluslep1_  = newtree->Branch("transverseMass_ivfPluslep1", &trMass_ivfPluslep1_, "trMass_ivfPluslep1_/F");
    
    //=============================== sv mass correction ===========================================//
    Float_t sv_mass_correction_;
    
    TBranch* branch_sv_mass_correction_ = newtree->Branch("sv_mass_correction", &sv_mass_correction_, "sv_mass_correction_/F");

    
    //=============================== rpc variables ===========================================//
    Int_t mu_cmb_nDof_newtree_, mu_rpc_nDof_newtree_;
    Double_t mu_cmb_time_newtree_, mu_rpc_time_newtree_;
    Double_t mu_cmb_timeErr_newtree_,mu_rpc_timeErr_newtree_;
    
    TBranch* branch_rpc_nDof_ = newtree->Branch("rpc_nDof_", &mu_rpc_nDof_newtree_, "mu_rpc_nDof_newtree_/I");
    TBranch* branch_rpc_time_ = newtree->Branch("rpc_time_", &mu_rpc_time_newtree_, "mu_rpc_time_newtree_/D");
    TBranch* branch_rpc_timeErr_ = newtree->Branch("rpc_timeErr_", &mu_rpc_timeErr_newtree_, "mu_rpc_timeErr_newtree_/D");
    TBranch* branch_cmb_nDof_ = newtree->Branch("cmb_nDof_", &mu_cmb_nDof_newtree_, "mu_cmb_nDof_newtree_/I");
    TBranch* branch_cmb_time_ = newtree->Branch("cmb_time_", &mu_cmb_time_newtree_, "mu_cmb_time_newtree_/D");
    TBranch* branch_cmb_timeErr_ = newtree->Branch("cmb_timeErr_", &mu_cmb_timeErr_newtree_, "mu_cmb_timeErr_newtree_/D");
    
    
    
    
    //======================= Start the running over input branches ==========================================//
    for (int i=0;i<oldtree->GetEntriesFast(); i++) {
      if (i%10000==0) cout<<i<<endl;
      oldtree->GetEntry(i);
      
      if (passIsoMu24==0) continue;  // cut on the trigger!
      
      /*
	unsigned pu = -1;
	for(unsigned i=0; i<PU_Weight->size(); ++i){
	pu = i;
	}*/
      
      Float_t   minPt_prompt = -1000;
      unsigned  FirstMuon = -1;
      for(unsigned i=0; i<mu_isTightMuon_->size(); i++){
	if (mu_isTightMuon_->at(i)==0. || deltaBeta_->at(i)>isoCut || mu_pt_->at(i) < 24 || abs(mu_eta_->at(i)) > 2.4) continue;  
	if (mu_pt_->at(i) > minPt_prompt){
	  minPt_prompt=mu_pt_->at(i);
	  FirstMuon=i;	  
	}
      }	
      //here I save only the info of the prompt muon -- for the sv muon we can clone the branch as it is..
      Float_t   minPt_second = -1000;
      unsigned  SecondMuon = -1;
      int count = 0;
      for(unsigned i=0; i<mu_isLooseMuon_->size(); i++){
	if(i == FirstMuon) continue;
	if(mu_isLooseMuon_->size()== 1) continue;
	if(mu_isLooseMuon_->at(i)==0 || mu_pt_->at(i) < 5  || abs(mu_eta_->at(i)) > 2.4) continue;
	count++;
	if (mu_pt_->at(i)>minPt_second){
	minPt_second=mu_pt_->at(i);
	SecondMuon=i;
	}
      }    
      unsigned  SecondVertex = -1;
      for(unsigned i=0; i<sv_mu_pt_->size(); i++){
	if(sv_mu_pt_->at(i) <  0) continue;
	SecondVertex=i;
      }
      
      if(FirstMuon != -1 && SecondMuon != -1 && SecondVertex != -1){
	
	//pile up weight
      //pu_weight = PU_Weight->at(pu);
	
     //prompt muon 
	mu_promptPt_ = mu_pt_->at(FirstMuon);
	mu_promptEta_ = mu_eta_->at(FirstMuon);
	mu_promptPhi_ = mu_phi_->at(FirstMuon);
	mu_promptCharge_ = mu_charge_->at(FirstMuon);
	mu_promptE_ = mu_en_->at(FirstMuon);
	mu_promptEt_ = mu_et_->at(FirstMuon);
	mu_promptRhoIso_   = mu_rhoIso_->at(FirstMuon);
	mu_promptTrackiso_ = mu_trackiso_->at(FirstMuon);
	mu_promptPfSumChHadPt_  = mu_pfSumChargedHadronPt_->at(FirstMuon);   
	mu_promptPfSumNHadEt_   = mu_pfSumNeutralHadronEt_->at(FirstMuon);  
	mu_promptPFSumPhotonEt_ = mu_PFSumPhotonEt_->at(FirstMuon);  
	mu_promptPfSumPUPt_ = mu_pfSumPUPt_->at(FirstMuon);   
	mu_promptEmIso_     = mu_emIso_->at(FirstMuon);  
	mu_promptHadIso_    = mu_hadIso_->at(FirstMuon); 
	//mu_promptNormalizedChi2_ = mu_normalizedChi2_->at(FirstMuon);
	mu_promptDPToverPT_      = mu_dPToverPTTunePMuonBestTrack_->at(FirstMuon); 
	mu_promptAbsdxy_      = mu_absdxyTunePMuonBestTrack_->at(FirstMuon);   
	mu_promptAbsdxyError_ = mu_absdxyErrorTunePMuonBestTrack_->at(FirstMuon);
	mu_promptAbsdxySig_   = mu_absdxySigTunePMuonBestTrack_->at(FirstMuon); 
	mu_promptAbsdz_ = mu_absdzTunePMuonBestTrack_->at(FirstMuon);  
	mu_promptAbsdz_ = mu_absdzErrorTunePMuonBestTrack_->at(FirstMuon);
	mu_promptAbsdzSig_      = mu_absdzSigTunePMuonBestTrack_->at(FirstMuon);  
	mu_promptRecoDeltaBeta_ = deltaBeta_->at(FirstMuon);  
	mu_promptRecoiso_   = mu_recoiso_->at(FirstMuon);  
	
	//mu_promptDirection_ = mu_STATofDirection_->at(FirstMuon); 
	mu_promptNDof_      = mu_STATofNDof_->at(FirstMuon);  
	mu_promptTimeAtIpInOut_    = mu_STATofTimeAtIpInOut_->at(FirstMuon);
	mu_promptTimeAtIpInOutErr_ = mu_STATofTimeAtIpInOutErr_->at(FirstMuon);
	mu_promptTimeAtIpOutIn_   = mu_STATofTimeAtIpOutIn_->at(FirstMuon); 
	mu_promptTimeAtIpOutInErr_ = mu_STATofTimeAtIpOutInErr_->at(FirstMuon);
	
	//mu_promptMatchedStations_  = mu_numberOfMatchedStations_->at(FirstMuon); 
	//mu_promptValidPixelHits_   = mu_numberOfValidPixelHits_->at(FirstMuon);  
	//mu_promptTrackQuality_     = mu_TrackQuality_->at(FirstMuon);   
	//mu_promptInrTrackQuality_  = mu_InnerTrackQuality_->at(FirstMuon);
	//mu_promptValidMuonHits_    = mu_numberOfValidMuonHits_->at(FirstMuon); 
	//mu_promptGenMatch_         = mu_FirstGenMatch_->at(FirstMuon); 
	
	//non_prompt muon
	mu_secondPt_ = mu_pt_->at(SecondMuon);
	mu_secondEta_ = mu_eta_->at(SecondMuon);
	mu_secondPhi_ = mu_phi_->at(SecondMuon);
	mu_secondCharge_ = mu_charge_->at(SecondMuon);
	mu_secondE_ = mu_en_->at(SecondMuon);
	mu_secondEt_ = mu_et_->at(SecondMuon);
	mu_secondRhoIso_   = mu_rhoIso_->at(SecondMuon);
	mu_secondTrackiso_ = mu_trackiso_->at(SecondMuon);
	mu_secondPfSumChHadPt_  = mu_pfSumChargedHadronPt_->at(SecondMuon);   
	mu_secondPfSumNHadEt_   = mu_pfSumNeutralHadronEt_->at(SecondMuon);  
	mu_secondPFSumPhotonEt_ = mu_PFSumPhotonEt_->at(SecondMuon);  
	mu_secondPfSumPUPt_ = mu_pfSumPUPt_->at(SecondMuon);   
	mu_secondEmIso_     = mu_emIso_->at(SecondMuon);  
	mu_secondHadIso_    = mu_hadIso_->at(SecondMuon); 
	//mu_secondNormalizedChi2_ = mu_normalizedChi2_->at(SecondMuon);
	mu_secondDPToverPT_      = mu_dPToverPTTunePMuonBestTrack_->at(SecondMuon); 
	mu_secondAbsdxy_      = mu_absdxyTunePMuonBestTrack_->at(SecondMuon);   
	mu_secondAbsdxyError_ = mu_absdxyErrorTunePMuonBestTrack_->at(SecondMuon);
	mu_secondAbsdxySig_   = mu_absdxySigTunePMuonBestTrack_->at(SecondMuon); 
	mu_secondAbsdz_ = mu_absdzTunePMuonBestTrack_->at(SecondMuon);  
	mu_secondAbsdz_ = mu_absdzErrorTunePMuonBestTrack_->at(SecondMuon);
	mu_secondAbsdzSig_      = mu_absdzSigTunePMuonBestTrack_->at(SecondMuon);  
	mu_secondRecoDeltaBeta_ = deltaBeta_->at(SecondMuon);  
	mu_secondRecoiso_   = mu_recoiso_->at(SecondMuon);  
	
	//mu_secondDirection_ = mu_STATofDirection_->at(SecondMuon); 
	mu_secondNDof_      = mu_STATofNDof_->at(SecondMuon);  
	mu_secondTimeAtIpInOut_    = mu_STATofTimeAtIpInOut_->at(SecondMuon);
	mu_secondTimeAtIpInOutErr_ = mu_STATofTimeAtIpInOutErr_->at(SecondMuon);
	mu_secondTimeAtIpOutIn_    = mu_STATofTimeAtIpOutIn_->at(SecondMuon); 
	mu_secondTimeAtIpOutInErr_ = mu_STATofTimeAtIpOutInErr_->at(SecondMuon);
	
	//mu_secondMatchedStations_  = mu_numberOfMatchedStations_->at(SecondMuon); 
	//mu_secondValidPixelHits_   = mu_numberOfValidPixelHits_->at(SecondMuon);  
	//mu_secondTrackQuality_     = mu_TrackQuality_->at(SecondMuon);   
	//mu_secondInrTrackQuality_  = mu_InnerTrackQuality_->at(SecondMuon);
	//mu_secondValidMuonHits_    = mu_numberOfValidMuonHits_->at(SecondMuon); 
	//mu_secondGenMatch_         = mu_SecondGenMatch_->at(SecondMuon);
	
	
	float R = sqrt(((mu_secondEta_-mu_promptEta_)*(mu_secondEta_-mu_promptEta_))+((mu_secondPhi_-mu_promptPhi_)*(mu_secondPhi_-mu_promptPhi_)));
	float charged =  mu_pfSumChargedHadronPt_->at(SecondMuon);
	float neutral =  mu_pfSumNeutralHadronEt_->at(SecondMuon);
	float sumPhotonEt = mu_PFSumPhotonEt_->at(SecondMuon);
	float pileup =  mu_pfSumPUPt_->at(SecondMuon); 
	
	double deltaBetaR3 = (charged + std::max(0.0, neutral+sumPhotonEt-0.5*pileup))/mu_secondPt_;
	
	TLorentzVector Mu1;
	TLorentzVector Mu2; 
	Mu1.SetPtEtaPhiE(mu_promptPt_,mu_promptEta_,mu_promptPhi_,mu_promptE_);
	Mu2.SetPtEtaPhiE(mu_secondPt_,mu_secondEta_,mu_secondPhi_,mu_secondE_);
	float DiMuMass = (Mu1 + Mu2).M();
	
	mu_DeltaBetaR3_ = deltaBetaR3;
	mu_DiMuMass_ = DiMuMass;
	mu_Size_ = count;
	mu_DeltaR_ = R ;
    
	//secondary vertex info   
	//sv_TrackSize_new_ =  sv_mu_TrackSize_->at(SecondVertex);
	sv_LXYSig_new_ =  sv_mu_LXYSig_->at(SecondVertex);
	sv_LXYZSig_new_ =  sv_mu_LXYZSig_->at(SecondVertex);
	sv_LX_new_ =  sv_mu_LX_->at(SecondVertex);
	sv_LY_new_ =  sv_mu_LY_->at(SecondVertex);
	sv_LZ_new_ =  sv_mu_LZ_->at(SecondVertex);
	sv_LXY_new_ =  sv_mu_LXY_->at(SecondVertex);
	sv_LXYZ_new_ =  sv_mu_LXYZ_->at(SecondVertex);
	sv_mass_new_ =  sv_mu_mass_->at(SecondVertex);
	sv_eta_new_ =  sv_mu_eta_->at(SecondVertex);
	sv_phi_new_ =  sv_mu_phi_->at(SecondVertex);
	sv_pt_new_ =  sv_mu_pt_->at(SecondVertex);
	sv_p_new_ =  sv_mu_p_->at(SecondVertex);
	sv_en_new_ =  sv_mu_en_->at(SecondVertex);
	sv_Beta_new_ =  sv_mu_Beta_->at(SecondVertex);
	sv_Gamma_new_ =  sv_mu_Gamma_->at(SecondVertex);
	sv_CTau0_new_ =  sv_mu_CTau0_->at(SecondVertex);
	sv_NDof_new_ =  sv_mu_NDof_->at(SecondVertex);
	sv_Chi2_new_ =  sv_mu_Chi2_->at(SecondVertex);
	sv_Angle3D_new_ =  sv_mu_Angle3D_->at(SecondVertex);
	sv_Angle2D_new_ =  sv_mu_Angle2D_->at(SecondVertex);
	sv_tracks_Sumcharge_new_ =  sv_mu_tracks_Sumcharge_->at(SecondVertex);
	sv_tracks_Sumpt_new_ =  sv_mu_tracks_Sumpt_->at(SecondVertex);
	//sv_match_dxyz_new_ =  sv_mu_match_dxyz_->at(SecondVertex);
    
	//jet info
	unsigned  jets = -1;
	Float_t pt_jets = -1000;
	for(unsigned i=0; i<jet_pt->size(); i++){
	  if(jet_pt->at(i) > pt_jets) {
	    pt_jets = jet_pt->at(i);
	    jets = i;
	  }
	}
	
	jet_charge_ = jet_charge->at(jets);
	jet_et_ = jet_et->at(jets);
	jet_pt_ = jet_pt->at(jets);
	jet_eta_ = jet_eta->at(jets); 
	jet_phi_ = jet_phi->at(jets);
	jet_theta_ = jet_theta->at(jets);
	jet_en_ = jet_en->at(jets);
	jet_chEmEn_ = jet_chargedEmEnergy->at(jets); 
	jet_NEmEnFraction_ = jet_neutralEmEnergyFraction->at(jets);
	jet_chHadEn_ = jet_chargedHadronEnergy->at(jets);
	jet_NHadEnFraction_ = jet_neutralHadronEnergy->at(jets);
	jet_chMuEn_ = jet_chargedMuEnergy->at(jets);
	jet_chMuEnFraction_ = jet_chargedMuEnergyFraction->at(jets);
	jet_numberOfDaughters_ = jet_numberOfDaughters->at(jets);
	jet_muonEnergy_ = jet_muonEnergy->at(jets);
	jet_muonEnergyFraction_ = jet_muonEnergyFraction->at(jets);
	jet_muonMultiplicity_ = jet_muonMultiplicity->at(jets);
	jet_neutralEmEnergy_ = jet_neutralEmEnergy->at(jets);
	jet_neutralHadronEnergy_ = jet_neutralHadronEnergy->at(jets); 
	jet_NHadMultiplicity_ = jet_neutralHadronMultiplicity->at(jets);
	jet_NMultiplicity_ = jet_neutralMultiplicity->at(jets);
	
	
    //missing energy Info
	pfMet_et_    =  pfMet_et;
	pfMet_pt_    =  pfMet_pt;
	pfMet_phi_   =  pfMet_phi;
	pfMet_en_    =  pfMet_en;
	pfMet_sumEt_ =  pfMet_sumEt;
	caloMet_pt_  =  caloMet_pt;
	caloMet_phi_ =  caloMet_phi;
	/*     
	//bjet info 
	unsigned  bjet_bbb = -1;
	Float_t bjet_disc_bbb = 1000;
	for(unsigned i=0; i<jet_btag_pfDeepCSV_bbb_discriminator->size(); i++){
	if(jet_btag_pfDeepCSV_bbb_discriminator->at(i) < bjet_disc_bbb) {
	bjet_disc_bbb = jet_btag_pfDeepCSV_bbb_discriminator->at(i);
        bjet_bbb = i;
	//jet_btag_bbb_discriminator_ = jet_btag_pfDeepCSV_bbb_discriminator->at(i);
	}
    }

    unsigned  bjet_bb = -1;
    Float_t bjet_disc_bb = 1000;
    for(unsigned i=0; i<jet_btag_pfDeepCSV_bb_discriminator->size(); i++){
      if(jet_btag_pfDeepCSV_bb_discriminator->at(i) < bjet_disc_bb) {
        bjet_disc_bb = jet_btag_pfDeepCSV_bb_discriminator->at(i);
	bjet_bb = i;
	//jet_btag_bb_discriminator_ = jet_btag_pfDeepCSV_bb_discriminator->at(bjet_bb);
      }
    }
    
    unsigned  bjet_bc = -1;
    Float_t bjet_disc_bc = 1000;
    for(unsigned i=0; i<jet_btag_pfDeepCSV_bc_discriminator->size(); i++){
      if(jet_btag_pfDeepCSV_bc_discriminator->at(i) < bjet_disc_bc) {
        bjet_disc_bc = jet_btag_pfDeepCSV_bc_discriminator->at(i);
        bjet_bc = i;
	//jet_btag_bc_discriminator_ = jet_btag_pfDeepCSV_bc_discriminator->at(bjet_bc);
      }
    }


    jet_btag_pt_  = jet_btag_pt->at(bjet_bbb);
    jet_btag_eta_ = jet_btag_eta->at(bjet_bbb);
    jet_btag_phi_ = jet_btag_phi->at(bjet_bbb);
    jet_btag_flavor_ = jet_btag_flavor->at(bjet_bbb);
    jet_btag_bbb_discriminator_ = jet_btag_pfDeepCSV_bbb_discriminator->at(bjet_bbb);
    jet_btag_bb_discriminator_ = jet_btag_pfDeepCSV_bb_discriminator->at(bjet_bb);
    jet_btag_bc_discriminator_ = jet_btag_pfDeepCSV_bc_discriminator->at(bjet_bc);
    */
    
    /*
    for(unsigned i=0; i<jet_btag_pfDeepCSV_bc_discriminator->size(); i++){
      if(jet_btag_pfDeepCSV_bc_discriminator->at(i) > -1){  
      jet_btag_bc_discriminator_ = jet_btag_pfDeepCSV_bc_discriminator->at(i);
      jet_btag_pt_  = jet_btag_pt->at(i);
      jet_btag_eta_ = jet_btag_eta->at(i);
      jet_btag_phi_ = jet_btag_phi->at(i);
      }
    }

    for(unsigned i=0; i<jet_btag_pfDeepCSV_bbb_discriminator->size(); i++){
      if(jet_btag_pfDeepCSV_bbb_discriminator->at(i) > -1){
	jet_btag_bbb_discriminator_ = jet_btag_pfDeepCSV_bc_discriminator->at(i);
      }
    }

    for(unsigned i=0; i<jet_btag_pfDeepCSV_bb_discriminator->size(); i++){
      if(jet_btag_pfDeepCSV_bb_discriminator->at(i) > -1){
	jet_btag_bb_discriminator_ = jet_btag_pfDeepCSV_bc_discriminator->at(i);
      }
    }
    */
    
	//transverse Mass info


	TLorentzVector ivf;
	ivf.SetPtEtaPhiE(sv_pt_new_,sv_eta_new_,sv_phi_new_,sv_en_new_);
	TLorentzVector ivfPluslep1 = ivf + Mu1;


	float transverse_mass_ivf = sqrt(pow(ivf.Pt() + pfMet_pt, 2) - pow(ivf.Px() + pfMet_px, 2) - pow(ivf.Py() + pfMet_py, 2));//ivf transverse mass
	float transverse_mass_lep1 = sqrt(pow(Mu1.Pt() + pfMet_pt, 2) - pow(Mu1.Px() + pfMet_px, 2) - pow(Mu1.Py() + pfMet_py, 2));//prompt lepton transverse mass
	float transverse_mass_ivfPluslep1 = sqrt(pow(ivfPluslep1.Pt() + pfMet_pt, 2) - pow(ivfPluslep1.Px() + pfMet_px, 2) - pow(ivfPluslep1.Py() + pfMet_py, 2));//ivf + prompt lepton transverse mass

	  //mass correction
	TVector3 dir_sv;

	double ivf_mass = ivf.M();
	dir_sv.SetXYZ(sv_LX_new_, sv_LY_new_, sv_LZ_new_);	
	double vertexPt2 = dir_sv.Cross(ivf.Vect()).Mag2() / dir_sv.Mag2();
	double mass_corrected = std::sqrt(ivf_mass * ivf_mass + vertexPt2) + std::sqrt(vertexPt2);


	trMass_ivf_ = transverse_mass_ivf;
	trMass_lep1_ = transverse_mass_lep1;
	trMass_ivfPluslep1_ = transverse_mass_ivfPluslep1;

	//sv mass correction
	sv_mass_correction_ = mass_corrected;

	//rpc info
	mu_rpc_nDof_newtree_ = mu_rpc_nDof_;
	mu_rpc_time_newtree_ = mu_rpc_time_;
	mu_rpc_timeErr_newtree_ = mu_rpc_timeErr_;
	mu_cmb_nDof_newtree_ = mu_cmb_nDof_;
	mu_cmb_time_newtree_ = mu_cmb_time_;
	mu_cmb_timeErr_newtree_ = mu_cmb_timeErr_;

	newtree->Fill();
      }
      
    }
     
  newtree->Print();
  newtree->AutoSave();

  delete newfile;
  }//chiusura del while eventuale
  return 0;
}


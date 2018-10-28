//**********************************************************************************************************************************
// Remove some branches + selects the events + add variables -- for muonic channel
//***************************************** To Compile******************************************************************************
// g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit CloneTree.C -o CloneTree.exe
//**********************************************************************************************************************************

#ifndef __CINT__
#include "RooGlobalFunc.h"
//------------------------------------------------
     
#endif
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooStats/SPlot.h"
#include <vector>
#include <string>
#include <iostream>
#include "RooRandom.h"
#include "RooMinuit.h"
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
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooArgusBG.h"
#include "TString.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMappedCategory.h"
#include "RooCmdArg.h"
#include "RooChebychev.h"
#include "RooUnblindUniform.h"
#include "RooUnblindPrecision.h"
#include "RooExponential.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooSimWSTool.h"
#include "RooWorkspace.h"
#include <TLatex.h>
#include "RooFit.h"
#include "RooConstVar.h"
#include "RooSimPdfBuilder.h"
#include "RooStringVar.h"
#include "TText.h"
#include "TPaveText.h"
#include "RooHist.h"
#include "TLorentzVector.h"


using namespace std;
using namespace RooFit;
using namespace RooStats;


int main(){
  

  //Get old file, old tree and set top branch address
  //TFile *oldfile = new TFile("/eos/cms/store/group/phys_exotica/HNL/Background/crab_Analysis_WZToLLLNu/Background_Analysis.root");
  //TTree *oldtree = (TTree*)oldfile->Get("HeavyNeutralLepton/tree_");


  //////selection cuts

  Float_t isoCut = 0.15;
  bool isMC = true;
  //if I want to use a TChain.....
  //cout<< "starting..."<<endl;
  TChain * oldtree = new TChain("HeavyNeutralLepton/tree_","");

  oldtree->Add("/eos/cms/store/group/phys_exotica/HNL/Background/crab_Analysis_WZToLLLNu_ext1/Background_Analysis.root/HeavyNeutralLepton/tree_");

  TFile *newfile = new TFile("skimmedSignale.root","recreate");
  //Create a new file + a clone of old tree in new file 
  //TTree *newtree = oldtree->CloneTree(0); 
  TTree *newtree  = new TTree("newtree","Analysis Tree");
  //cout<<"cloning done"<<endl;

  // Long64_t nentries = oldtree->GetEntries();
  Long64_t nentries = oldtree->GetEntriesFast();
  cout  <<nentries<<endl;

//======================= Old Tree Variables ==========================================// 
  // These are the variables I cut on 
  Bool_t passIsoMu24All;
  Bool_t passIsoMu27All;

  oldtree->SetBranchAddress("passIsoMu24All",&passIsoMu24All);
  oldtree->SetBranchAddress("passIsoMu27All",&passIsoMu27All);

  vector<Float_t>   *PU_Weight = 0;
  oldtree->SetBranchAddress("PU_Weight",&PU_Weight);

  vector<Float_t>   *mu_isTightMuon = 0;
  vector<Float_t>   *mu_isLoose = 0;
  vector<Float_t>   *mu_pt = 0;
  vector<Float_t>   *mu_eta =0;
  vector<Float_t>   *mu_phi=0; 
  vector<Float_t>   *mu_charge=0; 
  vector<Float_t>   *mu_en=0;
  vector<Float_t>   *mu_et=0;
  vector<Float_t>   *deltaBeta=0;
  vector<Float_t>   *mu_rhoIso = 0 ;
  vector<Float_t>   *mu_trackiso = 0 ;
  vector<Float_t>   *mu_pfSumChargedHadronPt = 0 ;
  vector<Float_t>   *mu_pfSumNeutralHadronEt = 0 ;
  vector<Float_t>   *mu_PFSumPhotonEt = 0 ;
  vector<Float_t>   *mu_pfSumPUPt = 0 ;
  vector<Float_t>   *mu_emIso = 0 ;
  vector<Float_t>   *mu_hadIso = 0 ;
  vector<Float_t>   *mu_normalizedChi2 = 0 ;
  vector<Float_t>   *mu_dPToverPTTunePMuonBestTrack = 0 ;
  vector<Float_t>   *mu_absdxyTunePMuonBestTrack = 0 ;
  vector<Float_t>   *mu_absdxyErrorTunePMuonBestTrack = 0 ;
  vector<Float_t>   *mu_absdxySigTunePMuonBestTrack = 0 ;
  vector<Float_t>   *mu_absdzTunePMuonBestTrack = 0 ;
  vector<Float_t>   *mu_absdzErrorTunePMuonBestTrack = 0 ;
  vector<Float_t>   *mu_absdzSigTunePMuonBestTrack = 0 ;
  vector<Float_t>   *mu_recoiso = 0 ;
  vector<Float_t>   *mu_STATofDirection = 0 ;
  vector<Float_t>   *mu_STATofNDof = 0 ;
  vector<Float_t>   *mu_STATofTimeAtIpInOut = 0 ;
  vector<Float_t>   *mu_STATofTimeAtIpInOutErr = 0 ;
  vector<Float_t>   *mu_STATofTimeAtIpOutIn = 0 ;
  vector<Float_t>   *mu_STATofTimeAtIpOutInErr = 0 ;
  vector<Int_t>     *mu_numberOfValidMuonHits = 0 ;
  vector<Int_t>     *mu_numberOfMatchedStations = 0 ;
  vector<Int_t>     *mu_numberOfValidPixelHits = 0 ;
  vector<Int_t>     *mu_TrackQuality = 0 ;
  vector<Int_t>     *mu_InnerTrackQuality = 0 ;
  vector<Int_t>     *mu_FirstGenMatch = 0;
  vector<Int_t>     *mu_SecondGenMatch = 0;

  vector<Int_t>   *sv_mu_TrackSize = 0 ;
  vector<Float_t> *sv_mu_LXYSig = 0 ;
  vector<Float_t> *sv_mu_LXYZSig = 0 ;
  vector<Float_t> *sv_mu_LXY = 0 ;
  vector<Float_t> *sv_mu_LXYZ = 0 ;
  vector<Float_t> *sv_mu_mass = 0 ;
  vector<Float_t> *sv_mu_eta = 0 ;
  vector<Float_t> *sv_mu_phi = 0 ;
  vector<Float_t> *sv_mu_pt = 0 ;
  vector<Float_t> *sv_mu_p = 0 ;
  vector<Float_t> *sv_mu_Beta = 0 ;
  vector<Float_t> *sv_mu_Gamma = 0 ;
  vector<Float_t> *sv_mu_CTau0 = 0 ;
  vector<Float_t> *sv_mu_NDof = 0 ;
  vector<Float_t> *sv_mu_Chi2 = 0 ;
  vector<Float_t> *sv_mu_Angle3D = 0 ;
  vector<Float_t> *sv_mu_Angle2D = 0 ;
  vector<Int_t>   *sv_mu_tracks_Sumcharge = 0 ;
  vector<Float_t> *sv_mu_tracks_Sumpt = 0 ;
  vector<Float_t> *sv_mu_match = 0 ;

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
  Float_t   pfMet_phi;
  Float_t   pfMet_en;
  Float_t   pfMet_sumEt;
  Float_t   caloMet_pt;
  Float_t   caloMet_phi;

  vector<Float_t>   *jet_btag_pt = 0;
  vector<Float_t>   *jet_btag_eta = 0;
  vector<Float_t>   *jet_btag_phi = 0;
  vector<Int_t>     *jet_btag_flavor = 0;
  vector<Float_t>   *jet_btag_pfCSVv2IVF_discriminator = 0;

 
  oldtree->SetBranchAddress("mu_isTightMuon",&mu_isTightMuon);
  oldtree->SetBranchAddress("mu_isLoose",&mu_isLoose);
  oldtree->SetBranchAddress("mu_pt",&mu_pt);
  oldtree->SetBranchAddress("mu_eta",&mu_eta);
  oldtree->SetBranchAddress("mu_phi",&mu_phi);
  oldtree->SetBranchAddress("mu_charge",&mu_charge);
  oldtree->SetBranchAddress("mu_en",&mu_en);
  oldtree->SetBranchAddress("mu_et",&mu_et);
  oldtree->SetBranchAddress("mu_recoDeltaBeta",&deltaBeta);
  oldtree->SetBranchAddress("mu_rhoIso", &mu_rhoIso);
  oldtree->SetBranchAddress("mu_trackiso", &mu_trackiso);
  oldtree->SetBranchAddress("mu_pfSumChargedHadronPt", &mu_pfSumChargedHadronPt);
  oldtree->SetBranchAddress("mu_pfSumNeutralHadronEt", &mu_pfSumNeutralHadronEt);
  oldtree->SetBranchAddress("mu_PFSumPhotonEt", &mu_PFSumPhotonEt);
  oldtree->SetBranchAddress("mu_pfSumPUPt", &mu_pfSumPUPt);
  oldtree->SetBranchAddress("mu_emIso", &mu_emIso);
  oldtree->SetBranchAddress("mu_hadIso", &mu_hadIso);
  oldtree->SetBranchAddress("mu_normalizedChi2", &mu_normalizedChi2);
  oldtree->SetBranchAddress("mu_dPToverPTTunePMuonBestTrack", &mu_dPToverPTTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_absdxyTunePMuonBestTrack", &mu_absdxyTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_absdxyErrorTunePMuonBestTrack", &mu_absdxyErrorTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_absdxySigTunePMuonBestTrack", &mu_absdxySigTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_absdzTunePMuonBestTrack", &mu_absdzTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_absdzErrorTunePMuonBestTrack", &mu_absdzErrorTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_absdzSigTunePMuonBestTrack", &mu_absdzSigTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_recoiso", &mu_recoiso);
  oldtree->SetBranchAddress("mu_STATofDirection", &mu_STATofDirection);
  oldtree->SetBranchAddress("mu_STATofNDof", &mu_STATofNDof);
  oldtree->SetBranchAddress("mu_STATofTimeAtIpInOut", &mu_STATofTimeAtIpInOut);
  oldtree->SetBranchAddress("mu_STATofTimeAtIpInOutErr", &mu_STATofTimeAtIpInOutErr);
  oldtree->SetBranchAddress("mu_STATofTimeAtIpOutIn", &mu_STATofTimeAtIpOutIn);
  oldtree->SetBranchAddress("mu_STATofTimeAtIpOutInErr", &mu_STATofTimeAtIpOutInErr);
  oldtree->SetBranchAddress("mu_numberOfValidMuonHits", &mu_numberOfValidMuonHits);
  oldtree->SetBranchAddress("mu_numberOfMatchedStations", &mu_numberOfMatchedStations);
  oldtree->SetBranchAddress("mu_numberOfValidPixelHits", &mu_numberOfValidPixelHits);
  oldtree->SetBranchAddress("mu_TrackQuality", &mu_TrackQuality);
  oldtree->SetBranchAddress("mu_InnerTrackQuality", &mu_InnerTrackQuality);
  oldtree->SetBranchAddress("mu_FirstGenMatch", &mu_FirstGenMatch);
  oldtree->SetBranchAddress("mu_SecondGenMatch", &mu_SecondGenMatch);

  oldtree->SetBranchAddress("sv_mu_TrackSize", &sv_mu_TrackSize);
  oldtree->SetBranchAddress("sv_mu_LXYSig", &sv_mu_LXYSig);
  oldtree->SetBranchAddress("sv_mu_LXYZSig", &sv_mu_LXYZSig);
  oldtree->SetBranchAddress("sv_mu_LXY", &sv_mu_LXY);
  oldtree->SetBranchAddress("sv_mu_LXYZ", &sv_mu_LXYZ);
  oldtree->SetBranchAddress("sv_mu_mass", &sv_mu_mass);
  oldtree->SetBranchAddress("sv_mu_eta", &sv_mu_eta);   
  oldtree->SetBranchAddress("sv_mu_phi", &sv_mu_phi);  
  oldtree->SetBranchAddress("sv_mu_pt", &sv_mu_pt);  
  oldtree->SetBranchAddress("sv_mu_p", &sv_mu_p);  
  oldtree->SetBranchAddress("sv_mu_Beta", &sv_mu_Beta);
  oldtree->SetBranchAddress("sv_mu_Gamma", &sv_mu_Gamma); 
  oldtree->SetBranchAddress("sv_mu_CTau0", &sv_mu_CTau0);  
  oldtree->SetBranchAddress("sv_mu_NDof", &sv_mu_NDof);   
  oldtree->SetBranchAddress("sv_mu_Chi2", &sv_mu_Chi2);   
  oldtree->SetBranchAddress("sv_mu_Angle3D", &sv_mu_Angle3D);
  oldtree->SetBranchAddress("sv_mu_Angle2D", &sv_mu_Angle2D);  
  oldtree->SetBranchAddress("sv_mu_tracks_Sumcharge", &sv_mu_tracks_Sumcharge); 
  oldtree->SetBranchAddress("sv_mu_tracks_Sumpt", &sv_mu_tracks_Sumpt);  
  oldtree->SetBranchAddress("sv_mu_match", &sv_mu_match);

  oldtree->SetBranchAddress("jet_btag_pt", &jet_btag_pt);
  oldtree->SetBranchAddress("jet_btag_eta", &jet_btag_eta);
  oldtree->SetBranchAddress("jet_btag_phi", &jet_btag_phi);
  oldtree->SetBranchAddress("jet_btag_flavor", &jet_btag_flavor);
  oldtree->SetBranchAddress("jet_btag_pfCSVv2IVF_discriminator", &jet_btag_pfCSVv2IVF_discriminator);

  oldtree->SetBranchAddress("pfMet_et", &pfMet_et);
  oldtree->SetBranchAddress("pfMet_pt", &pfMet_pt);
  oldtree->SetBranchAddress("pfMet_phi", &pfMet_phi);
  oldtree->SetBranchAddress("pfMet_en", &pfMet_en);
  oldtree->SetBranchAddress("pfMet_sumEt", &pfMet_sumEt);
  oldtree->SetBranchAddress("caloMet_pt", &caloMet_pt);
  oldtree->SetBranchAddress("caloMet_phi", &caloMet_phi);

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

  //Throw away some branches -- thos are all the  one I don't want to keep in the ntuples
  //IMPORTANT: we cannot throw away the branches with the variable we are using in the loop!
  oldtree->SetBranchStatus("ele_*",    0);
  //oldtree->SetBranchStatus("mu_ST*", 0);

//================================ PU Weight ============================================// 
Float_t pu_weight;

TBranch* branch_pu_weight = newtree->Branch("pu_weight",&pu_weight,"pu_weight/F");
//======================= First Muon Variables ==========================================//   
Float_t mu_promptPt,mu_promptEta, mu_promptPhi, mu_promptCharge, mu_promptEt, mu_promptE;

Float_t mu_promptRhoIso,        mu_promptTrackiso,       mu_promptPfSumChHadPt,
	mu_promptPFSumPhotonEt, mu_promptPfSumPUPt,      mu_promptEmIso,
	mu_promptHadIso,        mu_promptNormalizedChi2, mu_promptDPToverPT,
	mu_promptAbsdxy,        mu_promptAbsdxyError,    mu_promptAbsdxySig,
	mu_promptAbsdz,         mu_promptAbsdzError,     mu_promptAbsdzSig,
	mu_promptRecoDeltaBeta, mu_promptRecoiso,        mu_promptDirection,
	mu_promptNDof,          mu_promptTimeAtIpInOut,  mu_promptTimeAtIpInOutErr,
	mu_promptTimeAtIpOutIn, mu_promptTimeAtIpOutInErr,  mu_promptPfSumNHadEt;

 Int_t   mu_promptValidMuonHits, mu_promptMatchedStations, mu_promptValidPixelHits, 
         mu_promptTrackQuality,  mu_promptInrTrackQuality, mu_promptGenMatch;


TBranch* branch_mu_promptPt = newtree->Branch("mu_promptPt",&mu_promptPt,"mu_promptPt/F");
TBranch* branch_mu_promptEta = newtree->Branch("mu_promptEta",&mu_promptEta,"mu_promptEta/F");
TBranch* branch_mu_promptPhi = newtree->Branch("mu_promptPhi",&mu_promptPhi,"mu_promptPhi/F");
TBranch* branch_mu_promptCharge = newtree->Branch("mu_promptCharge",&mu_promptCharge,"mu_promptCharge/F");
TBranch* branch_mu_promptE = newtree->Branch("mu_promptE",&mu_promptE,"mu_promptE/F");
TBranch* branch_mu_promptEt = newtree->Branch("mu_promptEt",&mu_promptEt,"mu_promptEt/F");
TBranch* branch_mu_promptRhoIso = newtree->Branch("mu_promptRhoIso", &mu_promptRhoIso, "mu_promptRhoIso/F");
TBranch* branch_mu_promptTrackiso = newtree->Branch("mu_promptTrackiso", &mu_promptTrackiso, "mu_promptTrackiso/F");
TBranch* branch_mu_promptPfSumChHadPt = newtree->Branch("mu_promptPfSumChHadPt", &mu_promptPfSumChHadPt, "mu_promptPfSumChHadPt/F");
TBranch* branch_mu_promptPfSumNHadEt = newtree->Branch("mu_promptPfSumNHadEt", &mu_promptPfSumNHadEt, "mu_promptPfSumNHadEt/F");
TBranch* branch_mu_promptPFSumPhotonEt = newtree->Branch("mu_promptPFSumPhotonEt", &mu_promptPFSumPhotonEt, "mu_promptPFSumPhotonEt/F");
TBranch* branch_mu_promptPfSumPUPt = newtree->Branch("mu_promptPfSumPUPt", &mu_promptPfSumPUPt, "mu_promptPfSumPUPt/F");
TBranch* branch_mu_promptEmIso = newtree->Branch("mu_promptEmIso", &mu_promptEmIso, "mu_promptEmIso/F");
TBranch* branch_mu_promptHadIso = newtree->Branch("mu_promptHadIso", &mu_promptHadIso, "mu_promptHadIso/F");
TBranch* branch_mu_promptNormalizedChi2 = newtree->Branch("mu_promptNormalizedChi2", &mu_promptNormalizedChi2, "mu_promptNormalizedChi2/F");
TBranch* branch_mu_promptDPToverPT = newtree->Branch("mu_promptDPToverPT", &mu_promptDPToverPT, "mu_promptDPToverPT/F");
TBranch* branch_mu_promptAbsdxy = newtree->Branch("mu_promptAbsdxy", &mu_promptAbsdxy, "mu_promptAbsdxy/F");
TBranch* branch_mu_promptAbsdxyError = newtree->Branch("mu_promptAbsdxyError", &mu_promptAbsdxyError, "mu_promptAbsdxyError/F");
TBranch* branch_mu_promptAbsdxySig = newtree->Branch("mu_promptAbsdxySig", &mu_promptAbsdxySig, "mu_promptAbsdxySig/F");
TBranch* branch_mu_promptAbsdz = newtree->Branch("mu_promptAbsdz", &mu_promptAbsdz, "mu_promptAbsdz/F");
TBranch* branch_mu_promptAbsdzError = newtree->Branch("mu_promptAbsdzError", &mu_promptAbsdzError, "mu_promptAbsdzError/F");
TBranch* branch_mu_promptAbsdzSig = newtree->Branch("mu_promptAbsdzSig", &mu_promptAbsdz, "mu_promptAbsdzSig/F");
TBranch* branch_mu_promptRecoDeltaBeta = newtree->Branch("mu_promptRecoDeltaBeta", &mu_promptRecoDeltaBeta, "mu_promptRecoDeltaBeta/F");
TBranch* branch_mu_promptRecoiso = newtree->Branch("mu_promptRecoiso", &mu_promptRecoiso, "mu_promptRecoiso/F");
TBranch* branch_mu_promptDirection = newtree->Branch("mu_promptDirection", &mu_promptDirection, "mu_promptDirection/F");
TBranch* branch_mu_promptNDof = newtree->Branch("mu_promptNDof", &mu_promptNDof, "mu_promptNDof/F");
TBranch* branch_mu_promptTimeAtIpInOut = newtree->Branch("mu_promptTimeAtIpInOut", &mu_promptTimeAtIpInOut, "mu_promptTimeAtIpInOut/F");
TBranch* branch_mu_promptTimeAtIpInOutErr = newtree->Branch("mu_promptTimeAtIpInOutErr", &mu_promptTimeAtIpInOutErr, "mu_promptTimeAtIpInOutErr/F");
TBranch* branch_mu_promptTimeAtIpOutIn = newtree->Branch("mu_promptTimeAtIpOutIn", &mu_promptTimeAtIpOutIn, "mu_promptTimeAtIpOutIn/F");
TBranch* branch_mu_promptTimeAtIpOutInErr = newtree->Branch("mu_promptTimeAtIpOutInErr", &mu_promptTimeAtIpOutInErr, "mu_promptTimeAtIpOutInErr/F");
TBranch* branch_mu_promptMatchedStations = newtree->Branch("mu_promptMatchedStations", &mu_promptMatchedStations, "mu_promptMatchedStations/I");
TBranch* branch_mu_promptValidPixelHits = newtree->Branch("mu_promptValidPixelHits", &mu_promptValidPixelHits, "mu_promptValidPixelHits/I");
TBranch* branch_mu_promptTrackQuality = newtree->Branch("mu_promptTrackQuality", &mu_promptTrackQuality, "mu_promptTrackQuality/I");
TBranch* branch_mu_promptInrTrackQuality = newtree->Branch("mu_promptInrTrackQuality", &mu_promptInrTrackQuality, "mu_promptInrTrackQuality/I");
TBranch* branch_mu_promptValidMuonHits = newtree->Branch("mu_promptValidMuonHits", &mu_promptValidMuonHits, "mu_promptValidMuonHits/I");
TBranch* branch_mu_promptGenMatch = newtree->Branch("mu_promptGenMatch", &mu_promptGenMatch, "mu_promptGenMatch/I");


//======================= Second Muon Variables ==========================================//
Float_t mu_secondPt,mu_secondEta, mu_secondPhi, mu_secondCharge, mu_secondEt, mu_secondE;

Float_t mu_secondRhoIso,        mu_secondTrackiso,  mu_secondPfSumChHadPt,
        mu_secondPFSumPhotonEt, mu_secondPfSumPUPt,      mu_secondEmIso,        
        mu_secondHadIso,        mu_secondNormalizedChi2, mu_secondDPToverPT,   
        mu_secondAbsdxy,        mu_secondAbsdxyError,    mu_secondAbsdxySig, 
        mu_secondAbsdz,         mu_secondAbsdzError,     mu_secondAbsdzSig,
        mu_secondRecoDeltaBeta, mu_secondRecoiso,        mu_secondDirection,
        mu_secondNDof,          mu_secondTimeAtIpInOut,  mu_secondTimeAtIpInOutErr, 
        mu_secondTimeAtIpOutIn, mu_secondTimeAtIpOutInErr,  mu_secondPfSumNHadEt;

 Float_t mu_DeltaBetaR3, mu_DiMuMass, mu_Size, mu_DeltaR;

 Int_t   mu_secondValidMuonHits, mu_secondMatchedStations, mu_secondValidPixelHits, 
         mu_secondTrackQuality,  mu_secondInrTrackQuality, mu_secondGenMatch;

TBranch* branch_mu_secondPt  = newtree->Branch("mu_secondPt",&mu_secondPt,"mu_secondPt/F");
TBranch* branch_mu_secondEta = newtree->Branch("mu_secondEta",&mu_secondEta,"mu_secondEta/F");
TBranch* branch_mu_secondPhi = newtree->Branch("mu_secondPhi",&mu_secondPhi,"mu_secondPhi/F");
TBranch* branch_mu_secondCharge = newtree->Branch("mu_secondCharge",&mu_secondCharge,"mu_secondCharge/F");
TBranch* branch_mu_secondE  = newtree->Branch("mu_secondE",&mu_secondE,"mu_secondE/F");
TBranch* branch_mu_secondEt = newtree->Branch("mu_secondEt",&mu_secondEt,"mu_secondEt/F");
TBranch* branch_mu_secondRhoIso   = newtree->Branch("mu_secondRhoIso", &mu_secondRhoIso, "mu_secondRhoIso/F");
TBranch* branch_mu_secondTrackiso = newtree->Branch("mu_secondTrackiso", &mu_secondTrackiso, "mu_secondTrackiso/F");
TBranch* branch_mu_secondPfSumChHadPt  = newtree->Branch("mu_secondPfSumChHadPt", &mu_secondPfSumChHadPt, "mu_secondPfSumChHadPt/F");
TBranch* branch_mu_secondPfSumNHadEt   = newtree->Branch("mu_secondPfSumNHadEt", &mu_secondPfSumNHadEt, "mu_secondPfSumNHadEt/F");
TBranch* branch_mu_secondPFSumPhotonEt = newtree->Branch("mu_secondPFSumPhotonEt", &mu_secondPFSumPhotonEt, "mu_secondPFSumPhotonEt/F");
TBranch* branch_mu_secondPfSumPUPt = newtree->Branch("mu_secondPfSumPUPt", &mu_secondPfSumPUPt, "mu_secondPfSumPUPt/F");
TBranch* branch_mu_secondEmIso  = newtree->Branch("mu_secondEmIso", &mu_secondEmIso, "mu_secondEmIso/F");
TBranch* branch_mu_secondHadIso = newtree->Branch("mu_secondHadIso", &mu_secondHadIso, "mu_secondHadIso/F");
TBranch* branch_mu_secondNormalizedChi2 = newtree->Branch("mu_secondNormalizedChi2", &mu_secondNormalizedChi2, "mu_secondNormalizedChi2/F");
TBranch* branch_mu_secondDPToverPT = newtree->Branch("mu_secondDPToverPT", &mu_secondDPToverPT, "mu_secondDPToverPT/F");
TBranch* branch_mu_secondAbsdxy = newtree->Branch("mu_secondAbsdxy", &mu_secondAbsdxy, "mu_secondAbsdxy/F");
TBranch* branch_mu_secondAbsdxyError = newtree->Branch("mu_secondAbsdxyError", &mu_secondAbsdxyError, "mu_secondAbsdxyError/F");
TBranch* branch_mu_secondAbsdxySig   = newtree->Branch("mu_secondAbsdxySig", &mu_secondAbsdxySig, "mu_secondAbsdxySig/F");
TBranch* branch_mu_secondAbsdz       = newtree->Branch("mu_secondAbsdz", &mu_secondAbsdz, "mu_secondAbsdz/F");
TBranch* branch_mu_secondAbsdzError  = newtree->Branch("mu_secondAbsdzError", &mu_secondAbsdzError, "mu_secondAbsdzError/F");
TBranch* branch_mu_secondAbsdzSig = newtree->Branch("mu_secondAbsdzSig", &mu_secondAbsdz, "mu_secondAbsdzSig/F");
TBranch* branch_mu_secondRecoDeltaBeta = newtree->Branch("mu_secondRecoDeltaBeta", &mu_secondRecoDeltaBeta, "mu_secondRecoDeltaBeta/F");
TBranch* branch_mu_secondRecoiso   = newtree->Branch("mu_secondRecoiso", &mu_secondRecoiso, "mu_secondRecoiso/F");
TBranch* branch_mu_secondDirection = newtree->Branch("mu_secondDirection", &mu_secondDirection, "mu_secondDirection/F");
TBranch* branch_mu_secondNDof = newtree->Branch("mu_secondNDof", &mu_secondNDof, "mu_secondNDof/F");
TBranch* branch_mu_secondTimeAtIpInOut    = newtree->Branch("mu_secondTimeAtIpInOut", &mu_secondTimeAtIpInOut, "mu_secondTimeAtIpInOut/F");
TBranch* branch_mu_secondTimeAtIpInOutErr = newtree->Branch("mu_secondTimeAtIpInOutErr", &mu_secondTimeAtIpInOutErr, "mu_secondTimeAtIpInOutErr/F");
TBranch* branch_mu_secondTimeAtIpOutIn = newtree->Branch("mu_secondTimeAtIpOutIn", &mu_secondTimeAtIpOutIn, "mu_secondTimeAtIpOutIn/F");
TBranch* branch_mu_secondTimeAtIpOutInErr = newtree->Branch("mu_secondTimeAtIpOutInErr", &mu_secondTimeAtIpOutInErr, "mu_secondTimeAtIpOutInErr/F");
TBranch* branch_mu_secondMatchedStations  = newtree->Branch("mu_secondMatchedStations", &mu_secondMatchedStations, "mu_secondMatchedStations/I");
TBranch* branch_mu_secondValidPixelHits   = newtree->Branch("mu_secondValidPixelHits", &mu_secondValidPixelHits, "mu_secondValidPixelHits/I");
TBranch* branch_mu_secondTrackQuality    = newtree->Branch("mu_secondTrackQuality", &mu_secondTrackQuality, "mu_secondTrackQuality/I");
TBranch* branch_mu_secondInrTrackQuality = newtree->Branch("mu_secondInrTrackQuality", &mu_secondInrTrackQuality, "mu_secondInrTrackQuality/I");
TBranch* branch_mu_secondValidMuonHits   = newtree->Branch("mu_secondValidMuonHits", &mu_secondValidMuonHits, "mu_secondValidMuonHits/I");
TBranch* branch_mu_secondGenMatch        = newtree->Branch("mu_secondGenMatch", &mu_secondGenMatch, "mu_secondGenMatch/I");

TBranch* branch_mu_DeltaBetaR3 = newtree->Branch("mu_DeltaBetaR3", &mu_DeltaBetaR3, "mu_DeltaBetaR3/F");
TBranch* branch_mu_DiMuMass    = newtree->Branch("mu_DiMuMass", &mu_DiMuMass, "mu_DiMuMass/F");
TBranch* branch_mu_DeltaR      = newtree->Branch("mu_DeltaR", &mu_DeltaR, "mu_DeltaR/F");
TBranch* branch_mu_Size        = newtree->Branch("mu_Size", &mu_Size, "mu_Size/F");

//======================= Second Vertex Variables ==========================================//
 Float_t  sv_LXYSig,  sv_LXYZSig, sv_LXY,  sv_LXYZ, sv_mass, 
          sv_eta,     sv_phi,     sv_pt,   sv_p,    sv_Beta, 
          sv_Gamma,   sv_CTau0,   sv_NDof, sv_Chi2, sv_Angle3D, 
          sv_Angle2D, sv_tracks_Sumpt,     sv_match;

 Int_t sv_TrackSize , sv_tracks_Sumcharge ;

 TBranch* branch_sv_TrackSize  = newtree->Branch("sv_TrackSize", &sv_TrackSize, "sv_TrackSize/I");
 TBranch* branch_sv_LXYSig     = newtree->Branch("sv_LXYSig", &sv_LXYSig, "sv_LXYSig/F");
 TBranch* branch_sv_LXYZSig    = newtree->Branch("sv_LXYZSig", &sv_LXYZSig, "sv_LXYZSig/F");
 TBranch* branch_sv_LXY        = newtree->Branch("sv_LXY", &sv_LXY, "sv_LXY/F");
 TBranch* branch_sv_LXYZ       = newtree->Branch("sv_LXYZ", &sv_LXYZ, "sv_LXYZ/F");
 TBranch* branch_sv_mass       = newtree->Branch("sv_mass", &sv_mass, "sv_mass/F");
 TBranch* branch_sv_eta        = newtree->Branch("sv_eta", &sv_eta, "sv_eta/F");
 TBranch* branch_sv_phi        = newtree->Branch("sv_phi", &sv_phi, "sv_phi/F");
 TBranch* branch_sv_pt         = newtree->Branch("sv_pt", &sv_pt, "sv_pt/F");
 TBranch* branch_sv_p          = newtree->Branch("sv_p", &sv_p, "sv_p/F");
 TBranch* branch_sv_Beta       = newtree->Branch("sv_Beta", &sv_Beta, "sv_Beta/F");
 TBranch* branch_sv_Gamma      = newtree->Branch("sv_Gamma", &sv_Gamma, "sv_Gamma/F");
 TBranch* branch_sv_CTau0      = newtree->Branch("sv_CTau0", &sv_CTau0, "sv_CTau0/F");
 TBranch* branch_sv_NDof       = newtree->Branch("sv_NDof", &sv_NDof, "sv_NDof/F");
 TBranch* branch_sv_Chi2       = newtree->Branch("sv_Chi2", &sv_Chi2, "sv_Chi2/F");
 TBranch* branch_sv_Angle3D    = newtree->Branch("sv_Angle3D", &sv_Angle3D, "sv_Angle3D/F");
 TBranch* branch_sv_Angle2D    = newtree->Branch("sv_Angle2D", &sv_Angle2D, "sv_Angle2D/F");
 TBranch* branch_sv_tracks_Sumcharge = newtree->Branch("sv_tracks_Sumcharge", &sv_tracks_Sumcharge, "sv_tracks_Sumcharge/I");
 TBranch* branch_sv_tracks_Sumpt     = newtree->Branch("sv_tracks_Sumpt", &sv_tracks_Sumpt, "sv_tracks_Sumpt/F");
 TBranch* branch_sv_match            = newtree->Branch("sv_match", &sv_match, "sv_match/F");
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
 Float_t jet_btag_pt_, jet_btag_eta_,  jet_btag_phi_, jet_btag_discriminator_;  
 Int_t   jet_btag_flavor_;

 TBranch* branch_jet_btag_pt_  = newtree->Branch("jet_btag_pt_", &jet_btag_pt_, "jet_btag_pt_/F");
 TBranch* branch_jet_btag_eta_  = newtree->Branch("jet_btag_eta_", &jet_btag_eta_, "jet_btag_eta_/F");
 TBranch* branch_jet_btag_phi_  = newtree->Branch("jet_btag_phi_", &jet_btag_phi_, "jet_btag_phi_/F");
 TBranch* branch_jet_btag_flavor_  = newtree->Branch("jet_btag_flavor_", &jet_btag_flavor_, "jet_btag_flavor_/F");
 TBranch* branch_jet_btag_discriminator_  = newtree->Branch("jet_btag_discriminator_", &jet_btag_discriminator_, "jet_btag_discriminator_/F");

//======================= Start the running over input branches ==========================================//
  for (int i=0;i<oldtree->GetEntriesFast(); i++) {
    if (i%10000==0) cout<<i<<endl;
    oldtree->GetEntry(i);

    if (passIsoMu24All==0 && passIsoMu27All == 0) continue;  // cut on the trigger!
    
    unsigned pu = -1;
    for(unsigned i=0; i<PU_Weight->size(); ++i){
      pu = i;
    }

    Float_t   minPt_prompt = -1000;
    unsigned  FirstMuon = -1;
    for(unsigned i=0; i<mu_isTightMuon->size(); i++){
      if (mu_isTightMuon->at(i)==0. || deltaBeta->at(i)>isoCut || mu_pt->at(i) < 24 || abs(mu_eta->at(i)) > 2.4) continue;  
      if (mu_pt->at(i) > minPt_prompt){
	minPt_prompt=mu_pt->at(i);
	FirstMuon=i;	  
      }
    }	
    //here I save only the info of the prompt muon -- for the sv muon we can clone the branch as it is..
    Float_t   minPt_second = -1000;
    unsigned  SecondMuon = -1;
    int count = 0;
    for(unsigned i=0; i<mu_isLoose->size(); i++){
      if(i == FirstMuon) continue;
      if(mu_isLoose->size()== 1) continue;
      if(mu_isLoose->at(i)==0 || mu_pt->at(i) < 5  || abs(mu_eta->at(i)) > 2.4) continue;
      count++;
      if (mu_pt->at(i)>minPt_second){
	minPt_second=mu_pt->at(i);
	SecondMuon=i;
      }
    }    
    unsigned  SecondVertex = -1;
    for(unsigned i=0; i<sv_mu_pt->size(); i++){
      if(sv_mu_pt->at(i) <  0) continue;
      SecondVertex=i;
    }

    if(FirstMuon != -1 && SecondMuon != -1 && SecondVertex != -1){

      //pile up weight
      pu_weight = PU_Weight->at(pu);

     //prompt muon 
    mu_promptPt = mu_pt->at(FirstMuon);
    mu_promptEta = mu_eta->at(FirstMuon);
    mu_promptPhi = mu_phi->at(FirstMuon);
    mu_promptCharge = mu_charge->at(FirstMuon);
    mu_promptE = mu_en->at(FirstMuon);
    mu_promptEt = mu_et->at(FirstMuon);
    mu_promptRhoIso   = mu_rhoIso->at(FirstMuon);
    mu_promptTrackiso = mu_trackiso->at(FirstMuon);
    mu_promptPfSumChHadPt  = mu_pfSumChargedHadronPt->at(FirstMuon);   
    mu_promptPfSumNHadEt   = mu_pfSumNeutralHadronEt->at(FirstMuon);  
    mu_promptPFSumPhotonEt = mu_PFSumPhotonEt->at(FirstMuon);  
    mu_promptPfSumPUPt = mu_pfSumPUPt->at(FirstMuon);   
    mu_promptEmIso     = mu_emIso->at(FirstMuon);  
    mu_promptHadIso    = mu_hadIso->at(FirstMuon); 
    mu_promptNormalizedChi2 = mu_normalizedChi2->at(FirstMuon);
    mu_promptDPToverPT      = mu_dPToverPTTunePMuonBestTrack->at(FirstMuon); 
    mu_promptAbsdxy      = mu_absdxyTunePMuonBestTrack->at(FirstMuon);   
    mu_promptAbsdxyError = mu_absdxyErrorTunePMuonBestTrack->at(FirstMuon);
    mu_promptAbsdxySig   = mu_absdxySigTunePMuonBestTrack->at(FirstMuon); 
    mu_promptAbsdz = mu_absdzTunePMuonBestTrack->at(FirstMuon);  
    mu_promptAbsdz = mu_absdzErrorTunePMuonBestTrack->at(FirstMuon);
    mu_promptAbsdzSig      = mu_absdzSigTunePMuonBestTrack->at(FirstMuon);  
    mu_promptRecoDeltaBeta = deltaBeta->at(FirstMuon);  
    mu_promptRecoiso   = mu_recoiso->at(FirstMuon);  

    mu_promptDirection = mu_STATofDirection->at(FirstMuon); 
    mu_promptNDof      = mu_STATofNDof->at(FirstMuon);  
    mu_promptTimeAtIpInOut    = mu_STATofTimeAtIpInOut->at(FirstMuon);
    mu_promptTimeAtIpInOutErr = mu_STATofTimeAtIpInOutErr->at(FirstMuon);
    mu_promptTimeAtIpOutIn    = mu_STATofTimeAtIpOutIn->at(FirstMuon); 
    mu_promptTimeAtIpOutInErr = mu_STATofTimeAtIpOutInErr->at(FirstMuon);

    mu_promptMatchedStations  = mu_numberOfMatchedStations->at(FirstMuon); 
    mu_promptValidPixelHits   = mu_numberOfValidPixelHits->at(FirstMuon);  
    mu_promptTrackQuality     = mu_TrackQuality->at(FirstMuon);   
    mu_promptInrTrackQuality  = mu_InnerTrackQuality->at(FirstMuon);
    mu_promptValidMuonHits    = mu_numberOfValidMuonHits->at(FirstMuon); 
    mu_promptGenMatch         = mu_FirstGenMatch->at(FirstMuon); 

    //non_prompt muon
    mu_secondPt = mu_pt->at(SecondMuon);
    mu_secondEta = mu_eta->at(SecondMuon);
    mu_secondPhi = mu_phi->at(SecondMuon);
    mu_secondCharge = mu_charge->at(SecondMuon);
    mu_secondE = mu_en->at(SecondMuon);
    mu_secondEt = mu_et->at(SecondMuon);
    mu_secondRhoIso   = mu_rhoIso->at(SecondMuon);
    mu_secondTrackiso = mu_trackiso->at(SecondMuon);
    mu_secondPfSumChHadPt  = mu_pfSumChargedHadronPt->at(SecondMuon);   
    mu_secondPfSumNHadEt   = mu_pfSumNeutralHadronEt->at(SecondMuon);  
    mu_secondPFSumPhotonEt = mu_PFSumPhotonEt->at(SecondMuon);  
    mu_secondPfSumPUPt = mu_pfSumPUPt->at(SecondMuon);   
    mu_secondEmIso     = mu_emIso->at(SecondMuon);  
    mu_secondHadIso    = mu_hadIso->at(SecondMuon); 
    mu_secondNormalizedChi2 = mu_normalizedChi2->at(SecondMuon);
    mu_secondDPToverPT      = mu_dPToverPTTunePMuonBestTrack->at(SecondMuon); 
    mu_secondAbsdxy      = mu_absdxyTunePMuonBestTrack->at(SecondMuon);   
    mu_secondAbsdxyError = mu_absdxyErrorTunePMuonBestTrack->at(SecondMuon);
    mu_secondAbsdxySig   = mu_absdxySigTunePMuonBestTrack->at(SecondMuon); 
    mu_secondAbsdz = mu_absdzTunePMuonBestTrack->at(SecondMuon);  
    mu_secondAbsdz = mu_absdzErrorTunePMuonBestTrack->at(SecondMuon);
    mu_secondAbsdzSig      = mu_absdzSigTunePMuonBestTrack->at(SecondMuon);  
    mu_secondRecoDeltaBeta = deltaBeta->at(SecondMuon);  
    mu_secondRecoiso   = mu_recoiso->at(SecondMuon);  

    mu_secondDirection = mu_STATofDirection->at(SecondMuon); 
    mu_secondNDof      = mu_STATofNDof->at(SecondMuon);  
    mu_secondTimeAtIpInOut    = mu_STATofTimeAtIpInOut->at(SecondMuon);
    mu_secondTimeAtIpInOutErr = mu_STATofTimeAtIpInOutErr->at(SecondMuon);
    mu_secondTimeAtIpOutIn    = mu_STATofTimeAtIpOutIn->at(SecondMuon); 
    mu_secondTimeAtIpOutInErr = mu_STATofTimeAtIpOutInErr->at(SecondMuon);

    mu_secondMatchedStations  = mu_numberOfMatchedStations->at(SecondMuon); 
    mu_secondValidPixelHits   = mu_numberOfValidPixelHits->at(SecondMuon);  
    mu_secondTrackQuality     = mu_TrackQuality->at(SecondMuon);   
    mu_secondInrTrackQuality  = mu_InnerTrackQuality->at(SecondMuon);
    mu_secondValidMuonHits    = mu_numberOfValidMuonHits->at(SecondMuon); 
    mu_secondGenMatch         = mu_SecondGenMatch->at(SecondMuon);


    float R = sqrt(((mu_secondEta-mu_promptEta)*(mu_secondEta-mu_promptEta))+((mu_secondPhi-mu_promptPhi)*(mu_secondPhi-mu_promptPhi)));
    float charged =  mu_pfSumChargedHadronPt->at(SecondMuon);
    float neutral =  mu_pfSumNeutralHadronEt->at(SecondMuon);
    float sumPhotonEt = mu_PFSumPhotonEt->at(SecondMuon);
    float pileup =  mu_pfSumPUPt->at(SecondMuon); 

    double deltaBetaR3 = (charged + std::max(0.0, neutral+sumPhotonEt-0.5*pileup))/mu_secondPt;

    TLorentzVector Mu1;
    TLorentzVector Mu2; 
    Mu1.SetPtEtaPhiE(mu_promptPt,mu_promptEta,mu_promptPhi,mu_promptE);
    Mu2.SetPtEtaPhiE(mu_secondPt,mu_secondEta,mu_secondPhi,mu_secondE);
    float DiMuMass = (Mu1 + Mu2).M();

    mu_DeltaBetaR3 = deltaBetaR3;
    mu_DiMuMass = DiMuMass;
    mu_Size = count;
    mu_DeltaR = R ;

    //secondary vertex info   
    sv_TrackSize =  sv_mu_TrackSize->at(SecondVertex);
    sv_LXYSig =  sv_mu_LXYSig->at(SecondVertex);
    sv_LXYZSig =  sv_mu_LXYZSig->at(SecondVertex);
    sv_LXY =  sv_mu_LXY->at(SecondVertex);
    sv_LXYZ =  sv_mu_LXYZ->at(SecondVertex);
    sv_mass =  sv_mu_mass->at(SecondVertex);
    sv_eta =  sv_mu_eta->at(SecondVertex);
    sv_phi =  sv_mu_phi->at(SecondVertex);
    sv_pt =  sv_mu_pt->at(SecondVertex);
    sv_p =  sv_mu_p->at(SecondVertex);
    sv_Beta =  sv_mu_Beta->at(SecondVertex);
    sv_Gamma =  sv_mu_Gamma->at(SecondVertex);
    sv_CTau0 =  sv_mu_CTau0->at(SecondVertex);
    sv_NDof =  sv_mu_NDof->at(SecondVertex);
    sv_Chi2 =  sv_mu_Chi2->at(SecondVertex);
    sv_Angle3D =  sv_mu_Angle3D->at(SecondVertex);
    sv_Angle2D =  sv_mu_Angle2D->at(SecondVertex);
    sv_tracks_Sumcharge =  sv_mu_tracks_Sumcharge->at(SecondVertex);
    sv_tracks_Sumpt =  sv_mu_tracks_Sumpt->at(SecondVertex);
    sv_match =  sv_mu_match->at(SecondVertex);

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

    //bjet info 
    unsigned  bjet = -1;
    Float_t bjet_disc = 1000;
    for(unsigned i=0; i<jet_btag_pfCSVv2IVF_discriminator->size(); i++){
      if(jet_btag_pfCSVv2IVF_discriminator->at(i) < bjet_disc) {
	bjet_disc = jet_btag_pfCSVv2IVF_discriminator->at(i);
        bjet = i;
      }
    }
      jet_btag_pt_  = jet_btag_pt->at(bjet);
      jet_btag_eta_ = jet_btag_eta->at(bjet);
      jet_btag_phi_ = jet_btag_phi->at(bjet);
      jet_btag_flavor_ = jet_btag_flavor->at(bjet);
      jet_btag_discriminator_ = jet_btag_pfCSVv2IVF_discriminator->at(bjet);

    newtree->Fill();
    }

  }
  
  newtree->Print();
  newtree->AutoSave();

  delete newfile;

  return 0;
}

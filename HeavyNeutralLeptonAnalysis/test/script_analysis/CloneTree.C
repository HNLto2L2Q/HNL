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


int main(int argc, char** argv){

  if(argc < 4) {
    std::cout<<"Arguments Missing!"<<std::endl;
    std::cout<<"Usage:  "<<argv[0]<<"  path/to/ntuples  input_filename  isMC"<<std::endl;
    return 0;
  }
  //Get old file, old tree and set top branch address
  //TFile *oldfile = new TFile("/eos/cms/store/group/phys_exotica/HNL/Background/crab_Analysis_WZToLLLNu/Background_Analysis.root");
  //TFile *oldfile = new TFile("/user/moanwar/Run2016/CMSSW_8_0_29/src/HNL/HeavyNeutralLeptonAnalysis/test/signal/HNL_M5_mu_2.95.root");
  //TTree *oldtree = (TTree*)oldfile->Get("HeavyNeutralLepton/tree_");
  std::string filepath = argv[1];
  std::string filename = argv[2];
  std::string   flagMC = argv[3];

  if( (flagMC.find("true") == string::npos && flagMC.find("True") == string::npos) &&
      (flagMC.find("false") == string::npos && flagMC.find("False") == string::npos) ){
        std::cout<<" Flag for specifing data or MC processing should be only [true/True] or [false/False]\n";
        return 0;
      }
  //////selection cuts

  Float_t isoCut = 0.15;
  bool isMC = (flagMC.find("true") != string::npos || flagMC.find("True") != string::npos);
  //if I want to use a TChain.....
  //cout<< "starting..."<<endl;
  TChain * oldtree = new TChain("HeavyNeutralLepton/tree_","");
  string tree ="HeavyNeutralLepton/tree_";
  //Background

  oldtree->Add(Form("%s/%s/%s",filepath.c_str(),filename.c_str(),tree.c_str()));

//=============================================================================================//

  TH1F* h_nTrueInteractions50      = new TH1F("nTrueInteractions50" , "nTrueInteractions 50 bin" , 50, 0., 50. );
  TH1F* h_nTrueInteractions100     = new TH1F("nTrueInteractions100" , "nTrueInteractions 100 bin" , 100, 0., 100. );

  TFile *newfile = new TFile(filename.c_str(),"recreate");

  //TFile *newfile = new TFile("skimmedSignale.root","recreate");
  //Create a new file + a clone of old tree in new file 
  //TTree *newtree = oldtree->CloneTree(0); 
  TTree *newtree  = new TTree("newtree","Analysis Tree");
  //cout<<"cloning done"<<endl;

  // Long64_t nentries = oldtree->GetEntries();
  Long64_t nentries = oldtree->GetEntries();
  cout  <<nentries<<endl;

//======================= Old Tree Variables ==========================================// 
  // These are the variables I cut on 
  Bool_t passIsoMu24All;
  Bool_t passIsoMu27All;
  Bool_t passMu17      ;
  Bool_t passMu17_Mu8_SameSign;
  Bool_t passMu20             ;
  Bool_t passMu17_Mu8         ;
  Bool_t passMu27_TkMu8       ;

  oldtree->SetBranchAddress("passIsoMu24All",&passIsoMu24All);
  oldtree->SetBranchAddress("passIsoMu27All",&passIsoMu27All);

  oldtree->SetBranchAddress("passMu17",&passMu17);
  oldtree->SetBranchAddress("passMu17_Mu8_SameSign",&passMu17_Mu8_SameSign);
  oldtree->SetBranchAddress("passMu20",&passMu20);
  oldtree->SetBranchAddress("passMu17_Mu8",&passMu17_Mu8);
  oldtree->SetBranchAddress("passMu27_TkMu8",&passMu27_TkMu8);

  vector<Float_t>   *PU_Weight = 0;
  vector<Float_t>   *PU_WeightUp = 0;
  vector<Float_t>   *PU_WeightDown = 0;
  vector<Float_t>   *npT = 0;

  oldtree->SetBranchAddress("PU_Weight",&PU_Weight);
  oldtree->SetBranchAddress("PU_WeightUp",&PU_WeightUp);
  oldtree->SetBranchAddress("PU_WeightDown",&PU_WeightDown);
  oldtree->SetBranchAddress("npT",&npT);

  Float_t pvX ;
  Float_t pvY ;
  Float_t pvZ ;
  Float_t pvXErr ;
  Float_t pvYErr ;
  Float_t pvZErr ;
  Float_t pvLxy  ;
  Float_t pvLxyz ;
  Float_t pvLxySig  ;
  Float_t pvLxyzSig ;
  Float_t pvChi2 ;
  Float_t pvSumPtSq ;
  Int_t numberPV;
  Int_t pvNTrack ;

  oldtree->SetBranchAddress("pvX" , &pvX); 
  oldtree->SetBranchAddress("pvY" , &pvY); 
  oldtree->SetBranchAddress("pvZ" , &pvZ); 
  oldtree->SetBranchAddress("pvXErr" , &pvXErr);
  oldtree->SetBranchAddress("pvYErr" , &pvYErr);
  oldtree->SetBranchAddress("pvZErr" , &pvZErr);
  oldtree->SetBranchAddress("pvLxy" , &pvLxy);
  oldtree->SetBranchAddress("pvLxyz" , &pvLxyz);
  oldtree->SetBranchAddress("pvLxySig" , &pvLxySig);
  oldtree->SetBranchAddress("pvLxyzSig" , &pvLxyzSig);
  oldtree->SetBranchAddress("pvChi2" , &pvChi2);
  oldtree->SetBranchAddress("pvNTrack" , &pvNTrack);
  oldtree->SetBranchAddress("pvSumPtSq" , &pvSumPtSq);
  oldtree->SetBranchAddress("numberPV" , &numberPV);

  vector<Float_t>   *mu_isTightMuon = 0;
  vector<Float_t>   *mu_isLoose = 0;
  vector<Float_t>   *mu_ptTunePMuonBestTrack = 0;
  vector<Float_t>   *mu_etaTunePMuonBestTrack =0;
  vector<Float_t>   *mu_phiTunePMuonBestTrack=0; 
  vector<Float_t>   *mu_chargeTunePMuonBestTrack=0; 
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
  vector<Float_t>   *mu_pxTunePMuonBestTrack = 0;
  vector<Float_t>   *mu_pyTunePMuonBestTrack = 0;
  vector<Float_t>   *mu_pzTunePMuonBestTrack = 0;
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
  vector<Int_t>     *mu_numberOfpixelLayersWithMeasurement = 0;
  vector<Int_t>     *mu_numberOftrackerLayersWithMeasurement = 0;
  vector<Int_t>     *mu_TrackQuality = 0 ;
  vector<Int_t>     *mu_InnerTrackQuality = 0 ;
  vector<Int_t>     *mu_FirstGenMatch = 0;
  vector<Int_t>     *mu_SecondGenMatch = 0;

  vector<Float_t>   *mu_InnerTrackValidFraction = 0;
  vector<Float_t>   *mu_segmentCompatibilityMuonBestTrack = 0;
  vector<Float_t>   *mu_trkKinkMuonBestTrack = 0;
  vector<Float_t>   *mu_chi2LocalPositionMuonBestTrack = 0;
  vector<Float_t>   *mu_isGlobalMuon = 0 ;

  vector<Float_t>   *mu_RPCTofDirection  = 0 ;
  vector<Float_t>   *mu_RPCTofNDof  = 0 ;
  vector<Float_t>   *mu_RPCTofTimeAtIpInOut  = 0 ;
  vector<Float_t>   *mu_RPCTofTimeAtIpInOutErr  = 0 ;
  vector<Float_t>   *mu_RPCTofTimeAtIpOutIn  = 0 ;
  vector<Float_t>   *mu_RPCTofTimeAtIpOutInErr  = 0 ;

  vector<Int_t>   *sv_mu_TrackSize = 0 ;
  vector<Float_t> *sv_mu_LXYSig = 0 ;
  vector<Float_t> *sv_mu_LXYZSig = 0 ;
  vector<Float_t> *sv_mu_LXY = 0 ;
  vector<Float_t> *sv_mu_LXYZ = 0 ;
  vector<Float_t> *sv_mu_dir_x = 0 ;
  vector<Float_t> *sv_mu_dir_y = 0 ;
  vector<Float_t> *sv_mu_dir_z = 0 ;
  vector<Float_t> *sv_mu_mass = 0 ;
  vector<Float_t> *sv_mu_eta = 0 ;
  vector<Float_t> *sv_mu_phi = 0 ;
  vector<Float_t> *sv_mu_pt = 0 ;
  vector<Float_t> *sv_mu_p = 0 ;
  vector<Float_t> *sv_mu_px = 0 ;
  vector<Float_t> *sv_mu_py = 0 ;
  vector<Float_t> *sv_mu_pz = 0 ;
  vector<Float_t> *sv_mu_energy = 0 ;
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

  vector<float> *sv_mu_Xpos = 0;
  vector<float> *sv_mu_Ypos = 0;
  vector<float> *sv_mu_Zpos = 0;
  vector<float> *sv_mu_xError = 0;
  vector<float> *sv_mu_yError = 0;
  vector<float> *sv_mu_zError = 0;
  vector<float> *sv_mu_pvX = 0;
  vector<float> *sv_mu_pvY = 0;
  vector<float> *sv_mu_pvZ = 0;
  vector<float> *sv_mu_pvXError = 0;
  vector<float> *sv_mu_pvYError = 0;
  vector<float> *sv_mu_pvZError = 0;
  // this wrong but leave it for now 
  vector<vector<int> >    *sv_mu_tracks_charge = 0;
  vector<vector<float> >  *sv_mu_tracks_eta = 0;
  vector<vector<float> >  *sv_mu_tracks_phi = 0;
  vector<vector<float> >  *sv_mu_tracks_pt  = 0;
  vector<vector<float> >  *sv_mu_tracks_dxySig = 0;
  vector<vector<float> >  *sv_mu_tracks_dxy = 0;
  vector<vector<float> >  *sv_mu_tracks_dxyz = 0;
  vector<vector<float> >  *sv_mu_tracks_en = 0;

  // this the correct 
  //vector<vector<float  > > *sv_mu_tracks_pt  = 0;

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
  vector<Float_t>   *jet_chargedMultiplicity = 0;
  vector<Float_t>   *jet_chargedEmEnergyFraction = 0;
  vector<Float_t>   *jet_chargedHadronEnergyFraction = 0;

  vector<Float_t>   *jet_CsvV2 = 0;
  vector<Float_t>   *jet_DeepCsv_udsg = 0;
  vector<Float_t>   *jet_DeepCsv_b = 0;
  vector<Float_t>   *jet_DeepCsv_c = 0;
  vector<Float_t>   *jet_DeepCsv_bb = 0;
  vector<Float_t>   *jet_HadronFlavor = 0;


  Float_t   pfMet_et;
  Float_t   pfMet_pt;
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
  vector<Float_t>   *jet_btag_pfCSVv2IVF_discriminator = 0;
  */
 
  oldtree->SetBranchAddress("mu_isTightMuon",&mu_isTightMuon);
  oldtree->SetBranchAddress("mu_isLoose",&mu_isLoose);
  oldtree->SetBranchAddress("mu_ptTunePMuonBestTrack",&mu_ptTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_etaTunePMuonBestTrack",&mu_etaTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_phiTunePMuonBestTrack",&mu_phiTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_chargeTunePMuonBestTrack",&mu_chargeTunePMuonBestTrack);
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
  oldtree->SetBranchAddress("mu_pxTunePMuonBestTrack" , &mu_pxTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_pyTunePMuonBestTrack" , &mu_pyTunePMuonBestTrack);
  oldtree->SetBranchAddress("mu_pzTunePMuonBestTrack" , &mu_pzTunePMuonBestTrack);
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
  oldtree->SetBranchAddress("mu_numberOfpixelLayersWithMeasurement", &mu_numberOfpixelLayersWithMeasurement);
  oldtree->SetBranchAddress("mu_numberOftrackerLayersWithMeasurement", &mu_numberOftrackerLayersWithMeasurement);
  oldtree->SetBranchAddress("mu_InnerTrackQuality", &mu_InnerTrackQuality);
  oldtree->SetBranchAddress("mu_InnerTrackValidFraction", &mu_InnerTrackValidFraction);
  oldtree->SetBranchAddress("mu_segmentCompatibilityMuonBestTrack", &mu_segmentCompatibilityMuonBestTrack);
  oldtree->SetBranchAddress("mu_trkKinkMuonBestTrack", &mu_trkKinkMuonBestTrack);
  oldtree->SetBranchAddress("mu_chi2LocalPositionMuonBestTrack", &mu_chi2LocalPositionMuonBestTrack);
  oldtree->SetBranchAddress("mu_isGlobalMuon", &mu_isGlobalMuon);

  oldtree->SetBranchAddress("mu_RPCTofDirection" , &mu_RPCTofDirection);
  oldtree->SetBranchAddress("mu_RPCTofNDof" , &mu_RPCTofNDof);
  oldtree->SetBranchAddress("mu_RPCTofTimeAtIpInOut" , &mu_RPCTofTimeAtIpInOut);
  oldtree->SetBranchAddress("mu_RPCTofTimeAtIpInOutErr" , &mu_RPCTofTimeAtIpInOutErr);
  oldtree->SetBranchAddress("mu_RPCTofTimeAtIpOutIn" , &mu_RPCTofTimeAtIpOutIn);
  oldtree->SetBranchAddress("mu_RPCTofTimeAtIpOutInErr" , &mu_RPCTofTimeAtIpOutInErr);

  oldtree->SetBranchAddress("mu_FirstGenMatch", &mu_FirstGenMatch);
  oldtree->SetBranchAddress("mu_SecondGenMatch", &mu_SecondGenMatch);

 
  oldtree->SetBranchAddress("sv_mu_TrackSize", &sv_mu_TrackSize);
  oldtree->SetBranchAddress("sv_mu_LXYSig", &sv_mu_LXYSig);
  oldtree->SetBranchAddress("sv_mu_LXYZSig", &sv_mu_LXYZSig);
  oldtree->SetBranchAddress("sv_mu_LXY", &sv_mu_LXY);
  oldtree->SetBranchAddress("sv_mu_LXYZ", &sv_mu_LXYZ);
  oldtree->SetBranchAddress("sv_mu_dir_x", &sv_mu_dir_x);
  oldtree->SetBranchAddress("sv_mu_dir_y", &sv_mu_dir_y);
  oldtree->SetBranchAddress("sv_mu_dir_z", &sv_mu_dir_z);
  oldtree->SetBranchAddress("sv_mu_mass", &sv_mu_mass);
  oldtree->SetBranchAddress("sv_mu_eta", &sv_mu_eta);   
  oldtree->SetBranchAddress("sv_mu_phi", &sv_mu_phi);  
  oldtree->SetBranchAddress("sv_mu_pt", &sv_mu_pt);  
  oldtree->SetBranchAddress("sv_mu_p", &sv_mu_p);  
  oldtree->SetBranchAddress("sv_mu_px", &sv_mu_px);
  oldtree->SetBranchAddress("sv_mu_py", &sv_mu_py);
  oldtree->SetBranchAddress("sv_mu_pz", &sv_mu_pz);
  oldtree->SetBranchAddress("sv_mu_energy", &sv_mu_energy);
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

  oldtree->SetBranchAddress("sv_mu_Xpos" , &sv_mu_Xpos);
  oldtree->SetBranchAddress("sv_mu_Ypos" , &sv_mu_Ypos);
  oldtree->SetBranchAddress("sv_mu_Zpos" , &sv_mu_Zpos);
  oldtree->SetBranchAddress("sv_mu_xError" , &sv_mu_xError);
  oldtree->SetBranchAddress("sv_mu_yError" , &sv_mu_yError);
  oldtree->SetBranchAddress("sv_mu_zError" , &sv_mu_zError);
  oldtree->SetBranchAddress("sv_mu_pvX" , &sv_mu_pvX);
  oldtree->SetBranchAddress("sv_mu_pvY" , &sv_mu_pvY);
  oldtree->SetBranchAddress("sv_mu_pvZ" , &sv_mu_pvZ);
  oldtree->SetBranchAddress("sv_mu_pvXError" , &sv_mu_pvXError);
  oldtree->SetBranchAddress("sv_mu_pvYError" , &sv_mu_pvYError);
  oldtree->SetBranchAddress("sv_mu_pvZError" , &sv_mu_pvZError);


  oldtree->SetBranchAddress("sv_mu_tracks_charge" , &sv_mu_tracks_charge);
  oldtree->SetBranchAddress("sv_mu_tracks_eta" , &sv_mu_tracks_eta);
  oldtree->SetBranchAddress("sv_mu_tracks_phi" , &sv_mu_tracks_phi);
  oldtree->SetBranchAddress("sv_mu_tracks_pt" , &sv_mu_tracks_pt);
  oldtree->SetBranchAddress("sv_mu_tracks_dxySig" , &sv_mu_tracks_dxySig);
  oldtree->SetBranchAddress("sv_mu_tracks_dxy" , &sv_mu_tracks_dxy);
  oldtree->SetBranchAddress("sv_mu_tracks_dxyz" , &sv_mu_tracks_dxyz);
  oldtree->SetBranchAddress("sv_mu_tracks_en" , &sv_mu_tracks_en);

  /*
  oldtree->SetBranchAddress("jet_btag_pt", &jet_btag_pt);
  oldtree->SetBranchAddress("jet_btag_eta", &jet_btag_eta);
  oldtree->SetBranchAddress("jet_btag_phi", &jet_btag_phi);
  oldtree->SetBranchAddress("jet_btag_flavor", &jet_btag_flavor);
  oldtree->SetBranchAddress("jet_btag_pfCSVv2IVF_discriminator", &jet_btag_pfCSVv2IVF_discriminator);
  */
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
  oldtree->SetBranchAddress("jet_chargedMultiplicity", &jet_chargedMultiplicity);
  oldtree->SetBranchAddress("jet_chargedHadronEnergyFraction" , &jet_chargedHadronEnergyFraction);
  oldtree->SetBranchAddress("jet_chargedEmEnergyFraction" ,&jet_chargedEmEnergyFraction);

  oldtree->SetBranchAddress("jet_CsvV2" ,&jet_CsvV2);
  oldtree->SetBranchAddress("jet_DeepCsv_udsg" ,&jet_DeepCsv_udsg);
  oldtree->SetBranchAddress("jet_DeepCsv_b" ,&jet_DeepCsv_b);
  oldtree->SetBranchAddress("jet_DeepCsv_c" ,&jet_DeepCsv_c);
  oldtree->SetBranchAddress("jet_DeepCsv_bb" ,&jet_DeepCsv_bb);
  oldtree->SetBranchAddress("jet_HadronFlavor" ,&jet_HadronFlavor);

  //Throw away some branches -- thos are all the  one I don't want to keep in the ntuples
  //IMPORTANT: we cannot throw away the branches with the variable we are using in the loop!
  //oldtree->SetBranchStatus("ele_*",    0);
  //oldtree->SetBranchStatus("mu_ST*", 0);
//================================ trigger info ============================================//   
  Int_t trig_IsoMu24All,trig_IsoMu27All,trig_Mu17,trig_Mu17_Mu8_SameSign,trig_Mu20,trig_Mu17_Mu8,trig_Mu27_TkMu8;

  TBranch* branch_trig_IsoMu24All          = newtree->Branch("trig_IsoMu24All", &trig_IsoMu24All, "trig_IsoMu24All/I");
  TBranch* branch_trig_IsoMu27All          = newtree->Branch("trig_IsoMu27All", &trig_IsoMu27All, "trig_IsoMu27All/I");
  TBranch* branch_trig_Mu17                = newtree->Branch("trig_Mu17", &trig_Mu17, "trig_Mu17/I");
  TBranch* branch_trig_Mu17_Mu8_SameSign   = newtree->Branch("trig_Mu17_Mu8_SameSign", &trig_Mu17_Mu8_SameSign, "trig_Mu17_Mu8_SameSign/I");
  TBranch* branch_trig_Mu20                = newtree->Branch("trig_Mu20", &trig_Mu20, "trig_Mu20/I");
  TBranch* branch_trig_Mu17_Mu8            = newtree->Branch("trig_Mu17_Mu8", &trig_Mu17_Mu8, "trig_Mu17_Mu8/I");
  TBranch* branch_trig_Mu27_TkMu8          = newtree->Branch("trig_Mu27_TkMu8", &trig_Mu27_TkMu8, "trig_Mu27_TkMu8/I");


//================================ pv info ============================================//
  Float_t pvX_ ,pvY_ , pvZ_ , pvXErr_ , pvYErr_,
         pvZErr_ , pvLxy_  , pvLxyz_ , pvLxySig_  ,
         pvLxyzSig_ , pvChi2_,  pvSumPtSq_ ;
  Int_t   pvNTrack_, numberPV_ ;

  TBranch* branch_pvX_        = newtree->Branch("pvX_", &pvX_, "pvX_/F");
  TBranch* branch_pvY_        = newtree->Branch("pvY_", &pvY_, "pvY_/F");
  TBranch* branch_pvZ_        = newtree->Branch("pvZ_", &pvZ_, "pvZ_/F");
  TBranch* branch_pvXErr_     = newtree->Branch("pvXErr_", &pvXErr_, "pvXErr_/F");
  TBranch* branch_pvYErr_     = newtree->Branch("pvYErr_", &pvYErr_, "pvYErr_/F");
  TBranch* branch_pvZErr_     = newtree->Branch("pvZErr_", &pvZErr_, "pvZErr_/F");
  TBranch* branch_pvLxy_      = newtree->Branch("pvLxy_", &pvLxy_, "pvLxy_/F");
  TBranch* branch_pvLxyz_     = newtree->Branch("pvLxyz_", &pvLxyz_, "pvLxyz_/F");
  TBranch* branch_pvLxySig_   = newtree->Branch("pvLxySig_", &pvLxySig_, "pvLxySig_/F");
  TBranch* branch_pvLxyzSig_  = newtree->Branch("pvLxyzSig_", &pvLxyzSig_, "pvLxyzSig_/F");
  TBranch* branch_pvChi2_     = newtree->Branch("pvChi2_", &pvChi2_, "pvChi2_/F");
  TBranch* branch_pvNTrack_   = newtree->Branch("pvNTrack_", &pvNTrack_, "pvNTrack/I");
  TBranch* branch_pvSumPtSq_  = newtree->Branch("pvSumPtSq_", &pvSumPtSq_, "pvSumPtSq_/F");
  TBranch* branch_numberPV_   = newtree->Branch("numberPV_", &numberPV_, "numberPV/I");

//================================ PU Weight ============================================// 
  Float_t pu_weight, pu_weightUp, pu_weightDown,TrueNumInteractions;

  TBranch* branch_pu_weight               = newtree->Branch("pu_weight",&pu_weight,"pu_weight/F");
  TBranch* branch_pu_weightUp             = newtree->Branch("pu_weightUp",&pu_weightUp,"pu_weightUp/F");
  TBranch* branch_pu_weightDown           = newtree->Branch("pu_weightDown",&pu_weightDown,"pu_weightDown/F");
  TBranch* branch_TrueNumInteractions     = newtree->Branch("TrueNumInteractions",&TrueNumInteractions,"TrueNumInteractions/F");

//======================= First Muon Variables ==========================================//   
Float_t mu_promptPt,mu_promptEta, mu_promptPhi, mu_promptCharge, mu_promptEt, mu_promptE;

Float_t mu_promptRhoIso,        mu_promptTrackiso,       mu_promptPfSumChHadPt,
	mu_promptPFSumPhotonEt, mu_promptPfSumPUPt,      mu_promptEmIso,
	mu_promptHadIso,        mu_promptNormalizedChi2, mu_promptDPToverPT,
	mu_promptAbsdxy,        mu_promptAbsdxyError,    mu_promptAbsdxySig,
	mu_promptAbsdz,         mu_promptAbsdzError,     mu_promptAbsdzSig,
	mu_promptRecoDeltaBeta, mu_promptRecoiso,        mu_promptDirection,
	mu_promptNDof,          mu_promptTimeAtIpInOut,  mu_promptTimeAtIpInOutErr,
        mu_promptTimeAtIpOutIn, mu_promptTimeAtIpOutInErr,  mu_promptPfSumNHadEt, mu_promptGlobalMuon;

Float_t mu_promptInnerTrackFraction, mu_promptSegmentCompatibility,
        mu_promptTrkKink,    mu_promptChi2LocalPosition;

 Int_t   mu_promptValidMuonHits, mu_promptMatchedStations, mu_promptValidPixelHits, 
         mu_promptTrackQuality,  mu_promptInrTrackQuality, mu_promptGenMatch,
         mu_promptPixelLayers ,  mu_promptTrackerLayers;

Float_t mu_promptRPCTofDirection,        mu_promptRPCTofNDof,          mu_promptRPCTofTimeAtIpInOut ,
        mu_promptRPCTofTimeAtIpInOutErr, mu_promptRPCTofTimeAtIpOutIn, mu_promptRPCTofTimeAtIpOutInErr;


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
TBranch* branch_mu_promptPixelLayers     = newtree->Branch("mu_promptPixelLayers", &mu_promptPixelLayers , "mu_promptPixelLayers/I");
TBranch* branch_mu_promptTrackerLayers   = newtree->Branch("mu_promptTrackerLayers", &mu_promptTrackerLayers, "mu_promptTrackerLayers/I");

TBranch* branch_mu_promptGenMatch = newtree->Branch("mu_promptGenMatch", &mu_promptGenMatch, "mu_promptGenMatch/I");

TBranch* branch_mu_promptInnerTrackFraction = newtree->Branch("mu_promptInnerTrackFraction", &mu_promptInnerTrackFraction, "mu_promptInnerTrackFraction/F");
TBranch* branch_mu_promptSegmentCompatibility = newtree->Branch("mu_promptSegmentCompatibility", &mu_promptSegmentCompatibility, "mu_promptSegmentCompatibility/F");
TBranch* branch_mu_promptTrkKink = newtree->Branch("mu_promptTrkKink", &mu_promptTrkKink, "mu_promptTrkKink/F");
TBranch* branch_mu_promptChi2LocalPosition= newtree->Branch("mu_promptChi2LocalPosition",&mu_promptChi2LocalPosition,"mu_promptChi2LocalPosition/F");

TBranch* branch_mu_promptGlobalMuon = newtree->Branch("mu_promptGlobalMuon", &mu_promptGlobalMuon , "mu_promptGlobalMuon/F");

 TBranch* branch_mu_promptRPCTofDirection = newtree->Branch("mu_promptRPCTofDirection", &mu_promptRPCTofDirection , "mu_promptRPCTofDirection/F");
 TBranch* branch_mu_promptRPCTofNDof      = newtree->Branch("mu_promptRPCTofNDof", &mu_promptRPCTofNDof , "mu_promptRPCTofNDof/F");
 TBranch* branch_mu_promptRPCTofTimeAtIpInOut    = newtree->Branch("mu_promptRPCTofTimeAtIpInOut", &mu_promptRPCTofTimeAtIpInOut , "mu_promptRPCTofTimeAtIpInOut/F");
 TBranch* branch_mu_promptRPCTofTimeAtIpInOutErr = newtree->Branch("mu_promptRPCTofTimeAtIpInOutErr", &mu_promptRPCTofTimeAtIpInOutErr , "mu_promptRPCTofTimeAtIpInOutErr/F");
 TBranch* branch_mu_promptRPCTofTimeAtIpOutIn    = newtree->Branch("mu_promptRPCTofTimeAtIpOutIn", &mu_promptRPCTofTimeAtIpOutIn , "mu_promptRPCTofTimeAtIpOutIn/F");
 TBranch* branch_mu_promptRPCTofTimeAtIpOutInErr = newtree->Branch("mu_promptRPCTofTimeAtIpOutInErr", &mu_promptRPCTofTimeAtIpOutInErr , "mu_promptRPCTofTimeAtIpOutInErr/F");

//======================= Second Muon Variables ==========================================//
Float_t mu_secondPt,mu_secondEta, mu_secondPhi, mu_secondCharge, mu_secondEt, mu_secondE;

Float_t mu_secondRhoIso,        mu_secondTrackiso,  mu_secondPfSumChHadPt,
        mu_secondPFSumPhotonEt, mu_secondPfSumPUPt,      mu_secondEmIso,        
        mu_secondHadIso,        mu_secondNormalizedChi2, mu_secondDPToverPT,   
        mu_secondAbsdxy,        mu_secondAbsdxyError,    mu_secondAbsdxySig, 
        mu_secondAbsdz,         mu_secondAbsdzError,     mu_secondAbsdzSig,
        mu_secondRecoDeltaBeta, mu_secondRecoiso,        mu_secondDirection,
        mu_secondNDof,          mu_secondTimeAtIpInOut,  mu_secondTimeAtIpInOutErr, 
        mu_secondTimeAtIpOutIn, mu_secondTimeAtIpOutInErr,  mu_secondPfSumNHadEt, mu_secondGlobalMuon;

Float_t mu_secondInnerTrackFraction, mu_secondSegmentCompatibility,
        mu_secondTrkKink,            mu_secondChi2LocalPosition;


 Float_t mu_DeltaBetaR3, mu_DiMuMass, mu_Size, mu_DeltaR, mu_mT, vtxmu_mass;

 Int_t   mu_secondValidMuonHits, mu_secondMatchedStations, mu_secondValidPixelHits, 
         mu_secondTrackQuality,  mu_secondInrTrackQuality, mu_secondGenMatch,
         mu_secondPixelLayers ,  mu_secondTrackerLayers, mu_secondIsTight , mu_nbLoose;

 Float_t mu_secondRPCTofDirection,        mu_secondRPCTofNDof,          mu_secondRPCTofTimeAtIpInOut ,
         mu_secondRPCTofTimeAtIpInOutErr, mu_secondRPCTofTimeAtIpOutIn, mu_secondRPCTofTimeAtIpOutInErr;

TBranch* branch_mu_secondIsTight = newtree->Branch("mu_secondIsTight",&mu_secondIsTight,"mu_secondIsTight/I");
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
TBranch* branch_mu_secondPixelLayers     = newtree->Branch("mu_secondPixelLayers",   &mu_secondPixelLayers ,  "mu_secondPixelLayers/I");
TBranch* branch_mu_secondTrackerLayers   = newtree->Branch("mu_secondTrackerLayers", &mu_secondTrackerLayers, "mu_secondTrackerLayers/I");
TBranch* branch_mu_secondGenMatch        = newtree->Branch("mu_secondGenMatch", &mu_secondGenMatch, "mu_secondGenMatch/I");
TBranch* branch_mu_secondInnerTrackFraction = newtree->Branch("mu_secondInnerTrackFraction", &mu_secondInnerTrackFraction, "mu_secondInnerTrackFraction/F");
TBranch* branch_mu_secondSegmentCompatibility = newtree->Branch("mu_secondSegmentCompatibility", &mu_secondSegmentCompatibility, "mu_secondSegmentCompatibility/F");
TBranch* branch_mu_secondGlobalMuon = newtree->Branch("mu_secondGlobalMuon", &mu_secondGlobalMuon , "mu_secondGlobalMuon/F");
TBranch* branch_mu_secondTrkKink = newtree->Branch("mu_secondTrkKink", &mu_secondTrkKink, "mu_secondTrkKink/F");
TBranch* branch_mu_secondChi2LocalPosition= newtree->Branch("mu_secondChi2LocalPosition",&mu_secondChi2LocalPosition,"mu_secondChi2LocalPosition/F");

TBranch* branch_mu_secondRPCTofDirection = newtree->Branch("mu_secondRPCTofDirection", &mu_secondRPCTofDirection , "mu_secondRPCTofDirection/F");
TBranch* branch_mu_secondRPCTofNDof      = newtree->Branch("mu_secondRPCTofNDof", &mu_secondRPCTofNDof , "mu_secondRPCTofNDof/F");
TBranch* branch_mu_secondRPCTofTimeAtIpInOut    = newtree->Branch("mu_secondRPCTofTimeAtIpInOut", &mu_secondRPCTofTimeAtIpInOut , "mu_secondRPCTofTimeAtIpInOut/F");
TBranch* branch_mu_secondRPCTofTimeAtIpInOutErr = newtree->Branch("mu_secondRPCTofTimeAtIpInOutErr", &mu_secondRPCTofTimeAtIpInOutErr , "mu_secondRPCTofTimeAtIpInOutErr/F");
TBranch* branch_mu_secondRPCTofTimeAtIpOutIn    = newtree->Branch("mu_secondRPCTofTimeAtIpOutIn", &mu_secondRPCTofTimeAtIpOutIn , "mu_secondRPCTofTimeAtIpOutIn/F");
TBranch* branch_mu_secondRPCTofTimeAtIpOutInErr = newtree->Branch("mu_secondRPCTofTimeAtIpOutInErr", &mu_secondRPCTofTimeAtIpOutInErr , "mu_secondRPCTofTimeAtIpOutInErr/F");

TBranch* branch_mu_DeltaBetaR3 = newtree->Branch("mu_DeltaBetaR3", &mu_DeltaBetaR3, "mu_DeltaBetaR3/F");
TBranch* branch_mu_DiMuMass    = newtree->Branch("mu_DiMuMass", &mu_DiMuMass, "mu_DiMuMass/F");
TBranch* branch_mu_DeltaR      = newtree->Branch("mu_DeltaR", &mu_DeltaR, "mu_DeltaR/F");
TBranch* branch_mu_Size        = newtree->Branch("mu_Size", &mu_Size, "mu_Size/F");
TBranch* branch_mu_nbLoose     = newtree->Branch("mu_nbLoose", &mu_nbLoose, "mu_nbLoose/F");
TBranch* branch_mu_mT          = newtree->Branch("mu_mT", &mu_mT, "mu_mT/F");
TBranch* branch_vtxmu_mass     = newtree->Branch("vtxmu_mass", &vtxmu_mass, "vtxmu_mass/F");

//======================= prompt Jet Variables ==========================================// 
 Float_t promptJet_charge_, promptJet_et_,    promptJet_pt_, promptJet_eta_, promptJet_NEmEnFraction_,
         promptJet_phi_,    promptJet_theta_, promptJet_en_, promptJet_chEmEn_, promptJet_NHadEnFraction_,
         promptJet_chMuEn_,         promptJet_chMuEnFraction_, promptJet_numberOfDaughters_, 
         promptJet_muonEnergyFraction_, promptJet_muonMultiplicity_, 
         promptJet_neutralEmEnergy_,  promptJet_neutralHadronEnergy_, promptJet_NHadMultiplicity_, 
         promptJet_NMultiplicity_, promptJet_chHadEn_, promptJet_muonEnergy_, promptJet_chargedMultiplicity_ , 
         promptJet_chargedEmEnergyFraction, promptJet_chargedHadronEnergyFraction;

 Float_t promptJet_CsvV2, promptJet_DeepCsv_udsg, promptJet_DeepCsv_b, promptJet_DeepCsv_c, promptJet_DeepCsv_bb, promptJet_HadronFlavor;

 TBranch* branch_promptJet_charge_  = newtree->Branch("promptJet_charge_", &promptJet_charge_, "promptJet_charge_/F");
 TBranch* branch_promptJet_et_      = newtree->Branch("promptJet_et_", &promptJet_et_, "promptJet_et_/F");
 TBranch* branch_promptJet_pt_      = newtree->Branch("promptJet_pt_", &promptJet_pt_, "promptJet_pt_/F");
 TBranch* branch_promptJet_eta_     = newtree->Branch("promptJet_eta_", &promptJet_eta_, "promptJet_eta_/F");
 TBranch* branch_promptJet_phi_     = newtree->Branch("promptJet_phi_", &promptJet_phi_, "promptJet_phi_/F");
 TBranch* branch_promptJet_theta_   = newtree->Branch("promptJet_theta_", &promptJet_theta_, "promptJet_theta_/F");
 TBranch* branch_promptJet_en_      = newtree->Branch("promptJet_en_", &promptJet_en_, "promptJet_en_/F");
 TBranch* branch_promptJet_chargedEmEnergy_  = newtree->Branch("promptJet_chEmEn_", &promptJet_chEmEn_, "promptJet_chEmEn_/F");
 TBranch* branch_promptJet_NEmEnFraction_    = newtree->Branch("promptJet_NEmEnFraction_", &promptJet_NEmEnFraction_, "promptJet_NEmEnFraction_/F");
 TBranch* branch_promptJet_chHadEn_          = newtree->Branch("promptJet_chHadEn_", &promptJet_chHadEn_, "promptJet_chHadEn_/F");
 TBranch* branch_promptJet_NHadEnFraction_   = newtree->Branch("promptJet_NHadEnFraction_", &promptJet_NHadEnFraction_, "promptJet_NHadEnFraction_/F");
 TBranch* branch_promptJet_chMuEn_          = newtree->Branch("promptJet_chMuEn_", &promptJet_chMuEn_, "promptJet_chMuEn_/F");
 TBranch* branch_promptJet_chMuEnFraction_  = newtree->Branch("promptJet_chMuEnFraction_", &promptJet_chMuEnFraction_, "promptJet_chMuEnFraction_/F");
 TBranch* branch_promptJet_numberOfDaughters_  = newtree->Branch("promptJet_numberOfDaughters_", &promptJet_numberOfDaughters_, "promptJet_numberOfDaughters_/F");
 TBranch* branch_promptJet_muonEnergy_           = newtree->Branch("promptJet_muonEnergy_", &promptJet_muonEnergy_, "promptJet_muonEnergy_/F");
 TBranch* branch_promptJet_muonEnergyFraction_   = newtree->Branch("promptJet_muonEnergyFraction_", &promptJet_muonEnergyFraction_, "promptJet_muonEnergyFraction_/F");
 TBranch* branch_promptJet_muonMultiplicity_     = newtree->Branch("promptJet_muonMultiplicity_", &promptJet_muonMultiplicity_, "promptJet_muonMultiplicity_/F");
 TBranch* branch_promptJet_neutralEmEnergy_      = newtree->Branch("promptJet_neutralEmEnergy_", &promptJet_neutralEmEnergy_, "promptJet_neutralEmEnergy_/F");
 TBranch* branch_promptJet_neutralHadronEnergy_  = newtree->Branch("promptJet_neutralHadronEnergy_", &promptJet_neutralHadronEnergy_, "promptJet_neutralHadronEnergy_/F");
 TBranch* branch_promptJet_NHadMultiplicity_  = newtree->Branch("promptJet_NHadMultiplicity_", &promptJet_NHadMultiplicity_, "promptJet_NHadMultiplicity_/F");
 TBranch* branch_promptJet_NMultiplicity_     = newtree->Branch("promptJet_NMultiplicity_", &promptJet_NMultiplicity_, "promptJet_NMultiplicity_/F");
 TBranch* branch_promptJet_chargedMultiplicity_    = newtree->Branch("promptJet_chargedMultiplicity_", &promptJet_chargedMultiplicity_, "promptJet_chargedMultiplicity_/F");
 TBranch* branch_promptJet_chargedEmEnergyFraction =  newtree->Branch("promptJet_chargedEmEnergyFraction", &promptJet_chargedEmEnergyFraction, "promptJet_chargedEmEnergyFraction/F");
 TBranch* branch_promptJet_chargedHadronEnergyFraction = newtree->Branch("promptJet_chargedHadronEnergyFraction", &promptJet_chargedHadronEnergyFraction, "promptJet_chargedHadronEnergyFraction/F");

 TBranch* branch_promptJet_CsvV2        =  newtree->Branch("promptJet_CsvV2", &promptJet_CsvV2, "promptJet_CsvV2/F");
 TBranch* branch_promptJet_DeepCsv_udsg =  newtree->Branch("promptJet_DeepCsv_udsg", &promptJet_DeepCsv_udsg, "promptJet_DeepCsv_udsg/F");
 TBranch* branch_promptJet_DeepCsv_b    =  newtree->Branch("promptJet_DeepCsv_b", &promptJet_DeepCsv_b, "promptJet_DeepCsv_b/F");
 TBranch* branch_promptJet_DeepCsv_c    =  newtree->Branch("promptJet_DeepCsv_c", &promptJet_DeepCsv_c,"promptJet_DeepCsv_c/F");
 TBranch* branch_promptJet_DeepCsv_bb   =  newtree->Branch("promptJet_DeepCsv_bb", &promptJet_DeepCsv_bb,"promptJet_DeepCsv_bb/F");
 TBranch* branch_promptJet_HadronFlavor =  newtree->Branch("promptJet_HadronFlavor", &promptJet_HadronFlavor,"promptJet_HadronFlavor/F");

 //======================= non prompt Jet Variables ==========================================// 
 Float_t non_promptJet_charge_, non_promptJet_et_,    non_promptJet_pt_, non_promptJet_eta_, non_promptJet_NEmEnFraction_,
         non_promptJet_phi_,    non_promptJet_theta_, non_promptJet_en_, non_promptJet_chEmEn_, non_promptJet_NHadEnFraction_,
         non_promptJet_chMuEn_,         non_promptJet_chMuEnFraction_, non_promptJet_numberOfDaughters_, 
         non_promptJet_muonEnergyFraction_, non_promptJet_muonMultiplicity_, 
         non_promptJet_neutralEmEnergy_,  non_promptJet_neutralHadronEnergy_, non_promptJet_NHadMultiplicity_, 
         non_promptJet_NMultiplicity_, non_promptJet_chHadEn_, non_promptJet_muonEnergy_, non_promptJet_chargedMultiplicity_ , 
         non_promptJet_chargedEmEnergyFraction, non_promptJet_chargedHadronEnergyFraction;

 Float_t non_promptJet_CsvV2,     non_promptJet_DeepCsv_udsg, non_promptJet_DeepCsv_b, 
         non_promptJet_DeepCsv_c, non_promptJet_DeepCsv_bb,   non_promptJet_HadronFlavor;


 Float_t  DiJet_Mass, MuJet_Mass;

 TBranch* branch_non_promptJet_charge_  = newtree->Branch("non_promptJet_charge_", &non_promptJet_charge_, "non_promptJet_charge_/F");
 TBranch* branch_non_promptJet_et_      = newtree->Branch("non_promptJet_et_", &non_promptJet_et_, "non_promptJet_et_/F");
 TBranch* branch_non_promptJet_pt_      = newtree->Branch("non_promptJet_pt_", &non_promptJet_pt_, "non_promptJet_pt_/F");
 TBranch* branch_non_promptJet_eta_     = newtree->Branch("non_promptJet_eta_", &non_promptJet_eta_, "non_promptJet_eta_/F");
 TBranch* branch_non_promptJet_phi_     = newtree->Branch("non_promptJet_phi_", &non_promptJet_phi_, "non_promptJet_phi_/F");
 TBranch* branch_non_promptJet_theta_   = newtree->Branch("non_promptJet_theta_", &non_promptJet_theta_, "non_promptJet_theta_/F");
 TBranch* branch_non_promptJet_en_      = newtree->Branch("non_promptJet_en_", &non_promptJet_en_, "non_promptJet_en_/F");
 TBranch* branch_non_promptJet_chargedEmEnergy_  = newtree->Branch("non_promptJet_chEmEn_", &non_promptJet_chEmEn_, "non_promptJet_chEmEn_/F");
 TBranch* branch_non_promptJet_NEmEnFraction_    = newtree->Branch("non_promptJet_NEmEnFraction_", &non_promptJet_NEmEnFraction_, "non_promptJet_NEmEnFraction_/F");
 TBranch* branch_non_promptJet_chHadEn_         = newtree->Branch("non_promptJet_chHadEn_", &non_promptJet_chHadEn_, "non_promptJet_chHadEn_/F");
 TBranch* branch_non_promptJet_NHadEnFraction_  = newtree->Branch("non_promptJet_NHadEnFraction_", &non_promptJet_NHadEnFraction_, "non_promptJet_NHadEnFraction_/F");
 TBranch* branch_non_promptJet_chMuEn_          = newtree->Branch("non_promptJet_chMuEn_", &non_promptJet_chMuEn_, "non_promptJet_chMuEn_/F");
 TBranch* branch_non_promptJet_chMuEnFraction_  = newtree->Branch("non_promptJet_chMuEnFraction_", &non_promptJet_chMuEnFraction_, "non_promptJet_chMuEnFraction_/F");
 TBranch* branch_non_promptJet_numberOfDaughters_    = newtree->Branch("non_promptJet_numberOfDaughters_", &non_promptJet_numberOfDaughters_, "non_promptJet_numberOfDaughters_/F");
 TBranch* branch_non_promptJet_muonEnergy_           = newtree->Branch("non_promptJet_muonEnergy_", &non_promptJet_muonEnergy_, "non_promptJet_muonEnergy_/F");
 TBranch* branch_non_promptJet_muonEnergyFraction_   = newtree->Branch("non_promptJet_muonEnergyFraction_", &non_promptJet_muonEnergyFraction_, "non_promptJet_muonEnergyFraction_/F");
 TBranch* branch_non_promptJet_muonMultiplicity_     = newtree->Branch("non_promptJet_muonMultiplicity_", &non_promptJet_muonMultiplicity_, "non_promptJet_muonMultiplicity_/F");
 TBranch* branch_non_promptJet_neutralEmEnergy_      = newtree->Branch("non_promptJet_neutralEmEnergy_", &non_promptJet_neutralEmEnergy_, "non_promptJet_neutralEmEnergy_/F");
 TBranch* branch_non_promptJet_neutralHadronEnergy_  = newtree->Branch("non_promptJet_neutralHadronEnergy_", &non_promptJet_neutralHadronEnergy_, "non_promptJet_neutralHadronEnergy_/F");
 TBranch* branch_non_promptJet_NHadMultiplicity_  = newtree->Branch("non_promptJet_NHadMultiplicity_", &non_promptJet_NHadMultiplicity_, "non_promptJet_NHadMultiplicity_/F");
 TBranch* branch_non_promptJet_NMultiplicity_     = newtree->Branch("non_promptJet_NMultiplicity_", &non_promptJet_NMultiplicity_, "non_promptJet_NMultiplicity_/F");
 TBranch* branch_non_promptJet_chargedMultiplicity_    = newtree->Branch("non_promptJet_chargedMultiplicity_", &non_promptJet_chargedMultiplicity_, "non_promptJet_chargedMultiplicity_/F");
 TBranch* branch_non_promptJet_chargedEmEnergyFraction =  newtree->Branch("non_promptJet_chargedEmEnergyFraction", &non_promptJet_chargedEmEnergyFraction, "non_promptJet_chargedEmEnergyFraction/F");
 TBranch* branch_non_promptJet_chargedHadronEnergyFraction = newtree->Branch("non_promptJet_chargedHadronEnergyFraction", &non_promptJet_chargedHadronEnergyFraction, "non_promptJet_chargedHadronEnergyFraction/F");

 TBranch* branch_non_promptJet_CsvV2        =  newtree->Branch("non_promptJet_CsvV2", &non_promptJet_CsvV2, "non_promptJet_CsvV2/F");
 TBranch* branch_non_promptJet_DeepCsv_udsg =  newtree->Branch("non_promptJet_DeepCsv_udsg", &non_promptJet_DeepCsv_udsg, "non_promptJet_DeepCsv_udsg/F");
 TBranch* branch_non_promptJet_DeepCsv_b    =  newtree->Branch("non_promptJet_DeepCsv_b", &non_promptJet_DeepCsv_b, "non_promptJet_DeepCsv_b/F");
 TBranch* branch_non_promptJet_DeepCsv_c    =  newtree->Branch("non_promptJet_DeepCsv_c", &non_promptJet_DeepCsv_c,"non_promptJet_DeepCsv_c/F");
 TBranch* branch_non_promptJet_DeepCsv_bb   =  newtree->Branch("non_promptJet_DeepCsv_bb", &non_promptJet_DeepCsv_bb,"non_promptJet_DeepCsv_bb/F");
 TBranch* branch_non_promptJet_HadronFlavor =  newtree->Branch("non_promptJet_HadronFlavor", &non_promptJet_HadronFlavor,"non_promptJet_HadronFlavor/F");

 TBranch* branch_DiJet_Mass = newtree->Branch("DiJet_Mass", &DiJet_Mass ,"DiJet_Mass/F");
 TBranch* branch_MuJet_Mass = newtree->Branch("MuJet_Mass", &MuJet_Mass ,"MuJet_Mass/F");


//======================= Second Vertex Variables ==========================================//
 Float_t  sv_LXYSig,  sv_LXYZSig, sv_LXY,  sv_LXYZ,  sv_mass, 
          sv_eta,     sv_phi,     sv_pt,   sv_p,     sv_Beta, 
          sv_Gamma,   sv_CTau0,   sv_NDof, sv_Chi2,  sv_Angle3D, 
          sv_Angle2D, sv_tracks_Sumpt,     sv_match, sv_energy,
          sv_px,      sv_py,      sv_pz,   sv_dir_x, sv_dir_y, sv_dir_z;

 Float_t sv_Xpos,  sv_Ypos,  sv_Zpos,  sv_xError,   sv_yError,   sv_zError,
         sv_pvX,   sv_pvY,   sv_pvZ,   sv_pvXError, sv_pvYError, sv_pvZError;
 Float_t sv_track_sumdxySig;

 Int_t sv_TrackSize , sv_tracks_Sumcharge ;  
 unsigned sv_inside;

 TBranch* branch_sv_inside     = newtree->Branch("sv_inside", &sv_inside, "sv_inside/I");
 TBranch* branch_sv_TrackSize  = newtree->Branch("sv_TrackSize", &sv_TrackSize, "sv_TrackSize/I");
 TBranch* branch_sv_LXYSig     = newtree->Branch("sv_LXYSig", &sv_LXYSig, "sv_LXYSig/F");
 TBranch* branch_sv_LXYZSig    = newtree->Branch("sv_LXYZSig", &sv_LXYZSig, "sv_LXYZSig/F");
 TBranch* branch_sv_LXY        = newtree->Branch("sv_LXY", &sv_LXY, "sv_LXY/F");
 TBranch* branch_sv_LXYZ       = newtree->Branch("sv_LXYZ", &sv_LXYZ, "sv_LXYZ/F");
 TBranch* branch_sv_dir_x      = newtree->Branch("sv_dir_x", &sv_dir_x, "sv_dir_x/F");
 TBranch* branch_sv_dir_y      = newtree->Branch("sv_dir_y", &sv_dir_y, "sv_dir_y/F");
 TBranch* branch_sv_dir_z      = newtree->Branch("sv_dir_z", &sv_dir_z, "sv_dir_z/F");
 TBranch* branch_sv_mass       = newtree->Branch("sv_mass", &sv_mass, "sv_mass/F");
 TBranch* branch_sv_eta        = newtree->Branch("sv_eta", &sv_eta, "sv_eta/F");
 TBranch* branch_sv_phi        = newtree->Branch("sv_phi", &sv_phi, "sv_phi/F");
 TBranch* branch_sv_pt         = newtree->Branch("sv_pt", &sv_pt, "sv_pt/F");
 TBranch* branch_sv_p          = newtree->Branch("sv_p", &sv_p, "sv_p/F");
 TBranch* branch_sv_px         = newtree->Branch("sv_px", &sv_px, "sv_px/F");
 TBranch* branch_sv_py         = newtree->Branch("sv_py", &sv_py, "sv_py/F");
 TBranch* branch_sv_pz         = newtree->Branch("sv_pz", &sv_pz, "sv_pz/F");
 TBranch* branch_sv_energy     = newtree->Branch("sv_energy", &sv_energy, "sv_energy/F");
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

 TBranch* branch_sv_Xpos     = newtree->Branch("sv_Xpos" , &sv_Xpos,"sv_Xpos/F");
 TBranch* branch_sv_Ypos     = newtree->Branch("sv_Ypos" , &sv_Ypos,"sv_Ypos/F");
 TBranch* branch_sv_Zpos     = newtree->Branch("sv_Zpos" , &sv_Zpos,"sv_Zpos/F");
 TBranch* branch_sv_xError   = newtree->Branch("sv_xError" , &sv_xError,"sv_xError/F");
 TBranch* branch_sv_yError   = newtree->Branch("sv_yError" , &sv_yError,"sv_yError/F");
 TBranch* branch_sv_zError   = newtree->Branch("sv_zError" , &sv_zError,"sv_zError/F");
 TBranch* branch_sv_pvX      = newtree->Branch("sv_pvX" , &sv_pvX,"sv_pvX/F");
 TBranch* branch_sv_pvY      = newtree->Branch("sv_pvY" , &sv_pvY,"sv_pvY/F");
 TBranch* branch_sv_pvZ      = newtree->Branch("sv_pvZ" , &sv_pvZ,"sv_pvZ/F");
 TBranch* branch_sv_pvXError = newtree->Branch("sv_pvXError" , &sv_pvXError,"sv_pvXError/F");
 TBranch* branch_sv_pvYError = newtree->Branch("sv_pvYError" , &sv_pvYError,"sv_pvYError/F");
 TBranch* branch_sv_pvZError = newtree->Branch("sv_pvZError" , &sv_pvZError,"sv_pvZError/F");

 TBranch* branch_sv_track_sumdxySig = newtree->Branch("sv_track_sumdxySig" , &sv_track_sumdxySig, "sv_track_sumdxySig/F");
 //================================= First Track in sv  ================================================// 
  Float_t  firstTrack_eta , firstTrack_phi , firstTrack_pt , firstTrack_dxySig , firstTrack_dxy , firstTrack_dxyz , firstTrack_en; 
  Int_t    firstTrack_charge;

  TBranch* branch_firstTrack_eta    = newtree->Branch("firstTrack_eta" , &firstTrack_eta,"firstTrack_eta/F");
  TBranch* branch_firstTrack_phi    = newtree->Branch("firstTrack_phi" , &firstTrack_phi,"firstTrack_phi/F");
  TBranch* branch_firstTrack_pt     = newtree->Branch("firstTrack_pt" , &firstTrack_pt,"firstTrack_pt/F");
  TBranch* branch_firstTrack_dxySig = newtree->Branch("firstTrack_dxySig" , &firstTrack_dxySig,"firstTrack_dxySig/F");
  TBranch* branch_firstTrack_dxy    = newtree->Branch("firstTrack_dxy" , &firstTrack_dxy,"firstTrack_dxy/F");
  TBranch* branch_firstTrack_dxyz   = newtree->Branch("firstTrack_dxyz" , &firstTrack_dxyz,"firstTrack_dxyz/F");
  TBranch* branch_firstTrack_charge = newtree->Branch("firstTrack_charge" , &firstTrack_charge,"firstTrack_charge/I");
  TBranch* branch_firstTrack_en     = newtree->Branch("firstTrack_en" , &firstTrack_en,"firstTrack_en/F");
 //================================= Second Track in sv ================================================// 
  Float_t  secondTrack_eta , secondTrack_phi , secondTrack_pt , secondTrack_dxySig , secondTrack_dxy , secondTrack_dxyz ,  secondTrack_en;
  Int_t    secondTrack_charge;

  TBranch* branch_secondTrack_eta    = newtree->Branch("secondTrack_eta" , &secondTrack_eta,"secondTrack_eta/F");
  TBranch* branch_secondTrack_phi    = newtree->Branch("secondTrack_phi" , &secondTrack_phi,"secondTrack_phi/F");
  TBranch* branch_secondTrack_pt     = newtree->Branch("secondTrack_pt" , &secondTrack_pt,"secondTrack_pt/F");
  TBranch* branch_secondTrack_dxySig = newtree->Branch("secondTrack_dxySig" , &secondTrack_dxySig,"secondTrack_dxySig/F");
  TBranch* branch_secondTrack_dxy    = newtree->Branch("secondTrack_dxy" , &secondTrack_dxy,"secondTrack_dxy/F");
  TBranch* branch_secondTrack_dxyz   = newtree->Branch("secondTrack_dxyz" , &secondTrack_dxyz,"secondTrack_dxyz/F");
  TBranch* branch_secondTrack_charge = newtree->Branch("secondTrack_charge" , &secondTrack_charge,"secondTrack_charge/I");
  TBranch* branch_secondTrack_en     = newtree->Branch("secondTrack_en" , &secondTrack_en,"secondTrack_en/F");
 //================================= First Jet Variables ==============================================//
 Float_t FirstJet_charge_, FirstJet_et_,    FirstJet_pt_, FirstJet_eta_, FirstJet_NEmEnFraction_,
         FirstJet_phi_,    FirstJet_theta_, FirstJet_en_, FirstJet_chEmEn_, FirstJet_NHadEnFraction_,
         FirstJet_chMuEn_,         FirstJet_chMuEnFraction_, FirstJet_numberOfDaughters_, 
         FirstJet_muonEnergyFraction_, FirstJet_muonMultiplicity_, 
         FirstJet_neutralEmEnergy_,  FirstJet_neutralHadronEnergy_, FirstJet_NHadMultiplicity_, 
         FirstJet_NMultiplicity_, FirstJet_chHadEn_, FirstJet_muonEnergy_, FirstJet_chargedMultiplicity_ , 
         FirstJet_chargedEmEnergyFraction, FirstJet_chargedHadronEnergyFraction;

 Float_t FirstJet_CsvV2, FirstJet_DeepCsv_udsg, FirstJet_DeepCsv_b, FirstJet_DeepCsv_c, FirstJet_DeepCsv_bb, FirstJet_HadronFlavor;


 Int_t   jets_size, bjet1_size, bjet2_size, bjet3_size, bjet4_size, bjet5_size, bjet6_size, bjet7_size, bjet8_size, bjet9_size, bjet_L,bjet_M,bjet_T;

 TBranch* branch_FirstJet_charge_  = newtree->Branch("FirstJet_charge_", &FirstJet_charge_, "FirstJet_charge_/F");
 TBranch* branch_FirstJet_et_  = newtree->Branch("FirstJet_et_", &FirstJet_et_, "FirstJet_et_/F");
 TBranch* branch_FirstJet_pt_  = newtree->Branch("FirstJet_pt_", &FirstJet_pt_, "FirstJet_pt_/F");
 TBranch* branch_FirstJet_eta_  = newtree->Branch("FirstJet_eta_", &FirstJet_eta_, "FirstJet_eta_/F");
 TBranch* branch_FirstJet_phi_  = newtree->Branch("FirstJet_phi_", &FirstJet_phi_, "FirstJet_phi_/F");
 TBranch* branch_FirstJet_theta_  = newtree->Branch("FirstJet_theta_", &FirstJet_theta_, "FirstJet_theta_/F");
 TBranch* branch_FirstJet_en_  = newtree->Branch("FirstJet_en_", &FirstJet_en_, "FirstJet_en_/F");
 TBranch* branch_FirstJet_chargedEmEnergy_  = newtree->Branch("FirstJet_chEmEn_", &FirstJet_chEmEn_, "FirstJet_chEmEn_/F");
 TBranch* branch_FirstJet_NEmEnFraction_  = newtree->Branch("FirstJet_NEmEnFraction_", &FirstJet_NEmEnFraction_, "FirstJet_NEmEnFraction_/F");
 TBranch* branch_FirstJet_chHadEn_  = newtree->Branch("FirstJet_chHadEn_", &FirstJet_chHadEn_, "FirstJet_chHadEn_/F");
 TBranch* branch_FirstJet_NHadEnFraction_  = newtree->Branch("FirstJet_NHadEnFraction_", &FirstJet_NHadEnFraction_, "FirstJet_NHadEnFraction_/F");
 TBranch* branch_FirstJet_chMuEn_  = newtree->Branch("FirstJet_chMuEn_", &FirstJet_chMuEn_, "FirstJet_chMuEn_/F");
 TBranch* branch_FirstJet_chMuEnFraction_  = newtree->Branch("FirstJet_chMuEnFraction_", &FirstJet_chMuEnFraction_, "FirstJet_chMuEnFraction_/F");
 TBranch* branch_FirstJet_numberOfDaughters_  = newtree->Branch("FirstJet_numberOfDaughters_", &FirstJet_numberOfDaughters_, "FirstJet_numberOfDaughters_/F");
 TBranch* branch_FirstJet_muonEnergy_  = newtree->Branch("FirstJet_muonEnergy_", &FirstJet_muonEnergy_, "FirstJet_muonEnergy_/F");
 TBranch* branch_FirstJet_muonEnergyFraction_  = newtree->Branch("FirstJet_muonEnergyFraction_", &FirstJet_muonEnergyFraction_, "FirstJet_muonEnergyFraction_/F");
 TBranch* branch_FirstJet_muonMultiplicity_  = newtree->Branch("FirstJet_muonMultiplicity_", &FirstJet_muonMultiplicity_, "FirstJet_muonMultiplicity_/F");
 TBranch* branch_FirstJet_neutralEmEnergy_  = newtree->Branch("FirstJet_neutralEmEnergy_", &FirstJet_neutralEmEnergy_, "FirstJet_neutralEmEnergy_/F");
 TBranch* branch_FirstJet_neutralHadronEnergy_  = newtree->Branch("FirstJet_neutralHadronEnergy_", &FirstJet_neutralHadronEnergy_, "FirstJet_neutralHadronEnergy_/F");
 TBranch* branch_FirstJet_NHadMultiplicity_  = newtree->Branch("FirstJet_NHadMultiplicity_", &FirstJet_NHadMultiplicity_, "FirstJet_NHadMultiplicity_/F");
 TBranch* branch_FirstJet_NMultiplicity_  = newtree->Branch("FirstJet_NMultiplicity_", &FirstJet_NMultiplicity_, "FirstJet_NMultiplicity_/F");
 TBranch* branch_FirstJet_chargedMultiplicity_   = newtree->Branch("FirstJet_chargedMultiplicity_", &FirstJet_chargedMultiplicity_, "FirstJet_chargedMultiplicity_/F");
 TBranch* branch_FirstJet_chargedEmEnergyFraction =  newtree->Branch("FirstJet_chargedEmEnergyFraction", &FirstJet_chargedEmEnergyFraction, "FirstJet_chargedEmEnergyFraction/F");
 TBranch* branch_FirstJet_chargedHadronEnergyFraction = newtree->Branch("FirstJet_chargedHadronEnergyFraction", &FirstJet_chargedHadronEnergyFraction, "FirstJet_chargedHadronEnergyFraction/F");

 TBranch* branch_jets_size =  newtree->Branch("jets_size", &jets_size, "jets_size/I");
 TBranch* branch_bjet1_size =  newtree->Branch("bjet1_size", &bjet1_size, "bjet1_size/I");
 TBranch* branch_bjet2_size =  newtree->Branch("bjet2_size", &bjet2_size, "bjet2_size/I");
 TBranch* branch_bjet3_size =  newtree->Branch("bjet3_size", &bjet3_size, "bjet3_size/I");
 TBranch* branch_bjet4_size =  newtree->Branch("bjet4_size", &bjet4_size, "bjet4_size/I");
 TBranch* branch_bjet5_size =  newtree->Branch("bjet5_size", &bjet5_size, "bjet5_size/I");
 TBranch* branch_bjet6_size =  newtree->Branch("bjet6_size", &bjet6_size, "bjet6_size/I");
 TBranch* branch_bjet7_size =  newtree->Branch("bjet7_size", &bjet7_size, "bjet7_size/I");
 TBranch* branch_bjet8_size =  newtree->Branch("bjet8_size", &bjet8_size, "bjet8_size/I");
 TBranch* branch_bjet9_size =  newtree->Branch("bjet9_size", &bjet9_size, "bjet9_size/I");
 TBranch* branch_bjet_L =  newtree->Branch("bjet_L", &bjet_L, "bjet_L/I");
 TBranch* branch_bjet_M =  newtree->Branch("bjet_M", &bjet_M, "bjet_M/I");
 TBranch* branch_bjet_T =  newtree->Branch("bjet_T", &bjet_T, "bjet_T/I");

 TBranch* branch_FirstJet_CsvV2        =  newtree->Branch("FirstJet_CsvV2", &FirstJet_CsvV2, "FirstJet_CsvV2/F");
 TBranch* branch_FirstJet_DeepCsv_udsg =  newtree->Branch("FirstJet_DeepCsv_udsg", &FirstJet_DeepCsv_udsg, "FirstJet_DeepCsv_udsg/F");
 TBranch* branch_FirstJet_DeepCsv_b    =  newtree->Branch("FirstJet_DeepCsv_b", &FirstJet_DeepCsv_b, "FirstJet_DeepCsv_b/F");
 TBranch* branch_FirstJet_DeepCsv_c    =  newtree->Branch("FirstJet_DeepCsv_c", &FirstJet_DeepCsv_c,"FirstJet_DeepCsv_c/F");
 TBranch* branch_FirstJet_DeepCsv_bb   =  newtree->Branch("FirstJet_DeepCsv_bb", &FirstJet_DeepCsv_bb,"FirstJet_DeepCsv_bb/F");
 TBranch* branch_FirstJet_HadronFlavor =  newtree->Branch("FirstJet_HadronFlavor", &FirstJet_HadronFlavor,"FirstJet_HadronFlavor/F");

 //================================= Second Jet Variables ==============================================//
 Float_t SecondJet_charge_, SecondJet_et_,    SecondJet_pt_, SecondJet_eta_, SecondJet_NEmEnFraction_,
         SecondJet_phi_,    SecondJet_theta_, SecondJet_en_, SecondJet_chEmEn_, SecondJet_NHadEnFraction_,
         SecondJet_chMuEn_,         SecondJet_chMuEnFraction_, SecondJet_numberOfDaughters_, 
         SecondJet_muonEnergyFraction_, SecondJet_muonMultiplicity_, 
         SecondJet_neutralEmEnergy_,  SecondJet_neutralHadronEnergy_, SecondJet_NHadMultiplicity_, 
         SecondJet_NMultiplicity_, SecondJet_chHadEn_, SecondJet_muonEnergy_, SecondJet_chargedMultiplicity_ ,
         SecondJet_chargedEmEnergyFraction, SecondJet_chargedHadronEnergyFraction;

 Float_t SecondJet_CsvV2, SecondJet_DeepCsv_udsg, SecondJet_DeepCsv_b, SecondJet_DeepCsv_c, SecondJet_DeepCsv_bb, SecondJet_HadronFlavor;

 Float_t FSjet_DR, FjetFmu_DR, FjetSmu_DR, SjetFmu_DR, SjetSmu_DR, FSjet_M, FjetFmu_M, FjetSmu_M, SjetFmu_M, SjetSmu_M, FSjetFmu_M, FSjetSmu_M, FSjetFSmu_M;


 TBranch* branch_SecondJet_charge_  = newtree->Branch("SecondJet_charge_", &SecondJet_charge_, "SecondJet_charge_/F");
 TBranch* branch_SecondJet_et_  = newtree->Branch("SecondJet_et_", &SecondJet_et_, "SecondJet_et_/F");
 TBranch* branch_SecondJet_pt_  = newtree->Branch("SecondJet_pt_", &SecondJet_pt_, "SecondJet_pt_/F");
 TBranch* branch_SecondJet_eta_  = newtree->Branch("SecondJet_eta_", &SecondJet_eta_, "SecondJet_eta_/F");
 TBranch* branch_SecondJet_phi_  = newtree->Branch("SecondJet_phi_", &SecondJet_phi_, "SecondJet_phi_/F");
 TBranch* branch_SecondJet_theta_  = newtree->Branch("SecondJet_theta_", &SecondJet_theta_, "SecondJet_theta_/F");
 TBranch* branch_SecondJet_en_  = newtree->Branch("SecondJet_en_", &SecondJet_en_, "SecondJet_en_/F");
 TBranch* branch_SecondJet_chargedEmEnergy_  = newtree->Branch("SecondJet_chEmEn_", &SecondJet_chEmEn_, "SecondJet_chEmEn_/F");
 TBranch* branch_SecondJet_NEmEnFraction_  = newtree->Branch("SecondJet_NEmEnFraction_", &SecondJet_NEmEnFraction_, "SecondJet_NEmEnFraction_/F");
 TBranch* branch_SecondJet_chHadEn_  = newtree->Branch("SecondJet_chHadEn_", &SecondJet_chHadEn_, "SecondJet_chHadEn_/F");
 TBranch* branch_SecondJet_NHadEnFraction_  = newtree->Branch("SecondJet_NHadEnFraction_", &SecondJet_NHadEnFraction_, "SecondJet_NHadEnFraction_/F");
 TBranch* branch_SecondJet_chMuEn_  = newtree->Branch("SecondJet_chMuEn_", &SecondJet_chMuEn_, "SecondJet_chMuEn_/F");
 TBranch* branch_SecondJet_chMuEnFraction_  = newtree->Branch("SecondJet_chMuEnFraction_", &SecondJet_chMuEnFraction_, "SecondJet_chMuEnFraction_/F");
 TBranch* branch_SecondJet_numberOfDaughters_  = newtree->Branch("SecondJet_numberOfDaughters_", &SecondJet_numberOfDaughters_, "SecondJet_numberOfDaughters_/F");
 TBranch* branch_SecondJet_muonEnergy_  = newtree->Branch("SecondJet_muonEnergy_", &SecondJet_muonEnergy_, "SecondJet_muonEnergy_/F");
 TBranch* branch_SecondJet_muonEnergyFraction_  = newtree->Branch("SecondJet_muonEnergyFraction_", &SecondJet_muonEnergyFraction_, "SecondJet_muonEnergyFraction_/F");
 TBranch* branch_SecondJet_muonMultiplicity_  = newtree->Branch("SecondJet_muonMultiplicity_", &SecondJet_muonMultiplicity_, "SecondJet_muonMultiplicity_/F");
 TBranch* branch_SecondJet_neutralEmEnergy_  = newtree->Branch("SecondJet_neutralEmEnergy_", &SecondJet_neutralEmEnergy_, "SecondJet_neutralEmEnergy_/F");
 TBranch* branch_SecondJet_neutralHadronEnergy_  = newtree->Branch("SecondJet_neutralHadronEnergy_", &SecondJet_neutralHadronEnergy_, "SecondJet_neutralHadronEnergy_/F");
 TBranch* branch_SecondJet_NHadMultiplicity_  = newtree->Branch("SecondJet_NHadMultiplicity_", &SecondJet_NHadMultiplicity_, "SecondJet_NHadMultiplicity_/F");
 TBranch* branch_SecondJet_NMultiplicity_  = newtree->Branch("SecondJet_NMultiplicity_", &SecondJet_NMultiplicity_, "SecondJet_NMultiplicity_/F");

 TBranch* branch_SecondJet_chargedMultiplicity_   = newtree->Branch("SecondJet_chargedMultiplicity_", &SecondJet_chargedMultiplicity_, "SecondJet_chargedMultiplicity_/F");
 TBranch* branch_SecondJet_chargedEmEnergyFraction =  newtree->Branch("SecondJet_chargedEmEnergyFraction", &SecondJet_chargedEmEnergyFraction, "SecondJet_chargedEmEnergyFraction/F");
 TBranch* branch_SecondJet_chargedHadronEnergyFraction = newtree->Branch("SecondJet_chargedHadronEnergyFraction", &SecondJet_chargedHadronEnergyFraction,"SecondJet_chargedHadronEnergyFraction/F");

 TBranch* branch_SecondJet_CsvV2        =  newtree->Branch("SecondJet_CsvV2", &SecondJet_CsvV2, "SecondJet_CsvV2/F");
 TBranch* branch_SecondJet_DeepCsv_udsg =  newtree->Branch("SecondJet_DeepCsv_udsg", &SecondJet_DeepCsv_udsg, "SecondJet_DeepCsv_udsg/F");
 TBranch* branch_SecondJet_DeepCsv_b    =  newtree->Branch("SecondJet_DeepCsv_b", &SecondJet_DeepCsv_b, "SecondJet_DeepCsv_b/F");
 TBranch* branch_SecondJet_DeepCsv_c    =  newtree->Branch("SecondJet_DeepCsv_c", &SecondJet_DeepCsv_c,"SecondJet_DeepCsv_c/F");
 TBranch* branch_SecondJet_DeepCsv_bb   =  newtree->Branch("SecondJet_DeepCsv_bb", &SecondJet_DeepCsv_bb,"SecondJet_DeepCsv_bb/F");
 TBranch* branch_SecondJet_HadronFlavor =  newtree->Branch("SecondJet_HadronFlavor", &SecondJet_HadronFlavor,"SecondJet_HadronFlavor/F");


 TBranch* branch_FSjet_DR =  newtree->Branch("FSjet_DR", &FSjet_DR ,"FSjet_DR/F");
 TBranch* branch_FjetFmu_DR =  newtree->Branch("FjetFmu_DR", &FjetFmu_DR ,"FjetFmu_DR/F");
 TBranch* branch_FjetSmu_DR =  newtree->Branch("FjetSmu_DR", &FjetSmu_DR ,"FjetSmu_DR/F");
 TBranch* branch_SjetFmu_DR =  newtree->Branch("SjetFmu_DR", &SjetFmu_DR ,"SjetFmu_DR/F");
 TBranch* branch_SjetSmu_DR =  newtree->Branch("SjetSmu_DR", &SjetSmu_DR ,"SjetSmu_DR/F");   
 TBranch* branch_FSjet_M=  newtree->Branch("FSjet_M", &FSjet_M,"FSjet_M/F");   
 TBranch* branch_FjetFmu_M=  newtree->Branch("FjetFmu_M", &FjetFmu_M,"FjetFmu_M/F");
 TBranch* branch_FjetSmu_M=  newtree->Branch("FjetSmu_M", &FjetSmu_M,"FjetSmu_M/F");
 TBranch* branch_SjetFmu_M=  newtree->Branch("SjetFmu_M", &SjetFmu_M,"SjetFmu_M/F");
 TBranch* branch_SjetSmu_M=  newtree->Branch("SjetSmu_M", &SjetSmu_M,"SjetSmu_M/F"); 
 TBranch* branch_FSjetFmu_M=  newtree->Branch("FSjetFmu_M", &FSjetFmu_M,"FSjetFmu_M/F");
 TBranch* branch_FSjetSmu_M=  newtree->Branch("FSjetSmu_M", &FSjetSmu_M,"FSjetSmu_M/F");   
 TBranch* branch_FSjetFSmu_M=  newtree->Branch("FSjetFSmu_M", &FSjetFSmu_M,"FSjetFSmu_M/F");
  
 //================================= Thrid Jet Variables ==============================================//
 Float_t ThirdJet_charge_, ThirdJet_et_,    ThirdJet_pt_, ThirdJet_eta_, ThirdJet_NEmEnFraction_,
         ThirdJet_phi_,    ThirdJet_theta_, ThirdJet_en_, ThirdJet_chEmEn_, ThirdJet_NHadEnFraction_,
         ThirdJet_chMuEn_,         ThirdJet_chMuEnFraction_, ThirdJet_numberOfDaughters_, 
         ThirdJet_muonEnergyFraction_, ThirdJet_muonMultiplicity_, 
         ThirdJet_neutralEmEnergy_,  ThirdJet_neutralHadronEnergy_, ThirdJet_NHadMultiplicity_, 
         ThirdJet_NMultiplicity_, ThirdJet_chHadEn_, ThirdJet_muonEnergy_, ThirdJet_chargedMultiplicity_ ,
         ThirdJet_chargedEmEnergyFraction, ThirdJet_chargedHadronEnergyFraction;

 Float_t ThirdJet_CsvV2, ThirdJet_DeepCsv_udsg, ThirdJet_DeepCsv_b, ThirdJet_DeepCsv_c, ThirdJet_DeepCsv_bb, ThirdJet_HadronFlavor;


 TBranch* branch_ThirdJet_charge_  = newtree->Branch("ThirdJet_charge_", &ThirdJet_charge_, "ThirdJet_charge_/F");
 TBranch* branch_ThirdJet_et_  = newtree->Branch("ThirdJet_et_", &ThirdJet_et_, "ThirdJet_et_/F");
 TBranch* branch_ThirdJet_pt_  = newtree->Branch("ThirdJet_pt_", &ThirdJet_pt_, "ThirdJet_pt_/F");
 TBranch* branch_ThirdJet_eta_  = newtree->Branch("ThirdJet_eta_", &ThirdJet_eta_, "ThirdJet_eta_/F");
 TBranch* branch_ThirdJet_phi_  = newtree->Branch("ThirdJet_phi_", &ThirdJet_phi_, "ThirdJet_phi_/F");
 TBranch* branch_ThirdJet_theta_  = newtree->Branch("ThirdJet_theta_", &ThirdJet_theta_, "ThirdJet_theta_/F");
 TBranch* branch_ThirdJet_en_  = newtree->Branch("ThirdJet_en_", &ThirdJet_en_, "ThirdJet_en_/F");
 TBranch* branch_ThirdJet_chargedEmEnergy_  = newtree->Branch("ThirdJet_chEmEn_", &ThirdJet_chEmEn_, "ThirdJet_chEmEn_/F");
 TBranch* branch_ThirdJet_NEmEnFraction_  = newtree->Branch("ThirdJet_NEmEnFraction_", &ThirdJet_NEmEnFraction_, "ThirdJet_NEmEnFraction_/F");
 TBranch* branch_ThirdJet_chHadEn_  = newtree->Branch("ThirdJet_chHadEn_", &ThirdJet_chHadEn_, "ThirdJet_chHadEn_/F");
 TBranch* branch_ThirdJet_NHadEnFraction_  = newtree->Branch("ThirdJet_NHadEnFraction_", &ThirdJet_NHadEnFraction_, "ThirdJet_NHadEnFraction_/F");
 TBranch* branch_ThirdJet_chMuEn_  = newtree->Branch("ThirdJet_chMuEn_", &ThirdJet_chMuEn_, "ThirdJet_chMuEn_/F");
 TBranch* branch_ThirdJet_chMuEnFraction_  = newtree->Branch("ThirdJet_chMuEnFraction_", &ThirdJet_chMuEnFraction_, "ThirdJet_chMuEnFraction_/F");
 TBranch* branch_ThirdJet_numberOfDaughters_  = newtree->Branch("ThirdJet_numberOfDaughters_", &ThirdJet_numberOfDaughters_, "ThirdJet_numberOfDaughters_/F");
 TBranch* branch_ThirdJet_muonEnergy_  = newtree->Branch("ThirdJet_muonEnergy_", &ThirdJet_muonEnergy_, "ThirdJet_muonEnergy_/F");
 TBranch* branch_ThirdJet_muonEnergyFraction_  = newtree->Branch("ThirdJet_muonEnergyFraction_", &ThirdJet_muonEnergyFraction_, "ThirdJet_muonEnergyFraction_/F");
 TBranch* branch_ThirdJet_muonMultiplicity_  = newtree->Branch("ThirdJet_muonMultiplicity_", &ThirdJet_muonMultiplicity_, "ThirdJet_muonMultiplicity_/F");
 TBranch* branch_ThirdJet_neutralEmEnergy_  = newtree->Branch("ThirdJet_neutralEmEnergy_", &ThirdJet_neutralEmEnergy_, "ThirdJet_neutralEmEnergy_/F");
 TBranch* branch_ThirdJet_neutralHadronEnergy_  = newtree->Branch("ThirdJet_neutralHadronEnergy_", &ThirdJet_neutralHadronEnergy_, "ThirdJet_neutralHadronEnergy_/F");
 TBranch* branch_ThirdJet_NHadMultiplicity_  = newtree->Branch("ThirdJet_NHadMultiplicity_", &ThirdJet_NHadMultiplicity_, "ThirdJet_NHadMultiplicity_/F");
 TBranch* branch_ThirdJet_NMultiplicity_  = newtree->Branch("ThirdJet_NMultiplicity_", &ThirdJet_NMultiplicity_, "ThirdJet_NMultiplicity_/F");
 TBranch* branch_ThirdJet_chargedMultiplicity_   = newtree->Branch("ThirdJet_chargedMultiplicity_", &ThirdJet_chargedMultiplicity_, "ThirdJet_chargedMultiplicity_/F");
 TBranch* branch_ThirdJet_chargedEmEnergyFraction =  newtree->Branch("ThirdJet_chargedEmEnergyFraction", &ThirdJet_chargedEmEnergyFraction, "ThirdJet_chargedEmEnergyFraction/F");
 TBranch* branch_ThirdJet_chargedHadronEnergyFraction = newtree->Branch("ThirdJet_chargedHadronEnergyFraction", &ThirdJet_chargedHadronEnergyFraction,"ThirdJet_chargedHadronEnergyFraction/F");
 
 TBranch* branch_ThirdJet_CsvV2        =  newtree->Branch("ThirdJet_CsvV2", &ThirdJet_CsvV2, "ThirdJet_CsvV2/F");
 TBranch* branch_ThirdJet_DeepCsv_udsg =  newtree->Branch("ThirdJet_DeepCsv_udsg", &ThirdJet_DeepCsv_udsg, "ThirdJet_DeepCsv_udsg/F");
 TBranch* branch_ThirdJet_DeepCsv_b    =  newtree->Branch("ThirdJet_DeepCsv_b", &ThirdJet_DeepCsv_b, "ThirdJet_DeepCsv_b/F");
 TBranch* branch_ThirdJet_DeepCsv_c    =  newtree->Branch("ThirdJet_DeepCsv_c", &ThirdJet_DeepCsv_c,"ThirdJet_DeepCsv_c/F");
 TBranch* branch_ThirdJet_DeepCsv_bb   =  newtree->Branch("ThirdJet_DeepCsv_bb", &ThirdJet_DeepCsv_bb,"ThirdJet_DeepCsv_bb/F");
 TBranch* branch_ThirdJet_HadronFlavor =  newtree->Branch("ThirdJet_HadronFlavor", &ThirdJet_HadronFlavor,"ThirdJet_HadronFlavor/F");
 
//============================ Missing Energy Variables =========================================//
 Float_t pfMet_et_, pfMet_pt_, pfMet_phi_, pfMet_en_, pfMet_sumEt_, caloMet_phi_,caloMet_pt_;

 TBranch* branch_pfMet_et_    = newtree->Branch("pfMet_et_", &pfMet_et_, "pfMet_et_/F");
 TBranch* branch_pfMet_pt_    = newtree->Branch("pfMet_pt_", &pfMet_pt_, "pfMet_pt_/F");
 TBranch* branch_pfMet_phi_   = newtree->Branch("pfMet_phi_", &pfMet_phi_, "pfMet_phi_/F");
 TBranch* branch_pfMet_en_    = newtree->Branch("pfMet_en_", &pfMet_en_, "pfMet_en_/F");
 TBranch* branch_pfMet_sumEt_ = newtree->Branch("pfMet_sumEt_", &pfMet_sumEt_, "pfMet_sumEt_/F");
 TBranch* branch_caloMet_pt_  = newtree->Branch("caloMet_pt_", &caloMet_pt_, "caloMet_pt_/F");
 TBranch* branch_caloMet_phi_ = newtree->Branch("caloMet_phi_", &caloMet_phi_, "caloMet_phi_/F");

 /*  
//=============================== first b-tagging Variables ===========================================//
 Float_t FirstJet_btag_pt_, FirstJet_btag_eta_,  FirstJet_btag_phi_, FirstJet_btag_discriminator_;
 Int_t   FirstJet_btag_flavor_, bJet_size;

 TBranch* branch_FirstJet_btag_pt_  = newtree->Branch("FirstJet_btag_pt_", &FirstJet_btag_pt_, "FirstJet_btag_pt_/F");
 TBranch* branch_FirstJet_btag_eta_  = newtree->Branch("FirstJet_btag_eta_", &FirstJet_btag_eta_, "FirstJet_btag_eta_/F");
 TBranch* branch_FirstJet_btag_phi_  = newtree->Branch("FirstJet_btag_phi_", &FirstJet_btag_phi_, "FirstJet_btag_phi_/F");
 TBranch* branch_FirstJet_btag_flavor_  = newtree->Branch("FirstJet_btag_flavor_", &FirstJet_btag_flavor_, "FirstJet_btag_flavor_/F");
 TBranch* branch_FirstJet_btag_discriminator_  = newtree->Branch("FirstJet_btag_discriminator_", &FirstJet_btag_discriminator_, "FirstJet_btag_discriminator_/F");

 TBranch* branch_bJet_size  = newtree->Branch("bJet_size", &bJet_size, "bJet_size/F");

 //=============================== second b-tagging Variables ===========================================// 
 Float_t SecondJet_btag_pt_, SecondJet_btag_eta_,  SecondJet_btag_phi_, SecondJet_btag_discriminator_;
 Int_t   SecondJet_btag_flavor_;

 TBranch* branch_SecondJet_btag_pt_  = newtree->Branch("SecondJet_btag_pt_", &SecondJet_btag_pt_, "SecondJet_btag_pt_/F");
 TBranch* branch_SecondJet_btag_eta_  = newtree->Branch("SecondJet_btag_eta_", &SecondJet_btag_eta_, "SecondJet_btag_eta_/F");
 TBranch* branch_SecondJet_btag_phi_  = newtree->Branch("SecondJet_btag_phi_", &SecondJet_btag_phi_, "SecondJet_btag_phi_/F");
 TBranch* branch_SecondJet_btag_flavor_  = newtree->Branch("SecondJet_btag_flavor_", &SecondJet_btag_flavor_, "SecondJet_btag_flavor_/F");
 TBranch* branch_SecondJet_btag_discriminator_  = newtree->Branch("SecondJet_btag_discriminator_", &SecondJet_btag_discriminator_, "SecondJet_btag_discriminator_/F");

 //=============================== third b-tagging Variables ===========================================//
 Float_t ThirdJet_btag_pt_, ThirdJet_btag_eta_,  ThirdJet_btag_phi_, ThirdJet_btag_discriminator_;
 Int_t   ThirdJet_btag_flavor_;

 TBranch* branch_ThirdJet_btag_pt_  = newtree->Branch("ThirdJet_btag_pt_", &ThirdJet_btag_pt_, "ThirdJet_btag_pt_/F");
 TBranch* branch_ThirdJet_btag_eta_  = newtree->Branch("ThirdJet_btag_eta_", &ThirdJet_btag_eta_, "ThirdJet_btag_eta_/F");
 TBranch* branch_ThirdJet_btag_phi_  = newtree->Branch("ThirdJet_btag_phi_", &ThirdJet_btag_phi_, "ThirdJet_btag_phi_/F");
 TBranch* branch_ThirdJet_btag_flavor_  = newtree->Branch("ThirdJet_btag_flavor_", &ThirdJet_btag_flavor_, "ThirdJet_btag_flavor_/F");
 TBranch* branch_ThirdJet_btag_discriminator_  = newtree->Branch("ThirdJet_btag_discriminator_", &ThirdJet_btag_discriminator_, "ThirdJet_btag_discriminator_/F");
*/
//======================= Start the running over input branches ==========================================//
 for (int i=0;i<oldtree->GetEntriesFast(); i++) {
    if (i%10000==0) cout<<i<<endl;
    oldtree->GetEntry(i);

    unsigned npt_ = -1;
    for(unsigned i=0; i<npT->size(); ++i){
      h_nTrueInteractions50->Fill(npT->at(i));
      h_nTrueInteractions100->Fill(npT->at(i));
      npt_ = i;
    }

    if (passIsoMu24All==0 && passIsoMu27All == 0) continue;  // cut on the trigger!

    unsigned pu = -1;
    for(unsigned i=0; i<PU_Weight->size(); ++i){
      pu = i;
    }

    unsigned pu_Up = -1;
    for(unsigned i=0; i<PU_WeightUp->size(); ++i){
      pu_Up = i;
    }
    unsigned pu_Down = -1;
    for(unsigned i=0; i<PU_WeightDown->size(); ++i){
      pu_Down = i;
    }

    Float_t   minPt_prompt = -1000;
    unsigned  FirstMuon = -1;
    for(unsigned i=0; i<mu_isTightMuon->size(); i++){
      if (mu_isTightMuon->at(i)==0.  || deltaBeta->at(i)>isoCut
	  || mu_ptTunePMuonBestTrack->at(i) < 25 || 
	  abs(mu_etaTunePMuonBestTrack->at(i)) > 2.4 ||
	  mu_absdxyTunePMuonBestTrack->at(i) > 0.02 ) continue;  
      if (mu_ptTunePMuonBestTrack->at(i) > minPt_prompt){
	minPt_prompt=mu_ptTunePMuonBestTrack->at(i);
	FirstMuon=i;	  
      }
    }
	
    //here I save only the info of the prompt muon -- for the sv muon we can clone the branch as it is..
    Float_t   minPt_second = -1000;
    unsigned  SecondMuon = -1;
    int count = 0;
    int nb_loose = 0;
    for(unsigned i=0; i<mu_isLoose->size(); i++){
      if(i == FirstMuon || FirstMuon == -1) continue;
      if(mu_isLoose->size() == 1) continue;
      if(mu_isLoose->at(i)  == 0 ) continue;
      nb_loose++;
      if(mu_ptTunePMuonBestTrack->at(i) < 5  
	 || abs(mu_etaTunePMuonBestTrack->at(i)) > 2.4 
	 || mu_ptTunePMuonBestTrack->at(i) > minPt_prompt 
	 || mu_absdxyTunePMuonBestTrack->at(i) < 0.02) continue;
      if(!(mu_isGlobalMuon->at(i) == 1 && 
	   mu_chi2LocalPositionMuonBestTrack->at(i) < 12 && 
	   mu_trkKinkMuonBestTrack->at(i) < 20 && 
	   mu_segmentCompatibilityMuonBestTrack->at(i) > 0.303) || 
	 !(mu_segmentCompatibilityMuonBestTrack->at(i) > 0.451)) continue;
    count++;
      if (mu_ptTunePMuonBestTrack->at(i)>minPt_second){
	minPt_second=mu_ptTunePMuonBestTrack->at(i);
	SecondMuon=i;
      }
    }    

    unsigned  SecondVertex = -1;
    Float_t   minPt_sv = -1000;
    if(SecondMuon != -1){
      for(unsigned i=0; i<sv_mu_pt->size(); i++){
	if( sv_mu_pt->at(i) < mu_ptTunePMuonBestTrack->at(SecondMuon)) continue;
	if(sv_mu_pt->at(i) > minPt_sv){
	  minPt_sv = sv_mu_pt->at(i);
	  SecondVertex=i;
	}
      }
    }

    // this the correct but didn't implented yet 
    unsigned  sv_FirstTrack = -1;
    Float_t  firstTrack_minPt = -1000;
    Float_t sumdxySig = 0;
    if(SecondVertex != -1 ){
      for(unsigned j=0; j<sv_mu_tracks_pt->at(SecondVertex).size(); j++){ 
	sumdxySig = +abs(sv_mu_tracks_dxySig->at(SecondVertex).at(j));
	if(sv_mu_tracks_pt->at(SecondVertex).at(j) > firstTrack_minPt){
	  firstTrack_minPt =sv_mu_tracks_pt->at(SecondVertex).at(j);
	  sv_FirstTrack = j;
	}
      }
    }


    // this the correct but didn't implented yet
    unsigned  sv_SecondTrack = -1;
    Float_t  secondTrack_minPt = -1000;
    if(SecondVertex != -1 ){
      for(unsigned j=0; j<sv_mu_tracks_pt->at(SecondVertex).size(); j++){
	if(j == sv_FirstTrack) continue;
	if(sv_mu_tracks_pt->at(SecondVertex).at(j) > secondTrack_minPt){
          secondTrack_minPt =sv_mu_tracks_pt->at(SecondVertex).at(j);
          sv_SecondTrack = j;
        }
      }
    }

    // first jet info
    unsigned  FirstJet = -1;    
    Float_t pt_1stjet = -1000;
    int jet_count = 0;
    int bjet1= 0;
    int bjet2= 0;
    int bjet3= 0;
    int bjet4= 0;
    int bjet5= 0;
    int bjet6= 0;
    int bjet7= 0;
    int bjet8= 0;
    int bjet9= 0;
    int bjet10= 0;
    int bjet11= 0;
    int bjet12= 0;

    for(unsigned i=0; i<jet_pt->size(); i++){
      if( (fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.99 
	   && jet_neutralEmEnergyFraction->at(i) <= 0.99 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
	  ||
	  ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0 
	    && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 ) 
	  ||
	  ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
	    && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2) 
	  ||	  
	  (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10)) {
	jet_count++;
	if(jet_pt->at(i) > pt_1stjet) {
	  pt_1stjet = jet_pt->at(i);
	  FirstJet = i;
	}

        if(jet_DeepCsv_b->at(i) > 0.1 ){ bjet1++; }
        if(jet_DeepCsv_b->at(i) > 0.2 ){ bjet2++; }
        if(jet_DeepCsv_b->at(i) > 0.3 ){ bjet3++; }
        if(jet_DeepCsv_b->at(i) > 0.4 ){ bjet4++; }
        if(jet_DeepCsv_b->at(i) > 0.5 ){ bjet5++; }
	if(jet_DeepCsv_b->at(i) > 0.6 ){ bjet6++; }
        if(jet_DeepCsv_b->at(i) > 0.7 ){ bjet7++; }
        if(jet_DeepCsv_b->at(i) > 0.8 ){ bjet8++; }
        if(jet_DeepCsv_b->at(i) > 0.9 ){ bjet9++; }
        if(jet_DeepCsv_b->at(i) > 0.2219 ){ bjet10++;}
        if(jet_DeepCsv_b->at(i) > 0.6324 ){ bjet11++;}
        if(jet_DeepCsv_b->at(i) > 0.8958 ){ bjet12++;}

      }
    }
   
    //second jet info
    unsigned  SecondJet = -1;
    Float_t pt_2ndjet = -1000;
    for(unsigned i=0; i<jet_pt->size(); i++){
      if(i == FirstJet) continue;
      if( (fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.99
           && jet_neutralEmEnergyFraction->at(i) <= 0.99 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
          ||
          ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0
            && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 )
          ||
          ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
            && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2)
	  ||
          (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10)) {
	if(jet_pt->at(i) > pt_2ndjet) {
	  pt_2ndjet = jet_pt->at(i);
	  SecondJet = i;
	}
      }
    }
    //third jet info
    unsigned  ThirdJet = -1;
    Float_t pt_3rdjet = -1000;
    for(unsigned i=0; i<jet_pt->size(); i++){
      if(i == FirstJet || i == SecondJet) continue;
      if( (fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.99
           && jet_neutralEmEnergyFraction->at(i) <= 0.99 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
          ||
          ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0
            && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 )
          ||
          ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
            && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2)
	  ||
          (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10)) {	
	if(jet_pt->at(i) > pt_3rdjet) {
	  pt_3rdjet = jet_pt->at(i);
	  ThirdJet = i;
	}
      }
    }
    if(FirstMuon != -1 && SecondMuon != -1 && SecondVertex != -1){

   trig_IsoMu24All = passIsoMu24All;
   trig_IsoMu27All = passIsoMu27All;
   trig_Mu17       = passMu17;
   trig_Mu17_Mu8_SameSign = passMu17_Mu8_SameSign;
   trig_Mu20       = passMu20;
   trig_Mu17_Mu8   = passMu17_Mu8; 
   trig_Mu27_TkMu8 = passMu27_TkMu8;

    //pile up weight
   pu_weight     = (isMC) ? PU_Weight->at(pu) : 0 ;
   pu_weightUp   = (isMC) ? PU_WeightUp->at(pu_Up) : 0 ;
   pu_weightDown = (isMC) ? PU_WeightDown->at(pu_Down) : 0 ;

   TrueNumInteractions = (isMC) ? npT->at(npt_) : 0;

   //pv info
   pvX_        =  pvX ;
   pvY_        =  pvY ;
   pvZ_        =  pvZ ;
   pvXErr_     =  pvXErr ;
   pvYErr_     =  pvYErr ;
   pvZErr_     =  pvZErr ;
   pvLxy_      =  pvLxy  ;
   pvLxyz_     =  pvLxyz ;
   pvLxySig_   =  pvLxySig  ;
   pvLxyzSig_  =  pvLxyzSig ;
   pvChi2_     =  pvChi2 ;
   pvSumPtSq_  =  pvSumPtSq ;
   numberPV_   =  numberPV;
   pvNTrack_   =  pvNTrack ;
   
    //prompt muon 
   mu_promptPt             = mu_ptTunePMuonBestTrack->at(FirstMuon);
   mu_promptEta            = mu_etaTunePMuonBestTrack->at(FirstMuon);
   mu_promptPhi            = mu_phiTunePMuonBestTrack->at(FirstMuon);
   mu_promptCharge         = mu_chargeTunePMuonBestTrack->at(FirstMuon);
   mu_promptE              = mu_en->at(FirstMuon);
   mu_promptEt             = mu_et->at(FirstMuon);
   mu_promptRhoIso         = mu_rhoIso->at(FirstMuon);
   mu_promptTrackiso       = mu_trackiso->at(FirstMuon);
   mu_promptPfSumChHadPt   = mu_pfSumChargedHadronPt->at(FirstMuon);   
   mu_promptPfSumNHadEt    = mu_pfSumNeutralHadronEt->at(FirstMuon);  
   mu_promptPFSumPhotonEt  = mu_PFSumPhotonEt->at(FirstMuon);  
   mu_promptPfSumPUPt      = mu_pfSumPUPt->at(FirstMuon);   
   mu_promptEmIso          = mu_emIso->at(FirstMuon);  
   mu_promptHadIso         = mu_hadIso->at(FirstMuon); 
   mu_promptNormalizedChi2 = mu_normalizedChi2->at(FirstMuon);
   mu_promptDPToverPT      = mu_dPToverPTTunePMuonBestTrack->at(FirstMuon); 
   mu_promptAbsdxy         = mu_absdxyTunePMuonBestTrack->at(FirstMuon);   
   mu_promptAbsdxyError    = mu_absdxyErrorTunePMuonBestTrack->at(FirstMuon);
   mu_promptAbsdxySig      = mu_absdxySigTunePMuonBestTrack->at(FirstMuon); 
   mu_promptAbsdz          = mu_absdzTunePMuonBestTrack->at(FirstMuon);  
   mu_promptAbsdz          = mu_absdzErrorTunePMuonBestTrack->at(FirstMuon);
   mu_promptAbsdzSig       = mu_absdzSigTunePMuonBestTrack->at(FirstMuon);  
   mu_promptRecoDeltaBeta  = deltaBeta->at(FirstMuon);  
   mu_promptRecoiso        = mu_recoiso->at(FirstMuon);  

   mu_promptDirection        = mu_STATofDirection->at(FirstMuon); 
   mu_promptNDof             = mu_STATofNDof->at(FirstMuon);  
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
   mu_promptPixelLayers      = mu_numberOfpixelLayersWithMeasurement->at(FirstMuon);
   mu_promptTrackerLayers    = mu_numberOftrackerLayersWithMeasurement->at(FirstMuon);

   mu_promptInnerTrackFraction   = mu_InnerTrackValidFraction->at(FirstMuon);
   mu_promptSegmentCompatibility = mu_segmentCompatibilityMuonBestTrack->at(FirstMuon);
   mu_promptTrkKink              = mu_trkKinkMuonBestTrack->at(FirstMuon);
   mu_promptChi2LocalPosition    = mu_chi2LocalPositionMuonBestTrack->at(FirstMuon);
   mu_promptGlobalMuon           = mu_isGlobalMuon->at(FirstMuon); 

   mu_promptRPCTofDirection        = mu_RPCTofDirection->at(FirstMuon);  
   mu_promptRPCTofNDof             = mu_RPCTofNDof->at(FirstMuon);
   mu_promptRPCTofTimeAtIpInOut    = mu_RPCTofTimeAtIpInOut->at(FirstMuon);
   mu_promptRPCTofTimeAtIpInOutErr = mu_RPCTofTimeAtIpInOutErr->at(FirstMuon);
   mu_promptRPCTofTimeAtIpOutIn    = mu_RPCTofTimeAtIpOutIn->at(FirstMuon);
   mu_promptRPCTofTimeAtIpOutInErr = mu_RPCTofTimeAtIpOutInErr->at(FirstMuon);

    //non_prompt muon
    mu_secondIsTight        = mu_isTightMuon->at(SecondMuon);
    mu_secondPt             = mu_ptTunePMuonBestTrack->at(SecondMuon);
    mu_secondEta            = mu_etaTunePMuonBestTrack->at(SecondMuon);
    mu_secondPhi            = mu_phiTunePMuonBestTrack->at(SecondMuon);
    mu_secondCharge         = mu_chargeTunePMuonBestTrack->at(SecondMuon);
    mu_secondE              = mu_en->at(SecondMuon);
    mu_secondEt             = mu_et->at(SecondMuon);
    mu_secondRhoIso         = mu_rhoIso->at(SecondMuon);
    mu_secondTrackiso       = mu_trackiso->at(SecondMuon);
    mu_secondPfSumChHadPt   = mu_pfSumChargedHadronPt->at(SecondMuon);   
    mu_secondPfSumNHadEt    = mu_pfSumNeutralHadronEt->at(SecondMuon);  
    mu_secondPFSumPhotonEt  = mu_PFSumPhotonEt->at(SecondMuon);  
    mu_secondPfSumPUPt      = mu_pfSumPUPt->at(SecondMuon);   
    mu_secondEmIso          = mu_emIso->at(SecondMuon);  
    mu_secondHadIso         = mu_hadIso->at(SecondMuon); 
    mu_secondNormalizedChi2 = mu_normalizedChi2->at(SecondMuon);
    mu_secondDPToverPT      = mu_dPToverPTTunePMuonBestTrack->at(SecondMuon); 
    mu_secondAbsdxy         = mu_absdxyTunePMuonBestTrack->at(SecondMuon);   
    mu_secondAbsdxyError    = mu_absdxyErrorTunePMuonBestTrack->at(SecondMuon);
    mu_secondAbsdxySig      = mu_absdxySigTunePMuonBestTrack->at(SecondMuon); 
    mu_secondAbsdz          = mu_absdzTunePMuonBestTrack->at(SecondMuon);  
    mu_secondAbsdzError     = mu_absdzErrorTunePMuonBestTrack->at(SecondMuon);
    mu_secondAbsdzSig       = mu_absdzSigTunePMuonBestTrack->at(SecondMuon);  
    mu_secondRecoDeltaBeta  = deltaBeta->at(SecondMuon);  
    mu_secondRecoiso        = mu_recoiso->at(SecondMuon);  

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
    mu_secondPixelLayers      = mu_numberOfpixelLayersWithMeasurement->at(SecondMuon);
    mu_secondTrackerLayers    = mu_numberOftrackerLayersWithMeasurement->at(SecondMuon);

    mu_secondInnerTrackFraction   = mu_InnerTrackValidFraction->at(FirstMuon);
    mu_secondSegmentCompatibility = mu_segmentCompatibilityMuonBestTrack->at(SecondMuon);
    mu_secondTrkKink              = mu_trkKinkMuonBestTrack->at(SecondMuon);
    mu_secondChi2LocalPosition    = mu_chi2LocalPositionMuonBestTrack->at(SecondMuon);
    mu_secondGlobalMuon           = mu_isGlobalMuon->at(SecondMuon);

    mu_secondRPCTofDirection        = mu_RPCTofDirection->at(FirstMuon);
    mu_secondRPCTofNDof             = mu_RPCTofNDof->at(FirstMuon);
    mu_secondRPCTofTimeAtIpInOut    = mu_RPCTofTimeAtIpInOut->at(FirstMuon);
    mu_secondRPCTofTimeAtIpInOutErr = mu_RPCTofTimeAtIpInOutErr->at(FirstMuon);
    mu_secondRPCTofTimeAtIpOutIn    = mu_RPCTofTimeAtIpOutIn->at(FirstMuon);
    mu_secondRPCTofTimeAtIpOutInErr = mu_RPCTofTimeAtIpOutInErr->at(FirstMuon);



    float R = sqrt(((mu_secondEta-mu_promptEta)*(mu_secondEta-mu_promptEta))+((mu_secondPhi-mu_promptPhi)*(mu_secondPhi-mu_promptPhi)));
    float charged =  mu_pfSumChargedHadronPt->at(SecondMuon);
    float neutral =  mu_pfSumNeutralHadronEt->at(SecondMuon);
    float sumPhotonEt = mu_PFSumPhotonEt->at(SecondMuon);
    float pileup =  mu_pfSumPUPt->at(SecondMuon); 

    double deltaBetaR3 = (charged + std::max(0.0, neutral+sumPhotonEt-0.5*pileup))/mu_secondPt;

    TLorentzVector Mu1;
    TLorentzVector Mu2; 
    TLorentzVector vtx;
    TLorentzVector Jet1;
    TLorentzVector Jet2;
    Mu1.SetPtEtaPhiE(mu_promptPt,mu_promptEta,mu_promptPhi,mu_promptE);
    Mu2.SetPtEtaPhiE(mu_secondPt,mu_secondEta,mu_secondPhi,mu_secondE);
    float DiMuMass = (Mu1 + Mu2).M();

    float px1 = mu_pxTunePMuonBestTrack->at(FirstMuon);
    float px2 = mu_pxTunePMuonBestTrack->at(SecondMuon);
    float py1 = mu_pyTunePMuonBestTrack->at(FirstMuon);
    float py2 = mu_pyTunePMuonBestTrack->at(SecondMuon);
    float mT  = sqrt(pow(mu_promptPt + mu_secondPt, 2) - pow(px1 + px2, 2) - pow(py1 + py2, 2));

    mu_DeltaBetaR3 = deltaBetaR3;
    mu_DiMuMass    = DiMuMass;
    mu_Size        = count;
    mu_DeltaR      = R ;
    mu_mT          = mT ;
    mu_nbLoose     = nb_loose;

    unsigned  mu_promptJet = -1;
    Float_t pt_j1 = -1000;
    for(unsigned i=0; i<jet_pt->size(); i++){
      float etaj1 = jet_eta->at(i);
      float phij1 = jet_phi->at(i);
      float R = sqrt(((mu_promptEta-etaj1)*(mu_promptEta-etaj1))+((mu_promptPhi-phij1)*(mu_promptPhi-phij1)));
      if( (fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.99
           && jet_neutralEmEnergyFraction->at(i) <= 0.99 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
          ||
          ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0
            && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 )
          ||
          ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
            && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2)
          ||
          (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10)) {
	if(R < 0.5 && jet_pt->at(i) > pt_j1){
	  pt_j1 = jet_pt->at(i);
	  mu_promptJet = i;
	}
      }
    }

    promptJet_charge_ = (mu_promptJet != -1) ? jet_charge->at(mu_promptJet) : -999;
    promptJet_et_     = (mu_promptJet != -1) ? jet_et->at(mu_promptJet)     : -999;
    promptJet_pt_     = (mu_promptJet != -1) ? jet_pt->at(mu_promptJet)     : -999;
    promptJet_eta_    = (mu_promptJet != -1) ? jet_eta->at(mu_promptJet)    : -999; 
    promptJet_phi_    = (mu_promptJet != -1) ? jet_phi->at(mu_promptJet)    : -999;
    promptJet_theta_  = (mu_promptJet != -1) ? jet_theta->at(mu_promptJet)  : -999;
    promptJet_en_     = (mu_promptJet != -1) ? jet_en->at(mu_promptJet)     : -999;
    promptJet_chEmEn_ = (mu_promptJet != -1) ? jet_chargedEmEnergy->at(mu_promptJet)                      : -999; 
    promptJet_NEmEnFraction_  = (mu_promptJet != -1) ? jet_neutralEmEnergyFraction->at(mu_promptJet)      : -999;
    promptJet_chHadEn_        = (mu_promptJet != -1) ? jet_chargedHadronEnergy->at(mu_promptJet)          : -999;
    promptJet_NHadEnFraction_ = (mu_promptJet != -1) ? jet_neutralHadronEnergyFraction->at(mu_promptJet)  : -999;
    promptJet_chMuEn_         = (mu_promptJet != -1) ? jet_chargedMuEnergy->at(mu_promptJet)              : -999;
    promptJet_chMuEnFraction_      = (mu_promptJet != -1) ? jet_chargedMuEnergyFraction->at(mu_promptJet) : -999;
    promptJet_numberOfDaughters_   = (mu_promptJet != -1) ? jet_numberOfDaughters->at(mu_promptJet)       : -999;
    promptJet_muonEnergy_          = (mu_promptJet != -1) ? jet_muonEnergy->at(mu_promptJet)              : -999;
    promptJet_muonEnergyFraction_  = (mu_promptJet != -1) ? jet_muonEnergyFraction->at(mu_promptJet)      : -999;
    promptJet_muonMultiplicity_    = (mu_promptJet != -1) ? jet_muonMultiplicity->at(mu_promptJet)        : -999;
    promptJet_neutralEmEnergy_     = (mu_promptJet != -1) ? jet_neutralEmEnergy->at(mu_promptJet)         : -999;
    promptJet_neutralHadronEnergy_ = (mu_promptJet != -1) ? jet_neutralHadronEnergy->at(mu_promptJet)     : -999; 
    promptJet_NHadMultiplicity_    = (mu_promptJet != -1) ? jet_neutralHadronMultiplicity->at(mu_promptJet)  : -999;
    promptJet_NMultiplicity_       = (mu_promptJet != -1) ? jet_neutralMultiplicity->at(mu_promptJet)        : -999;
    promptJet_chargedMultiplicity_        = (mu_promptJet != -1) ? jet_chargedMultiplicity->at(mu_promptJet) : -999;
    promptJet_chargedEmEnergyFraction     = (mu_promptJet != -1) ? jet_chargedEmEnergyFraction->at(mu_promptJet)     : -999;
    promptJet_chargedHadronEnergyFraction = (mu_promptJet != -1) ? jet_chargedHadronEnergyFraction->at(mu_promptJet) : -999;

    promptJet_CsvV2          = (mu_promptJet != -1) ? jet_CsvV2->at(mu_promptJet)        : -999;
    promptJet_DeepCsv_udsg   = (mu_promptJet != -1) ? jet_DeepCsv_udsg->at(mu_promptJet) : -999;
    promptJet_DeepCsv_b      = (mu_promptJet != -1) ? jet_DeepCsv_b->at(mu_promptJet)    : -999;
    promptJet_DeepCsv_c      = (mu_promptJet != -1) ? jet_DeepCsv_c->at(mu_promptJet)    : -999;
    promptJet_DeepCsv_bb     = (mu_promptJet != -1) ? jet_DeepCsv_bb->at(mu_promptJet)   : -999;
    promptJet_HadronFlavor   = (mu_promptJet != -1) ? jet_HadronFlavor->at(mu_promptJet) : -999;

    unsigned  mu_secondJet = -1;
    Float_t pt_j2 = -1000;
    for(unsigned i=0; i<jet_pt->size(); i++){
      float etaj2 = jet_eta->at(i);
      float phij2 = jet_phi->at(i);
      float R = sqrt(((mu_secondEta-etaj2)*(mu_secondEta-etaj2))+((mu_secondPhi-phij2)*(mu_secondPhi-phij2)));

      if( ((fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.99
	    && jet_neutralEmEnergyFraction->at(i) <= 0.99 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
	   ||
	   ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0
	     && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 )
	   ||
	   ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
	     && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2)
	   ||
	   (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10) ) && i != mu_promptJet) {
	if(R < 0.8 && jet_pt->at(i) > pt_j2){
	  pt_j2 = jet_pt->at(i);
	  mu_secondJet = i;
	}
      }
    }

    non_promptJet_charge_ = (mu_secondJet != -1) ? jet_charge->at(mu_secondJet) : -999;
    non_promptJet_et_     = (mu_secondJet != -1) ? jet_et->at(mu_secondJet)     : -999;
    non_promptJet_pt_     = (mu_secondJet != -1) ? jet_pt->at(mu_secondJet)     : -999;
    non_promptJet_eta_    = (mu_secondJet != -1) ? jet_eta->at(mu_secondJet)    : -999; 
    non_promptJet_phi_    = (mu_secondJet != -1) ? jet_phi->at(mu_secondJet)    : -999;
    non_promptJet_theta_  = (mu_secondJet != -1) ? jet_theta->at(mu_secondJet)  : -999;
    non_promptJet_en_     = (mu_secondJet != -1) ? jet_en->at(mu_secondJet)     : -999;
    non_promptJet_chEmEn_ = (mu_secondJet != -1) ? jet_chargedEmEnergy->at(mu_secondJet)                      : -999; 
    non_promptJet_NEmEnFraction_  = (mu_secondJet != -1) ? jet_neutralEmEnergyFraction->at(mu_secondJet)      : -999;
    non_promptJet_chHadEn_        = (mu_secondJet != -1) ? jet_chargedHadronEnergy->at(mu_secondJet)          : -999;
    non_promptJet_NHadEnFraction_ = (mu_secondJet != -1) ? jet_neutralHadronEnergyFraction->at(mu_secondJet)  : -999;
    non_promptJet_chMuEn_         = (mu_secondJet != -1) ? jet_chargedMuEnergy->at(mu_secondJet)              : -999;
    non_promptJet_chMuEnFraction_      = (mu_secondJet != -1) ? jet_chargedMuEnergyFraction->at(mu_secondJet) : -999;
    non_promptJet_numberOfDaughters_   = (mu_secondJet != -1) ? jet_numberOfDaughters->at(mu_secondJet)       : -999;
    non_promptJet_muonEnergy_          = (mu_secondJet != -1) ? jet_muonEnergy->at(mu_secondJet)              : -999;
    non_promptJet_muonEnergyFraction_  = (mu_secondJet != -1) ? jet_muonEnergyFraction->at(mu_secondJet)      : -999;
    non_promptJet_muonMultiplicity_    = (mu_secondJet != -1) ? jet_muonMultiplicity->at(mu_secondJet)        : -999;
    non_promptJet_neutralEmEnergy_     = (mu_secondJet != -1) ? jet_neutralEmEnergy->at(mu_secondJet)         : -999;
    non_promptJet_neutralHadronEnergy_ = (mu_secondJet != -1) ? jet_neutralHadronEnergy->at(mu_secondJet)     : -999; 
    non_promptJet_NHadMultiplicity_    = (mu_secondJet != -1) ? jet_neutralHadronMultiplicity->at(mu_secondJet)  : -999;
    non_promptJet_NMultiplicity_       = (mu_secondJet != -1) ? jet_neutralMultiplicity->at(mu_secondJet)        : -999;
    non_promptJet_chargedMultiplicity_        = (mu_secondJet != -1) ? jet_chargedMultiplicity->at(mu_secondJet) : -999;
    non_promptJet_chargedEmEnergyFraction     = (mu_secondJet != -1) ? jet_chargedEmEnergyFraction->at(mu_secondJet)     : -999;
    non_promptJet_chargedHadronEnergyFraction = (mu_secondJet != -1) ? jet_chargedHadronEnergyFraction->at(mu_secondJet) : -999;

    non_promptJet_CsvV2          = (mu_secondJet != -1) ? jet_CsvV2->at(mu_secondJet)        : -999;
    non_promptJet_DeepCsv_udsg   = (mu_secondJet != -1) ? jet_DeepCsv_udsg->at(mu_secondJet) : -999;
    non_promptJet_DeepCsv_b      = (mu_secondJet != -1) ? jet_DeepCsv_b->at(mu_secondJet)    : -999;
    non_promptJet_DeepCsv_c      = (mu_secondJet != -1) ? jet_DeepCsv_c->at(mu_secondJet)    : -999;
    non_promptJet_DeepCsv_bb     = (mu_secondJet != -1) ? jet_DeepCsv_bb->at(mu_secondJet)   : -999;
    non_promptJet_HadronFlavor   = (mu_secondJet != -1) ? jet_HadronFlavor->at(mu_secondJet) : -999;


    Jet1.SetPtEtaPhiE(promptJet_pt_,promptJet_eta_,promptJet_phi_,promptJet_en_);
    Jet2.SetPtEtaPhiE(non_promptJet_pt_,non_promptJet_eta_,non_promptJet_phi_,non_promptJet_en_);

    float DiJetMass = (mu_secondJet != -1 && mu_promptJet != -1) ? (Jet1 + Jet2).M()  : -999;
    float MuJetMass = (mu_secondJet != -1) ? (Mu1 + Jet2).M() : -999 ;

    DiJet_Mass = DiJetMass;
    MuJet_Mass = MuJetMass;

    //secondary vertex info   
    sv_inside     =  (SecondVertex != -1) ? SecondVertex                      : -999;
    sv_TrackSize  =  (SecondVertex != -1) ? sv_mu_TrackSize->at(SecondVertex) : -999;
    sv_LXYSig     =  (SecondVertex != -1) ? sv_mu_LXYSig->at(SecondVertex)    : -999;
    sv_LXYZSig    =  (SecondVertex != -1) ? sv_mu_LXYZSig->at(SecondVertex)   : -999;
    sv_LXY        =  (SecondVertex != -1) ? sv_mu_LXY->at(SecondVertex)       : -999;
    sv_LXYZ       =  (SecondVertex != -1) ? sv_mu_LXYZ->at(SecondVertex)      : -999;
    sv_dir_x      =  (SecondVertex != -1) ? sv_mu_dir_x->at(SecondVertex)     : -999;
    sv_dir_y      =  (SecondVertex != -1) ? sv_mu_dir_y->at(SecondVertex)     : -999;
    sv_dir_z      =  (SecondVertex != -1) ? sv_mu_dir_z->at(SecondVertex)     : -999;
    sv_mass       =  (SecondVertex != -1) ? sv_mu_mass->at(SecondVertex)      : -999;
    sv_eta        =  (SecondVertex != -1) ? sv_mu_eta->at(SecondVertex)       : -999;
    sv_phi        =  (SecondVertex != -1) ? sv_mu_phi->at(SecondVertex)       : -999;
    sv_pt         =  (SecondVertex != -1) ? sv_mu_pt->at(SecondVertex)        : -999;
    sv_p          =  (SecondVertex != -1) ? sv_mu_p->at(SecondVertex)         : -999;
    sv_px         =  (SecondVertex != -1) ? sv_mu_px->at(SecondVertex)        : -999;
    sv_py         =  (SecondVertex != -1) ? sv_mu_py->at(SecondVertex)        : -999;
    sv_pz         =  (SecondVertex != -1) ? sv_mu_pz->at(SecondVertex)        : -999;
    sv_energy     =  (SecondVertex != -1) ? sv_mu_energy->at(SecondVertex)    : -999;
    sv_Beta       =  (SecondVertex != -1) ? sv_mu_Beta->at(SecondVertex)      : -999;
    sv_Gamma      =  (SecondVertex != -1) ? sv_mu_Gamma->at(SecondVertex)     : -999;
    sv_CTau0      =  (SecondVertex != -1) ? sv_mu_CTau0->at(SecondVertex)     : -999;
    sv_NDof       =  (SecondVertex != -1) ? sv_mu_NDof->at(SecondVertex)      : -999;
    sv_Chi2       =  (SecondVertex != -1) ? sv_mu_Chi2->at(SecondVertex)      : -999;
    sv_Angle3D    =  (SecondVertex != -1) ?   sv_mu_Angle3D->at(SecondVertex) : -999;
    sv_Angle2D    =  (SecondVertex != -1) ? sv_mu_Angle2D->at(SecondVertex)   : -999;
    sv_tracks_Sumcharge =   (SecondVertex != -1) ? sv_mu_tracks_Sumcharge->at(SecondVertex): -999;
    sv_tracks_Sumpt     =   (SecondVertex != -1) ? sv_mu_tracks_Sumpt->at(SecondVertex)    : -999;
    sv_match      =  (SecondVertex != -1) ? sv_mu_match->at(SecondVertex)     : -999;
    sv_track_sumdxySig  = (SecondVertex != -1) ? sumdxySig  : -999;

    vtx.SetPtEtaPhiE(sv_pt,sv_eta,sv_phi,sv_energy);

    float vtxmumass =  (SecondVertex != -1) ? (Mu1 + vtx).M() : -999 ;

    vtxmu_mass = vtxmumass;

    sv_Xpos     = (SecondVertex != -1) ? sv_mu_Xpos->at(SecondVertex)  : -999;
    sv_Ypos     = (SecondVertex != -1) ? sv_mu_Ypos->at(SecondVertex)  : -999;
    sv_Zpos     = (SecondVertex != -1) ? sv_mu_Zpos->at(SecondVertex)  : -999;
    sv_xError   = (SecondVertex != -1) ? sv_mu_xError->at(SecondVertex)  : -999;
    sv_yError   = (SecondVertex != -1) ? sv_mu_yError->at(SecondVertex)  : -999;
    sv_zError   = (SecondVertex != -1) ? sv_mu_zError->at(SecondVertex)  : -999;
    sv_pvX      = (SecondVertex != -1) ? sv_mu_pvX->at(SecondVertex)  : -999;
    sv_pvY      = (SecondVertex != -1) ? sv_mu_pvY->at(SecondVertex)  : -999;
    sv_pvZ      = (SecondVertex != -1) ? sv_mu_pvZ->at(SecondVertex)  : -999;
    sv_pvXError = (SecondVertex != -1) ? sv_mu_pvXError->at(SecondVertex)  : -999;
    sv_pvYError = (SecondVertex != -1) ? sv_mu_pvYError->at(SecondVertex)  : -999;
    sv_pvZError = (SecondVertex != -1) ? sv_mu_pvZError->at(SecondVertex)  : -999;


    //first track in sv info
    firstTrack_eta    =  (SecondVertex != -1) ? sv_mu_tracks_eta->at(SecondVertex).at(sv_FirstTrack)    : -999; 
    firstTrack_phi    =  (SecondVertex != -1) ? sv_mu_tracks_phi->at(SecondVertex).at(sv_FirstTrack)    : -999;
    firstTrack_pt     =  (SecondVertex != -1) ? sv_mu_tracks_pt->at(SecondVertex).at(sv_FirstTrack)     : -999;
    firstTrack_dxySig =  (SecondVertex != -1) ? sv_mu_tracks_dxySig->at(SecondVertex).at(sv_FirstTrack) : -999;
    firstTrack_dxy    =  (SecondVertex != -1) ? sv_mu_tracks_dxy->at(SecondVertex).at(sv_FirstTrack)    : -999;
    firstTrack_dxyz   =  (SecondVertex != -1) ? sv_mu_tracks_dxyz->at(SecondVertex).at(sv_FirstTrack)   : -999;
    firstTrack_charge =  (SecondVertex != -1) ? sv_mu_tracks_charge->at(SecondVertex).at(sv_FirstTrack) : -999; 
    firstTrack_en     =  (SecondVertex != -1) ? sv_mu_tracks_en->at(SecondVertex).at(sv_FirstTrack)     : -999;

    //second track in sv info
    secondTrack_eta    =  (SecondVertex != -1) ? sv_mu_tracks_eta->at(SecondVertex).at(sv_SecondTrack)   : -999;
    secondTrack_phi    =  (SecondVertex != -1) ? sv_mu_tracks_phi->at(SecondVertex).at(sv_SecondTrack)   : -999;
    secondTrack_pt     =  (SecondVertex != -1) ? sv_mu_tracks_pt->at(SecondVertex).at(sv_SecondTrack)    : -999;
    secondTrack_dxySig =  (SecondVertex != -1) ? sv_mu_tracks_dxySig->at(SecondVertex).at(sv_SecondTrack): -999;
    secondTrack_dxy    =  (SecondVertex != -1) ? sv_mu_tracks_dxy->at(SecondVertex).at(sv_SecondTrack)   : -999;
    secondTrack_dxyz   =  (SecondVertex != -1) ? sv_mu_tracks_dxyz->at(SecondVertex).at(sv_SecondTrack)  : -999;
    secondTrack_charge =  (SecondVertex != -1) ? sv_mu_tracks_charge->at(SecondVertex).at(sv_SecondTrack): -999;
    secondTrack_en     =  (SecondVertex != -1) ? sv_mu_tracks_en->at(SecondVertex).at(sv_SecondTrack)    : -999;

    //leading jet info
    FirstJet_charge_ = (FirstJet != -1) ? jet_charge->at(FirstJet) : -999;
    FirstJet_et_     = (FirstJet != -1) ? jet_et->at(FirstJet)     : -999;
    FirstJet_pt_     = (FirstJet != -1) ? jet_pt->at(FirstJet)     : -999;
    FirstJet_eta_    = (FirstJet != -1) ? jet_eta->at(FirstJet)    : -999; 
    FirstJet_phi_    = (FirstJet != -1) ? jet_phi->at(FirstJet)    : -999;
    FirstJet_theta_  = (FirstJet != -1) ? jet_theta->at(FirstJet)  : -999;
    FirstJet_en_     = (FirstJet != -1) ? jet_en->at(FirstJet)     : -999;
    FirstJet_chEmEn_ = (FirstJet != -1) ? jet_chargedEmEnergy->at(FirstJet)                      : -999; 
    FirstJet_NEmEnFraction_  = (FirstJet != -1) ? jet_neutralEmEnergyFraction->at(FirstJet)      : -999;
    FirstJet_chHadEn_        = (FirstJet != -1) ? jet_chargedHadronEnergy->at(FirstJet)          : -999;
    FirstJet_NHadEnFraction_ = (FirstJet != -1) ? jet_neutralHadronEnergyFraction->at(FirstJet)  : -999;
    FirstJet_chMuEn_         = (FirstJet != -1) ? jet_chargedMuEnergy->at(FirstJet)              : -999;
    FirstJet_chMuEnFraction_      = (FirstJet != -1) ? jet_chargedMuEnergyFraction->at(FirstJet) : -999;
    FirstJet_numberOfDaughters_   = (FirstJet != -1) ? jet_numberOfDaughters->at(FirstJet)       : -999;
    FirstJet_muonEnergy_          = (FirstJet != -1) ? jet_muonEnergy->at(FirstJet)              : -999;
    FirstJet_muonEnergyFraction_  = (FirstJet != -1) ? jet_muonEnergyFraction->at(FirstJet)      : -999;
    FirstJet_muonMultiplicity_    = (FirstJet != -1) ? jet_muonMultiplicity->at(FirstJet)        : -999;
    FirstJet_neutralEmEnergy_     = (FirstJet != -1) ? jet_neutralEmEnergy->at(FirstJet)         : -999;
    FirstJet_neutralHadronEnergy_ = (FirstJet != -1) ? jet_neutralHadronEnergy->at(FirstJet)     : -999; 
    FirstJet_NHadMultiplicity_    = (FirstJet != -1) ? jet_neutralHadronMultiplicity->at(FirstJet)  : -999;
    FirstJet_NMultiplicity_       = (FirstJet != -1) ? jet_neutralMultiplicity->at(FirstJet)        : -999;
    FirstJet_chargedMultiplicity_        = (FirstJet != -1) ? jet_chargedMultiplicity->at(FirstJet) : -999;
    FirstJet_chargedEmEnergyFraction     = (FirstJet != -1) ? jet_chargedEmEnergyFraction->at(FirstJet)     : -999;
    FirstJet_chargedHadronEnergyFraction = (FirstJet != -1) ? jet_chargedHadronEnergyFraction->at(FirstJet) : -999;

    FirstJet_CsvV2          = (FirstJet != -1) ? jet_CsvV2->at(FirstJet)        : -999;
    FirstJet_DeepCsv_udsg   = (FirstJet != -1) ? jet_DeepCsv_udsg->at(FirstJet) : -999;
    FirstJet_DeepCsv_b      = (FirstJet != -1) ? jet_DeepCsv_b->at(FirstJet)    : -999;
    FirstJet_DeepCsv_c      = (FirstJet != -1) ? jet_DeepCsv_c->at(FirstJet)    : -999;
    FirstJet_DeepCsv_bb     = (FirstJet != -1) ? jet_DeepCsv_bb->at(FirstJet)   : -999;
    FirstJet_HadronFlavor   = (FirstJet != -1) ? jet_HadronFlavor->at(FirstJet) : -999;

    jets_size = jet_count;
    bjet1_size = bjet1;
    bjet2_size = bjet2;
    bjet3_size = bjet3;
    bjet4_size = bjet4;
    bjet5_size = bjet5;
    bjet6_size = bjet6;
    bjet7_size = bjet7;
    bjet8_size = bjet8;
    bjet9_size = bjet9;
    bjet_L = bjet10;
    bjet_M = bjet11;
    bjet_T = bjet12;


    //second Jet info
    SecondJet_charge_ = (SecondJet != -1) ? jet_charge->at(SecondJet) : -999;
    SecondJet_et_     = (SecondJet != -1) ? jet_et->at(SecondJet) : -999;
    SecondJet_pt_     = (SecondJet != -1) ? jet_pt->at(SecondJet) : -999;
    SecondJet_eta_    = (SecondJet != -1) ? jet_eta->at(SecondJet) : -999; 
    SecondJet_phi_    = (SecondJet != -1) ? jet_phi->at(SecondJet) : -999;
    SecondJet_theta_  = (SecondJet != -1) ? jet_theta->at(SecondJet) : -999;
    SecondJet_en_     = (SecondJet != -1) ? jet_en->at(SecondJet) : -999;
    SecondJet_chEmEn_ = (SecondJet != -1) ? jet_chargedEmEnergy->at(SecondJet) : -999; 
    SecondJet_NEmEnFraction_     = (SecondJet != -1) ? jet_neutralEmEnergyFraction->at(SecondJet) : -999;
    SecondJet_chHadEn_           = (SecondJet != -1) ? jet_chargedHadronEnergy->at(SecondJet) : -999;
    SecondJet_NHadEnFraction_    = (SecondJet != -1) ? jet_neutralHadronEnergyFraction->at(SecondJet) : -999;
    SecondJet_chMuEn_            = (SecondJet != -1) ? jet_chargedMuEnergy->at(SecondJet) : -999;
    SecondJet_chMuEnFraction_    = (SecondJet != -1) ? jet_chargedMuEnergyFraction->at(SecondJet) : -999;
    SecondJet_numberOfDaughters_ = (SecondJet != -1) ? jet_numberOfDaughters->at(SecondJet) : -999;
    SecondJet_muonEnergy_ = (SecondJet != -1) ? jet_muonEnergy->at(SecondJet) : -999;
    SecondJet_muonEnergyFraction_  = (SecondJet != -1) ? jet_muonEnergyFraction->at(SecondJet) : -999;
    SecondJet_muonMultiplicity_    = (SecondJet != -1) ? jet_muonMultiplicity->at(SecondJet) : -999;
    SecondJet_neutralEmEnergy_     = (SecondJet != -1) ? jet_neutralEmEnergy->at(SecondJet) : -999;
    SecondJet_neutralHadronEnergy_ = (SecondJet != -1) ? jet_neutralHadronEnergy->at(SecondJet) : -999; 
    SecondJet_NHadMultiplicity_    = (SecondJet != -1) ? jet_neutralHadronMultiplicity->at(SecondJet) : -999;
    SecondJet_NMultiplicity_       = (SecondJet != -1) ? jet_neutralMultiplicity->at(SecondJet) : -999;
    SecondJet_chargedMultiplicity_ = (SecondJet != -1) ? jet_chargedMultiplicity->at(SecondJet) : -999;
    SecondJet_chargedEmEnergyFraction     = (SecondJet != -1) ? jet_chargedEmEnergyFraction->at(SecondJet) : -999;
    SecondJet_chargedHadronEnergyFraction =  (SecondJet != -1) ? jet_chargedHadronEnergyFraction->at(SecondJet) : -999;

    SecondJet_CsvV2          = (SecondJet != -1) ? jet_CsvV2->at(SecondJet)        : -999;
    SecondJet_DeepCsv_udsg   = (SecondJet != -1) ? jet_DeepCsv_udsg->at(SecondJet) : -999;
    SecondJet_DeepCsv_b      = (SecondJet != -1) ? jet_DeepCsv_b->at(SecondJet)    : -999;
    SecondJet_DeepCsv_c      = (SecondJet != -1) ? jet_DeepCsv_c->at(SecondJet)    : -999;
    SecondJet_DeepCsv_bb     = (SecondJet != -1) ? jet_DeepCsv_bb->at(SecondJet)   : -999;
    SecondJet_HadronFlavor   = (SecondJet != -1) ? jet_HadronFlavor->at(SecondJet) : -999;

    TLorentzVector FJet;
    TLorentzVector SJet;

    FJet.SetPtEtaPhiE(FirstJet_pt_,FirstJet_eta_,FirstJet_phi_,FirstJet_en_);
    SJet.SetPtEtaPhiE(SecondJet_pt_,SecondJet_eta_,SecondJet_phi_,SecondJet_en_);

    float R1 = sqrt(((FirstJet_eta_-SecondJet_eta_)*(FirstJet_eta_-SecondJet_eta_))+((FirstJet_phi_-SecondJet_phi_)*(FirstJet_phi_-SecondJet_phi_)));
    float R2 = sqrt(((mu_promptEta-FirstJet_eta_)*(mu_promptEta-FirstJet_eta_))+((mu_promptPhi-FirstJet_phi_)*(mu_promptPhi-FirstJet_phi_)));
    float R3 = sqrt(((mu_secondEta-FirstJet_eta_)*(mu_secondEta-FirstJet_eta_))+((mu_secondPhi-FirstJet_phi_)*(mu_secondPhi-FirstJet_phi_)));
    float R4 = sqrt(((mu_promptEta-SecondJet_eta_)*(mu_promptEta-SecondJet_eta_))+((mu_promptPhi-SecondJet_phi_)*(mu_promptPhi-SecondJet_phi_)));
    float R5 = sqrt(((mu_secondEta-SecondJet_eta_)*(mu_secondEta-SecondJet_eta_))+((mu_secondPhi-SecondJet_phi_)*(mu_secondPhi-SecondJet_phi_)));


    FSjet_DR   = (FirstJet != -1 && SecondJet != -1 ) ? R1 : -999;
    FjetFmu_DR = (FirstJet != -1  ) ? R2 : -999;
    FjetSmu_DR = (FirstJet != -1  ) ? R3 : -999;
    SjetFmu_DR = (SecondJet != -1) ? R4 : -999;
    SjetSmu_DR = (SecondJet != -1) ? R4 : -999;

    FSjet_M = (SecondJet != -1 && FirstJet != -1) ? (SJet + FJet).M() : -999 ;

    FjetFmu_M = (FirstJet != -1)  ? (Mu1 + FJet).M() : -999 ;
    FjetSmu_M = (FirstJet != -1)  ? (Mu2 + FJet).M() : -999 ;
    SjetFmu_M = (SecondJet != -1) ? (Mu1 + SJet).M() : -999 ;
    SjetSmu_M = (SecondJet != -1) ? (Mu2 + SJet).M() : -999 ;

    FSjetFmu_M = (SecondJet != -1 && FirstJet != -1) ? (SJet + FJet + Mu1).M() : -999 ;
    FSjetSmu_M = (SecondJet != -1 && FirstJet != -1) ? (SJet + FJet + Mu2).M() : -999 ;

    FSjetFSmu_M = (SecondJet != -1 && FirstJet != -1) ? (SJet + FJet + Mu1 + Mu2).M() : -999 ;


    //third Jet info
    ThirdJet_charge_ = (ThirdJet != -1) ? jet_charge->at(ThirdJet) : -999;
    ThirdJet_et_ = (ThirdJet != -1) ? jet_et->at(ThirdJet) : -999;
    ThirdJet_pt_ = (ThirdJet != -1) ? jet_pt->at(ThirdJet) : -999;
    ThirdJet_eta_ = (ThirdJet != -1) ? jet_eta->at(ThirdJet) : -999; 
    ThirdJet_phi_ = (ThirdJet != -1) ? jet_phi->at(ThirdJet) : -999;
    ThirdJet_theta_ = (ThirdJet != -1) ? jet_theta->at(ThirdJet) : -999;
    ThirdJet_en_ = (ThirdJet != -1) ? jet_en->at(ThirdJet) : -999;
    ThirdJet_chEmEn_ = (ThirdJet != -1) ? jet_chargedEmEnergy->at(ThirdJet) : -999; 
    ThirdJet_NEmEnFraction_ = (ThirdJet != -1) ? jet_neutralEmEnergyFraction->at(ThirdJet) : -999;
    ThirdJet_chHadEn_ = (ThirdJet != -1) ? jet_chargedHadronEnergy->at(ThirdJet) : -999;
    ThirdJet_NHadEnFraction_ = (ThirdJet != -1) ? jet_neutralHadronEnergyFraction->at(ThirdJet) : -999;
    ThirdJet_chMuEn_ = (ThirdJet != -1) ? jet_chargedMuEnergy->at(ThirdJet) : -999;
    ThirdJet_chMuEnFraction_ = (ThirdJet != -1) ? jet_chargedMuEnergyFraction->at(ThirdJet) : -999;
    ThirdJet_numberOfDaughters_ = (ThirdJet != -1) ? jet_numberOfDaughters->at(ThirdJet) : -999;
    ThirdJet_muonEnergy_ = (ThirdJet != -1) ? jet_muonEnergy->at(ThirdJet) : -999;
    ThirdJet_muonEnergyFraction_ = (ThirdJet != -1) ? jet_muonEnergyFraction->at(ThirdJet) : -999;
    ThirdJet_muonMultiplicity_ = (ThirdJet != -1) ? jet_muonMultiplicity->at(ThirdJet) : -999;
    ThirdJet_neutralEmEnergy_ = (ThirdJet != -1) ? jet_neutralEmEnergy->at(ThirdJet) : -999;
    ThirdJet_neutralHadronEnergy_ = (ThirdJet != -1) ? jet_neutralHadronEnergy->at(ThirdJet) : -999; 
    ThirdJet_NHadMultiplicity_ = (ThirdJet != -1) ? jet_neutralHadronMultiplicity->at(ThirdJet) : -999;
    ThirdJet_NMultiplicity_ = (ThirdJet != -1) ? jet_neutralMultiplicity->at(ThirdJet) : -999;
    ThirdJet_chargedMultiplicity_ = (ThirdJet != -1) ? jet_chargedMultiplicity->at(ThirdJet) : -999;
    ThirdJet_chargedEmEnergyFraction = (ThirdJet != -1) ? jet_chargedEmEnergyFraction->at(ThirdJet) : -999;
    ThirdJet_chargedHadronEnergyFraction =  (ThirdJet != -1) ? jet_chargedHadronEnergyFraction->at(ThirdJet) : -999;

    ThirdJet_CsvV2          = (ThirdJet != -1) ? jet_CsvV2->at(ThirdJet)        : -999;
    ThirdJet_DeepCsv_udsg   = (ThirdJet != -1) ? jet_DeepCsv_udsg->at(ThirdJet) : -999;
    ThirdJet_DeepCsv_b      = (ThirdJet != -1) ? jet_DeepCsv_b->at(ThirdJet)    : -999;
    ThirdJet_DeepCsv_c      = (ThirdJet != -1) ? jet_DeepCsv_c->at(ThirdJet)    : -999;
    ThirdJet_DeepCsv_bb     = (ThirdJet != -1) ? jet_DeepCsv_bb->at(ThirdJet)   : -999;
    ThirdJet_HadronFlavor   = (ThirdJet != -1) ? jet_HadronFlavor->at(ThirdJet) : -999;

    //missing energy Info
    pfMet_et_    =  pfMet_et;
    pfMet_pt_    =  pfMet_pt;
    pfMet_phi_   =  pfMet_phi;
    pfMet_en_    =  pfMet_en;
    pfMet_sumEt_ =  pfMet_sumEt;
    caloMet_pt_  =  caloMet_pt;
    caloMet_phi_ =  caloMet_phi;
    /*
    // first b info
    FirstJet_btag_pt_  = (first_bjet != -1 ) ? jet_btag_pt->at(first_bjet)   : -999;
    FirstJet_btag_eta_ = (first_bjet != -1 ) ? jet_btag_eta->at(first_bjet)  : -999;
    FirstJet_btag_phi_ = (first_bjet != -1 ) ? jet_btag_phi->at(first_bjet)  : -999;
    FirstJet_btag_flavor_ = (first_bjet != -1 ) ? jet_btag_flavor->at(first_bjet)  : -999;
    FirstJet_btag_discriminator_ = (first_bjet != -1 ) ? jet_btag_pfCSVv2IVF_discriminator->at(first_bjet)  : -999;

    bJet_size = (first_bjet != -1 ) ?  bjet_count : -999;

    //second b info
    SecondJet_btag_pt_  = (second_bjet != -1 ) ? jet_btag_pt->at(second_bjet)   : -999;
    SecondJet_btag_eta_ = (second_bjet != -1 ) ? jet_btag_eta->at(second_bjet)  : -999;
    SecondJet_btag_phi_ = (second_bjet != -1 ) ? jet_btag_phi->at(second_bjet)  : -999;
    SecondJet_btag_flavor_ = (second_bjet != -1 ) ? jet_btag_flavor->at(second_bjet)  : -999;
    SecondJet_btag_discriminator_ = (second_bjet != -1 ) ? jet_btag_pfCSVv2IVF_discriminator->at(second_bjet)  : -999;

    //third b info
    ThirdJet_btag_pt_  = (third_bjet != -1 ) ? jet_btag_pt->at(third_bjet)   : -999;
    ThirdJet_btag_eta_ = (third_bjet != -1 ) ? jet_btag_eta->at(third_bjet)  : -999;
    ThirdJet_btag_phi_ = (third_bjet != -1 ) ? jet_btag_phi->at(third_bjet)  : -999;
    ThirdJet_btag_flavor_ = (third_bjet != -1 ) ? jet_btag_flavor->at(third_bjet)  : -999;
    ThirdJet_btag_discriminator_ = (third_bjet != -1 ) ? jet_btag_pfCSVv2IVF_discriminator->at(third_bjet)  : -999;
    */
    newtree->Fill();
    }

  }
 h_nTrueInteractions50->Write();
 h_nTrueInteractions100->Write();

 newfile->Write();  

 newtree->Print();
 newtree->AutoSave();

  delete newfile;

  return 0;
}

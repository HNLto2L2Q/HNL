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

  std::string filepath = argv[1];
  std::string filename = argv[2];
  std::string   flagMC = argv[3];

  if( (flagMC.find("true") == string::npos && flagMC.find("True") == string::npos) &&
      (flagMC.find("false") == string::npos && flagMC.find("False") == string::npos) ){
    std::cout<<" Flag for specifing data or MC processing should be only [true/True] or [false/False]\n";
    return 0;
  }
  Float_t isoCut = 0.1;
  bool isMC = (flagMC.find("true") != string::npos || flagMC.find("True") != string::npos);
  TChain * oldtree = new TChain("HeavyNeutralLepton/tree_","");
  string tree ="HeavyNeutralLepton/tree_";
  oldtree->Add(Form("%s/%s/%s",filepath.c_str(),filename.c_str(),tree.c_str()));

  TH1F* h_nTrueInteractions50      = new TH1F("nTrueInteractions50" , "nTrueInteractions 50 bin" , 50, 0., 50. );
  TH1F* h_nTrueInteractions100     = new TH1F("nTrueInteractions100" , "nTrueInteractions 100 bin" , 100, 0., 100. );

  TFile *newfile = new TFile(filename.c_str(),"recreate");
  TTree *newtree  = new TTree("newtree","Analysis Tree");
  Long64_t nentries = oldtree->GetEntries();
  cout  <<nentries<<endl;

  /*
int main(){
  

  //Get old file, old tree and set top branch address
  //TFile *oldfile = new TFile("/eos/cms/store/group/phys_exotica/HNL/Background/crab_Analysis_WZToLLLNu/Background_Analysis.root");
  //TFile *oldfile = new TFile("/user/moanwar/Run2016/CMSSW_8_0_29/src/HNL/HeavyNeutralLeptonAnalysis/test/signal/HNL_M5_mu_2.95.root");
  //TTree *oldtree = (TTree*)oldfile->Get("HeavyNeutralLepton/tree_");


  //////selection cuts

  Float_t isoCut = 0.1;
  bool isMC = true;
  //if I want to use a TChain.....
  //cout<< "starting..."<<endl;
  TChain * oldtree = new TChain("HeavyNeutralLepton/tree_","");
  string tree ="HeavyNeutralLepton/tree_";
  //Background
  //string inputPath = "/pnfs/iihe/cms/store/user/moanwar/SamplesToSkimm_run5_muon";
  //string inputPath = "/user/moanwar/SamplesToSkimm_run6";
  //string fileName = "DYJetsToLL_M-10to50_madgraphMLM.root";
  //string fileName = "DYJetsToLL_M-50_madgraphMLM.root";
  //string fileName = "WJetsToLNu_TuneCUETP8M1_1.root";
  //string fileName = "WJetsToLNu_TuneCUETP8M1_2.root";
  //string fileName = "WJetsToLNu_TuneCUETP8M1_madgraph.root"; // the other version of wjets at madgraph *** recommeded start skimm from here
  //string fileName = "DYJetsToLL_M-10to50_TuneCUETP8M.root";
  //string fileName = "DYJetsToLL_M-50_TuneCUETP8M1.root"; 
  //string fileName = "TT_TuneCUETP8M2T4_13TeV_MuonChannel.root";
  //string fileName = "WJetsToLNu_TuneCUETP8M1.root"; 
  //string fileName = "GluGluHToZZTo4L_M125.root"; 
  //string fileName = "ST_s-channel_4f_leptonDecays.root"; 
  //string fileName = "ST_t-channel_antitop_4f_inclusiveDecays.root"; 
  //string fileName = "ST_t-channel_top_4f_inclusiveDecays.root"; 
  //string fileName = "ST_tW_antitop_5f_inclusiveDecays.root"; 
  //string fileName = "ST_tW_top_5f_inclusiveDecays.root"; 
  //string fileName = "TTGJets_TuneCUETP8M1.root"; 
  //string fileName = "ttWJets_13TeV_madgraphMLM.root"; 
  //string fileName = "TTWJetsToLNu_TuneCUETP8M1.root"; 
  //string fileName = "ttZJets_13TeV_madgraphMLM.root"; 
  //string fileName = "TTZToLLNuNu_M-10.root"; 
  //string fileName = "VHToNonbb_M125.root"; 
  //string fileName = "WGToLNuG_TuneCUETP8M1.root"; 
  //string fileName = "WpWpJJ_QCD_TuneCUETP8M1.root"; 
  //string fileName = "WWTo2L2Nu_13TeV.root"; 
  //string fileName = "WWToLNuQQ_13TeV-powheg.root"; 
  //string fileName = "WWW_4F_TuneCUETP8M1.root"; 
  //string fileName = "WWZ_TuneCUETP8M1.root"; 
  //string fileName = "WZTo1L3Nu_13TeV_amcatnloFXFX.root"; 
  //string fileName = "WZTo3LNu_TuneCUETP8M1.root"; 
  //string fileName = "WZToLNu2QorQQ2L_aTGC.root"; 
  //string fileName = "WZZ_TuneCUETP8M1.root"; 
  //string fileName = "ZGTo2LG_TuneCUETP8M1.root"; 
  //string fileName = "ZZTo2L2Nu_13TeV_powheg.root"; 
  //string fileName = "ZZTo2L2Q_13TeV_amcatnloFXFX.root"; 
  //string fileName = "ZZTo4L_13TeV-amcatnloFXFX.root"; 
  //string fileName = "ZZZ_TuneCUETP8M1.root"; 
  //string fileName = "QCD_Pt-120to170_Mu.root"; 
  //string fileName = "QCD_Pt-170to300_Mu.root"; 
  //string fileName = "QCD_Pt-20to30_Mu.root"; 
  //string fileName = "QCD_Pt-30to50_Mu.root"; 
  //string fileName = "QCD_Pt-50to80_Mu.root"; 
  //string fileName = "QCD_Pt-80to120_Mu.root"; 


  //signal
  string inputPath = "/user/moanwar/heavyneutrino/CMSSW_9_4_13/src/HNL/HeavyNeutralLeptonAnalysis/test/signal_samples_mu";
  //string fileName = "M_mu_M1_14.84.root";
  //string fileName = "M_mu_M2_158.00.root";  
  //string fileName = "M_mu_M2_38.87.root";   
  //string fileName = "M_mu_M2_64.78.root";   
  //string fileName = "M_mu_M3_22.82.root";   
  //string fileName = "M_mu_M3_32.60.root";   
  //string fileName = "M_mu_M3_45.64.root";   
  //string fileName = "M_mu_M1_74.22.root";   
  //string fileName = "M_mu_M2_31.50.root";   
  //string fileName = "M_mu_M2_41.00.root";   
  //string fileName = "M_mu_M2_97.17.root";   
  //string fileName = "M_mu_M3_23.15.root";   
  //string fileName = "M_mu_M3_45.55.root";   
  //string fileName = "M_mu_M3_76.07.root";   
  //string fileName = "M_mu_M5_14.77.root";
  //string fileName = "M_mu_M5_18.46.root";
  //string fileName = "M_mu_M5_24.61.root";
  //string fileName = "M_mu_M5_69.66.root";
  //string fileName = "M_mu_M5_92.88.root";
  //string fileName = "M_mu_M4_16.17.root";
  //string fileName = "M_mu_M4_24.25.root";
  //string fileName = "M_mu_M4_48.51.root";
  //string fileName = "M_mu_M4_57.47.root";
  //string fileName = "M_mu_M4_60.63.root";
  //string fileName = "M_mu_M4_76.87.root";
  //string fileName = "M_mu_M4_80.84.root";
  //string fileName = "M_mu_M8_1.25.root";
  //string fileName = "M_mu_M8_1.56.root";
  //string fileName = "M_mu_M8_1.89.root";
  //string fileName = "M_mu_M10_6.28.root";
  //string fileName = "M_mu_M10_6.90.root";
  //string fileName = "M_mu_M6_12.53.root";
  //string fileName = "M_mu_M6_13.72.root";
  //string fileName = "M_mu_M6_82.25.root";
  string fileName = "M_mu_M15.root";

  oldtree->Add(Form("%s/%s/%s",inputPath.c_str(),fileName.c_str(),tree.c_str()));


  TH1F* h_nTrueInteractions50      = new TH1F("nTrueInteractions50" , "nTrueInteractions 50 bin" , 50, 0., 50. );
  TH1F* h_nTrueInteractions100     = new TH1F("nTrueInteractions100" , "nTrueInteractions 100 bin" , 100, 0., 100. );

  //char fileName[256];
  //cout<<"Please enter the name of the output root file you want to create (yyy.root) : "<<endl;
  //cin.getline(fileName,256);

  string outputPath = " ~/run1_mu_MiniAODv3/";
  //string outputPath = "/user/moanwar/skimmed_test/";

  TFile *newfile = new TFile(Form("%s/%s",outputPath.c_str(),fileName.c_str()),"recreate");

  //TFile *newfile = new TFile(fileName,"recreate");

  //TFile *newfile = new TFile("skimmedSignale.root","recreate");
  //Create a new file + a clone of old tree in new file 
  //TTree *newtree = oldtree->CloneTree(0); 
  TTree *newtree  = new TTree("newtree","Analysis Tree");

  //cout<<"cloning done"<<endl;

  // Long64_t nentries = oldtree->GetEntries();
  Long64_t nentries = oldtree->GetEntries();
  cout  <<nentries<<endl;
  */
  //special case to be remove 


  /*

  Float_t lep2_gen_MomCTau0;
  oldtree->SetBranchAddress("lep2_gen_MomCTau0", &lep2_gen_MomCTau0);

  Float_t GenCtau ;

  TBranch* branch_GenCtau        = newtree->Branch("GenCtau", &GenCtau, "GenCtau/F");
  */

//======================= Old Tree Variables ==========================================// 
  // These are the variables I cut on 
  Bool_t passIsoMu24All;
  Bool_t passIsoMu27All;
  Bool_t passIsoMuTk24;
  Bool_t passIsoMu24;
  Bool_t passMu3_PFJet40;
  Bool_t passMu8_TrkIsoVVL;
  Bool_t passMu17_TrkIsoVVL;

  oldtree->SetBranchAddress("passIsoMu24" , &passIsoMu24);
  oldtree->SetBranchAddress("passIsoMuTk24" , &passIsoMuTk24);
  oldtree->SetBranchAddress("passIsoMu24All",&passIsoMu24All);
  oldtree->SetBranchAddress("passIsoMu27All",&passIsoMu27All);
  oldtree->SetBranchAddress("passMu3_PFJet40"   , &passMu3_PFJet40);
  oldtree->SetBranchAddress("passMu8_TrkIsoVVL" , &passMu8_TrkIsoVVL);
  oldtree->SetBranchAddress("passMu17_TrkIsoVVL", &passMu17_TrkIsoVVL);

  Float_t gen_weight ;
  Float_t lhe_weight ;
  Float_t lhe_ctau ;
  vector<Float_t>   *mu_DecayChain = 0;
  if(isMC){
    oldtree->SetBranchAddress("gen_weight" , &gen_weight);
    oldtree->SetBranchAddress("lhe_weight" , &lhe_weight);
    oldtree->SetBranchAddress("lhe_ctau"   , &lhe_ctau);
    oldtree->SetBranchAddress("mu_DecayChain"   , &mu_DecayChain);
  }


  vector<Float_t>   *npT = 0;
  oldtree->SetBranchAddress("npT",&npT);

  Float_t pvX ;
  Float_t pvY ;
  Float_t pvZ ;
  Float_t pvXErr ;
  Float_t pvYErr ;
  Float_t pvZErr ;
  Float_t pvLxy  ;
  Float_t pvLxyz ;
  Float_t pvLxySigma  ;
  Float_t pvLxyzSigma ;
  Float_t pvChi2 ;
  Float_t pvSumPtSq ;
  //int numberPV;
  //int pvNTrack ;

  oldtree->SetBranchAddress("pvX" , &pvX);
  oldtree->SetBranchAddress("pvY" , &pvY);
  oldtree->SetBranchAddress("pvZ" , &pvZ);
  oldtree->SetBranchAddress("pvXErr" , &pvXErr);
  oldtree->SetBranchAddress("pvYErr" , &pvYErr);
  oldtree->SetBranchAddress("pvZErr" , &pvZErr);
  oldtree->SetBranchAddress("pvLxy" , &pvLxy);
  oldtree->SetBranchAddress("pvLxyz" , &pvLxyz);
  oldtree->SetBranchAddress("pvLxySigma" , &pvLxySigma);
  oldtree->SetBranchAddress("pvLxyzSigma" , &pvLxyzSigma);
  oldtree->SetBranchAddress("pvChi2" , &pvChi2);
  //oldtree->SetBranchAddress("pvNTrack" , &pvNTrack);
  oldtree->SetBranchAddress("pvSumPtSq" , &pvSumPtSq);
  //oldtree->SetBranchAddress("numberPV" , &numberPV);


  vector<Float_t>   *mu_isTightMuon = 0;
  vector<Float_t>   *mu_isLooseMuon = 0;
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

  vector<Float_t>   *mu_3dIP = 0;
  vector<Float_t>   *mu_3dIPSig = 0 ;
  vector<Float_t>   *mu_2dIP = 0 ;
  vector<Float_t>   *mu_2dIPSig = 0 ;


   vector<Bool_t>  *sv_hasMuon = 0 ;
   vector<Int_t>   *sv_munTracks = 0 ;
   vector<Float_t> *sv_LxySig = 0 ;
   vector<Float_t> *sv_LxyzSig = 0 ;
   vector<Float_t> *sv_Lxy = 0 ;
   vector<Float_t> *sv_Lxyz = 0 ;
   vector<Float_t> *sv_mass = 0 ;
   vector<Float_t> *sv_eta = 0 ;
   vector<Float_t> *sv_phi = 0 ;
   vector<Float_t> *sv_pt = 0 ;
   vector<Float_t> *sv_p = 0 ;
   vector<Float_t> *sv_px = 0 ;
   vector<Float_t> *sv_py = 0 ;
   vector<Float_t> *sv_pz = 0 ;
   vector<Float_t> *sv_energy = 0 ;
   vector<Float_t> *sv_Beta = 0 ;
   vector<Float_t> *sv_Gamma = 0 ;
   vector<Float_t> *sv_CTau0 = 0 ;
   vector<Float_t> *sv_NDof = 0 ;
   vector<Float_t> *sv_Chi2 = 0 ;
   vector<Float_t> *sv_Angle3D = 0 ;
   vector<Float_t> *sv_Angle2D = 0 ;
   vector<Int_t>   *sv_tracks_Sumcharge = 0 ;
   vector<Float_t> *sv_tracks_Sumpt = 0 ;
   vector<Float_t> *sv_match_dxyz = 0 ;

   vector<float> *sv_X = 0;
   vector<float> *sv_Y = 0;
   vector<float> *sv_Z = 0;
   vector<float> *sv_xErr = 0;
   vector<float> *sv_yErr = 0;
   vector<float> *sv_zErr = 0;
   vector<float> *sv_lx = 0;
   vector<float> *sv_ly = 0;
   vector<float> *sv_lz = 0;
   // this wrong but leave it for now                                                                                                                                                                        
   vector<vector<int> >    *sv_tracks_charge = 0;
   vector<vector<float> >  *sv_tracks_eta = 0;
   vector<vector<float> >  *sv_tracks_phi = 0;
   vector<vector<float> >  *sv_tracks_pt  = 0;
   vector<vector<float> >  *sv_tracks_dxySig = 0;
   vector<vector<float> >  *sv_tracks_dxy = 0;
   vector<vector<float> >  *sv_tracks_dxyz = 0;
   vector<vector<float> >  *sv_tracks_p = 0;


   vector<Float_t>   *jetSmearedPt = 0;
   vector<Float_t>   *jet_ptuncorrected =0;
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
  Float_t   pfMet_px;
  Float_t   pfMet_py;
  Float_t   pfMet_pz;


  oldtree->SetBranchAddress("mu_isTightMuon",&mu_isTightMuon);
  oldtree->SetBranchAddress("mu_isLooseMuon",&mu_isLooseMuon);
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

  oldtree->SetBranchAddress("mu_3dIP",& mu_3dIP);
  oldtree->SetBranchAddress("mu_3dIPSig",& mu_3dIPSig);
  oldtree->SetBranchAddress("mu_2dIP",& mu_2dIP);
  oldtree->SetBranchAddress("mu_2dIPSig",& mu_2dIPSig);

  oldtree->SetBranchAddress("mu_FirstGenMatch", &mu_FirstGenMatch);
  oldtree->SetBranchAddress("mu_SecondGenMatch", &mu_SecondGenMatch);

   oldtree->SetBranchAddress("sv_munTracks", &sv_munTracks);
   oldtree->SetBranchAddress("sv_LxySig", &sv_LxySig);
   oldtree->SetBranchAddress("sv_LxyzSig", &sv_LxyzSig);
   oldtree->SetBranchAddress("sv_Lxy", &sv_Lxy);
   oldtree->SetBranchAddress("sv_Lxyz", &sv_Lxyz);
   oldtree->SetBranchAddress("sv_mass", &sv_mass);
   oldtree->SetBranchAddress("sv_eta", &sv_eta);
   oldtree->SetBranchAddress("sv_phi", &sv_phi);
   oldtree->SetBranchAddress("sv_pt", &sv_pt);
   oldtree->SetBranchAddress("sv_p", &sv_p);
   oldtree->SetBranchAddress("sv_px", &sv_px);
   oldtree->SetBranchAddress("sv_py", &sv_py);
   oldtree->SetBranchAddress("sv_pz", &sv_pz);
   oldtree->SetBranchAddress("sv_energy", &sv_energy);
   oldtree->SetBranchAddress("sv_Beta", &sv_Beta);
   oldtree->SetBranchAddress("sv_Gamma", &sv_Gamma);
   oldtree->SetBranchAddress("sv_CTau0", &sv_CTau0);
   oldtree->SetBranchAddress("sv_NDof", &sv_NDof);
   oldtree->SetBranchAddress("sv_Chi2", &sv_Chi2);
   oldtree->SetBranchAddress("sv_Angle3D", &sv_Angle3D);
   oldtree->SetBranchAddress("sv_Angle2D", &sv_Angle2D);
   oldtree->SetBranchAddress("sv_tracks_Sumcharge", &sv_tracks_Sumcharge);
   oldtree->SetBranchAddress("sv_tracks_Sumpt", &sv_tracks_Sumpt);
   oldtree->SetBranchAddress("sv_match_dxyz", &sv_match_dxyz);
   oldtree->SetBranchAddress("sv_X" , &sv_X);
   oldtree->SetBranchAddress("sv_Y" , &sv_Y);
   oldtree->SetBranchAddress("sv_Z" , &sv_Z);
   oldtree->SetBranchAddress("sv_xErr" , &sv_xErr);
   oldtree->SetBranchAddress("sv_yErr" , &sv_yErr);
   oldtree->SetBranchAddress("sv_zErr" , &sv_zErr);
   oldtree->SetBranchAddress("sv_lx" , &sv_lx);
   oldtree->SetBranchAddress("sv_ly" , &sv_ly);
   oldtree->SetBranchAddress("sv_lz" , &sv_lz);
   oldtree->SetBranchAddress("sv_tracks_charge" , &sv_tracks_charge);
   oldtree->SetBranchAddress("sv_tracks_eta" , &sv_tracks_eta);
   oldtree->SetBranchAddress("sv_tracks_phi" , &sv_tracks_phi);
   oldtree->SetBranchAddress("sv_tracks_pt" , &sv_tracks_pt);
   oldtree->SetBranchAddress("sv_tracks_dxySig" , &sv_tracks_dxySig);
   oldtree->SetBranchAddress("sv_tracks_dxy" , &sv_tracks_dxy);
   oldtree->SetBranchAddress("sv_tracks_dxyz" , &sv_tracks_dxyz);
   oldtree->SetBranchAddress("sv_tracks_p" , &sv_tracks_p);
   oldtree->SetBranchAddress("sv_hasMuon", & sv_hasMuon);


  oldtree->SetBranchAddress("pfMet_et", &pfMet_et);
  oldtree->SetBranchAddress("pfMet_pt", &pfMet_pt);
  oldtree->SetBranchAddress("pfMet_phi", &pfMet_phi);
  oldtree->SetBranchAddress("pfMet_en", &pfMet_en);
  oldtree->SetBranchAddress("pfMet_sumEt", &pfMet_sumEt);
  oldtree->SetBranchAddress("caloMet_pt", &caloMet_pt);
  oldtree->SetBranchAddress("caloMet_phi", &caloMet_phi);
  oldtree->SetBranchAddress("pfMet_px" , &pfMet_px);
  oldtree->SetBranchAddress("pfMet_py" , &pfMet_py);
  oldtree->SetBranchAddress("pfMet_pz" , &pfMet_pz);

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

  oldtree->SetBranchAddress("jetSmearedPt", &jetSmearedPt);

  oldtree->SetBranchAddress("jet_ptuncorrected", &jet_ptuncorrected);

  oldtree->SetBranchAddress("jet_CsvV2" ,&jet_CsvV2);
  oldtree->SetBranchAddress("jet_DeepCsv_udsg" ,&jet_DeepCsv_udsg);
  oldtree->SetBranchAddress("jet_DeepCsv_b" ,&jet_DeepCsv_b);
  oldtree->SetBranchAddress("jet_DeepCsv_c" ,&jet_DeepCsv_c);
  oldtree->SetBranchAddress("jet_DeepCsv_bb" ,&jet_DeepCsv_bb);
  oldtree->SetBranchAddress("jet_HadronFlavor" ,&jet_HadronFlavor);

  //systematic
  vector<Float_t>   *jetPt_JECUp= 0;
  vector<Float_t>   *jetPt_JECDown = 0;
  vector<Float_t>   *jetSmearedPt_unUp =0;
  vector<Float_t>   *jetSmearedPt_unDown = 0 ;
  vector<Float_t>   *prefire_weight=0;
  vector<Float_t>   *prefire_weightup=0;
  vector<Float_t>   *prefire_weightdown=0;

  oldtree->SetBranchAddress("jetPt_JECUp",&jetPt_JECUp);
  oldtree->SetBranchAddress("jetPt_JECDown",&jetPt_JECDown);
  oldtree->SetBranchAddress("jetSmearedPt_unUp",&jetSmearedPt_unUp);
  oldtree->SetBranchAddress("jetSmearedPt_unDown",&jetSmearedPt_unDown);
  oldtree->SetBranchAddress("prefire_weight",&prefire_weight);
  oldtree->SetBranchAddress("prefire_weightup",&prefire_weightup);
  oldtree->SetBranchAddress("prefire_weightdown",&prefire_weightdown);

  //Throw away some branches -- thos are all the  one I don't want to keep in the ntuples
  //IMPORTANT: we cannot throw away the branches with the variable we are using in the loop!
  oldtree->SetBranchStatus("ele_*",    0);
  //oldtree->SetBranchStatus("mu*", 0);

  //================================ pv info ============================================//                                                                                                                   
  Float_t pvX_ ,pvY_ , pvZ_ , pvXErr_ , pvYErr_, pvZErr_ , pvLxy_  , pvLxyz_ , pvLxySig_  , pvLxyzSig_ , pvChi2_,  pvSumPtSq_ ;
  //Int_t   pvNTrack_, numberPV_ ;

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
  //TBranch* branch_pvNTrack_   = newtree->Branch("pvNTrack_", &pvNTrack_, "pvNTrack/i");
  TBranch* branch_pvSumPtSq_  = newtree->Branch("pvSumPtSq_", &pvSumPtSq_, "pvSumPtSq_/F");
  //TBranch* branch_numberPV_   = newtree->Branch("numberPV_", &numberPV_, "numberPV/i");

  //================================ PU Weight ============================================//                                                                                                                 
  Float_t TrueNumInteractions;
  Float_t mc_gen_weight ;
  Float_t mc_lhe_weight ;
  Float_t mc_lhe_ctau ;

  TBranch* branch_TrueNumInteractions     = newtree->Branch("TrueNumInteractions",&TrueNumInteractions,"TrueNumInteractions/F");
  TBranch* branch_mc_gen_weight = newtree->Branch("mc_gen_weight",&mc_gen_weight,"mc_gen_weight/F");
  TBranch* branch_mc_lhe_weight = newtree->Branch("mc_lhe_weight",&mc_lhe_weight,"mc_lhe_weight/F");
  TBranch* branch_mc_lhe_ctau   = newtree->Branch("mc_lhe_ctau",&mc_lhe_ctau,"mc_lhe_ctau/F");


//======================= First Muon Variables ==========================================//   
Float_t mu_promptPt,mu_promptEta, mu_promptPhi, mu_promptCharge, mu_promptEt, mu_promptE;
 Float_t mu_prompt3dIP, mu_prompt3dIPSig, mu_prompt2dIP, mu_prompt2dIPSig ;

Float_t mu_promptRhoIso,        mu_promptTrackiso,       mu_promptPfSumChHadPt,
	mu_promptPFSumPhotonEt, mu_promptPfSumPUPt,      mu_promptEmIso,
	mu_promptHadIso,        mu_promptNormalizedChi2, mu_promptDPToverPT,
	mu_promptAbsdxy,        mu_promptAbsdxyError,    mu_promptAbsdxySig,
	mu_promptAbsdz,         mu_promptAbsdzError,     mu_promptAbsdzSig,
	mu_promptRecoDeltaBeta, mu_promptRecoiso,        mu_promptDirection,
	mu_promptNDof,          mu_promptTimeAtIpInOut,  mu_promptTimeAtIpInOutErr,
        mu_promptTimeAtIpOutIn, mu_promptTimeAtIpOutInErr,  mu_promptPfSumNHadEt, mu_promptGlobalMuon;

 Float_t mu_promptInnerTrackFraction, mu_promptSegmentCompatibility,mu_promptMotherID,
        mu_promptTrkKink,    mu_promptChi2LocalPosition ,RelIso_prompt;

 Int_t   mu_promptValidMuonHits, mu_promptMatchedStations, mu_promptValidPixelHits, 
         mu_promptTrackQuality,  mu_promptInrTrackQuality, mu_promptGenMatch,
         mu_promptPixelLayers ,  mu_promptTrackerLayers;

Float_t mu_promptRPCTofDirection,        mu_promptRPCTofNDof,          mu_promptRPCTofTimeAtIpInOut ,
        mu_promptRPCTofTimeAtIpInOutErr, mu_promptRPCTofTimeAtIpOutIn, mu_promptRPCTofTimeAtIpOutInErr;
 Int_t mu_promptsize;

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

 TBranch* branch_mu_promptsize = newtree->Branch("mu_promptsize", &mu_promptsize ,"mu_promptsize/I"); 

 TBranch* branch_RelIso_prompt  = newtree->Branch("RelIso_prompt", &RelIso_prompt ,"RelIso_prompt/F");

 TBranch* branch_mu_promptMotherID = newtree->Branch("mu_promptMotherID" ,&mu_promptMotherID ,"mu_promptMotherID/F");

 TBranch* branch_mu_prompt3dIP = newtree->Branch("mu_prompt3dIP" ,& mu_prompt3dIP ,"mu_prompt3dIP/F");
 TBranch* branch_mu_prompt3dIPSig = newtree->Branch("mu_prompt3dIPSig" ,& mu_prompt3dIPSig ,"mu_prompt3dIPSig/F");
 TBranch* branch_mu_prompt2dIP = newtree->Branch("mu_prompt2dIP" ,& mu_prompt2dIP ,"mu_prompt2dIP/F");
 TBranch* branch_mu_prompt2dIPSig  = newtree->Branch("mu_prompt2dIPSig" ,& mu_prompt2dIPSig ,"mu_prompt2dIPSig/F");

//======================= Second Muon Variables ==========================================//
Float_t mu_secondPt,mu_secondEta, mu_secondPhi, mu_secondCharge, mu_secondEt, mu_secondE;
 Float_t mu_second3dIP, mu_second3dIPSig, mu_second2dIP, mu_second2dIPSig ;

Float_t mu_secondRhoIso,        mu_secondTrackiso,  mu_secondPfSumChHadPt,
        mu_secondPFSumPhotonEt, mu_secondPfSumPUPt,      mu_secondEmIso,        
        mu_secondHadIso,        mu_secondNormalizedChi2, mu_secondDPToverPT,   
        mu_secondAbsdxy,        mu_secondAbsdxyError,    mu_secondAbsdxySig, 
        mu_secondAbsdz,         mu_secondAbsdzError,     mu_secondAbsdzSig,
        mu_secondRecoDeltaBeta, mu_secondRecoiso,        mu_secondDirection,
        mu_secondNDof,          mu_secondTimeAtIpInOut,  mu_secondTimeAtIpInOutErr, 
        mu_secondTimeAtIpOutIn, mu_secondTimeAtIpOutInErr,  mu_secondPfSumNHadEt, mu_secondGlobalMuon;

 Float_t mu_secondInnerTrackFraction, mu_secondSegmentCompatibility,mu_secondMotherID,
        mu_secondTrkKink,            mu_secondChi2LocalPosition, RelIso_second;


 Float_t mu_DeltaBetaR3, mu_DiMuMass, mu_DeltaR, mu_mT, vtxmu_mass,mu_DeltaPhi,mu_DeltaR_vec;

 Int_t   mu_secondValidMuonHits, mu_secondMatchedStations, mu_secondValidPixelHits, 
         mu_secondTrackQuality,  mu_secondInrTrackQuality, mu_secondGenMatch,
   mu_secondPixelLayers ,  mu_secondTrackerLayers, mu_secondIsTight , mu_nbLoose,mu_Size;

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
TBranch* branch_mu_DeltaR_vec  = newtree->Branch("mu_DeltaR_vec", &mu_DeltaR_vec, "mu_DeltaR_vec/F");
TBranch* branch_mu_DeltaPhi    = newtree->Branch("mu_DeltaPhi", &mu_DeltaPhi, "mu_DeltaPhi/F");
TBranch* branch_mu_Size        = newtree->Branch("mu_Size", &mu_Size, "mu_Size/I");
TBranch* branch_mu_nbLoose     = newtree->Branch("mu_nbLoose", &mu_nbLoose, "mu_nbLoose/I");
TBranch* branch_mu_mT          = newtree->Branch("mu_mT", &mu_mT, "mu_mT/F");
TBranch* branch_vtxmu_mass     = newtree->Branch("vtxmu_mass", &vtxmu_mass, "vtxmu_mass/F");   

 TBranch* branch_RelIso_second  = newtree->Branch("RelIso_second", &RelIso_second ,"RelIso_second/F");

 TBranch* branch_mu_secondMotherID = newtree->Branch("mu_secondMotherID" ,&mu_secondMotherID ,"mu_secondMotherID/F");

 TBranch* branch_mu_second3dIP     = newtree->Branch("mu_second3dIP" ,& mu_second3dIP ,"mu_second3dIP/F");
 TBranch* branch_mu_second3dIPSig  = newtree->Branch("mu_second3dIPSig" ,& mu_second3dIPSig ,"mu_second3dIPSig/F");
 TBranch* branch_mu_second2dIP     = newtree->Branch("mu_second2dIP" ,& mu_second2dIP ,"mu_second2dIP/F");
 TBranch* branch_mu_second2dIPSig  = newtree->Branch("mu_second2dIPSig" ,& mu_second2dIPSig ,"mu_second2dIPSig/F");

//======================= Tracks in Second Vertex Variables ==========================================//

 vector<float> sv_mass_n , sv_eta_n, sv_phi_n , sv_pt_n , sv_p_n , sv_px_n , sv_py_n , 
               sv_pz_n , sv_energy_n , sv_Beta_n , sv_Gamma_n , sv_CTau0_n , sv_NDof_n , 
               sv_Chi2_n , sv_Angle3D_n , sv_Angle2D_n;

 vector<float> sv_tracks_eta_n , sv_tracks_phi_n , sv_tracks_pt_n , sv_tracks_p_n , 
               sv_tracks_dxySig_n , sv_tracks_dxy_n , sv_tracks_dxyz_n;

 vector<int>   sv_tracks_charge_n;

TBranch* branch_sv_mass_n   =   newtree->Branch("sv_mass_n", &sv_mass_n);
TBranch* branch_sv_eta_n    =   newtree->Branch("sv_eta_n", &sv_eta_n);
TBranch* branch_sv_phi_n    =   newtree->Branch("sv_phi_n", &sv_phi_n);
TBranch* branch_sv_pt_n     =   newtree->Branch("sv_pt_n", &sv_pt_n);
TBranch* branch_sv_p_n      =   newtree->Branch("sv_p_n", &sv_p_n);
TBranch* branch_sv_px_n     =   newtree->Branch("sv_px_n", &sv_px_n);
TBranch* branch_sv_py_n     =   newtree->Branch("sv_py_n", &sv_py_n);
TBranch* branch_sv_pz_n     =   newtree->Branch("sv_pz_n", &sv_pz_n);
TBranch* branch_sv_energy_n =   newtree->Branch("sv_energy_n", &sv_energy_n);
TBranch* branch_sv_Beta_n   =   newtree->Branch("sv_Beta_n", &sv_Beta_n);
TBranch* branch_sv_Gamma_n  =   newtree->Branch("sv_Gamma_n", &sv_Gamma_n);
TBranch* branch_sv_CTau0_n  =   newtree->Branch("sv_CTau0_n", &sv_CTau0_n);
TBranch* branch_sv_NDof_n   =   newtree->Branch("sv_NDof_n", &sv_NDof_n);
TBranch* branch_sv_Chi2_n    =   newtree->Branch("sv_Chi2_n", &sv_Chi2_n);
TBranch* branch_sv_Angle3D_n =   newtree->Branch("sv_Angle3D_n", &sv_Angle3D_n);
TBranch* branch_sv_Angle2D_n =   newtree->Branch("sv_Angle2D_n", &sv_Angle2D_n);
TBranch* branch_sv_tracks_charge_n =   newtree->Branch("sv_tracks_charge_n", &sv_tracks_charge_n);
TBranch* branch_sv_tracks_eta_n    =   newtree->Branch("sv_tracks_eta_n", &sv_tracks_eta_n);
TBranch* branch_sv_tracks_phi_n    =   newtree->Branch("sv_tracks_phi_n", &sv_tracks_phi_n);
TBranch* branch_sv_tracks_pt_n     =   newtree->Branch("sv_tracks_pt_n", &sv_tracks_pt_n);
TBranch* branch_sv_tracks_p_n      =   newtree->Branch("sv_tracks_p_n", &sv_tracks_p_n);
TBranch* branch_sv_tracks_dxySig_n =   newtree->Branch("sv_tracks_dxySig_n", &sv_tracks_dxySig_n);
TBranch* branch_sv_tracks_dxy_n    =   newtree->Branch("sv_tracks_dxy_n", &sv_tracks_dxy_n);
TBranch* branch_sv_tracks_dxyz_n   =   newtree->Branch("sv_tracks_dxyz_n", &sv_tracks_dxyz_n);

 //======================= Second Vertex Variables ==========================================//  
 Float_t  sv_LXYSig_,  sv_LXYZSig_, sv_LXY_,  sv_LXYZ_,  sv_mass_,
   sv_eta_,     sv_phi_,     sv_pt_,   sv_p_,     sv_Beta_,
   sv_Gamma_,   sv_CTau0_,   sv_NDof_, sv_Chi2_,  sv_Angle3D_,
   sv_Angle2D_, sv_tracks_Sumpt_,     sv_match_, sv_energy_,
   sv_px_,      sv_py_,      sv_pz_,   sv_dir_x_, sv_dir_y_, sv_dir_z_;

 Float_t sv_Xpos_,  sv_Ypos_,  sv_Zpos_,  sv_xError_,   sv_yError_,   sv_zError_, sv_lx_, sv_ly_, sv_lz_;
 Float_t sv_track_sumdxySig_;
 Int_t sv_hasMuon_;
 Int_t sv_TrackSize_ , sv_tracks_Sumcharge_ ;
 unsigned sv_inside_;

 TBranch* branch_sv_inside     = newtree->Branch("sv_inside", &sv_inside_, "sv_inside/I");
 TBranch* branch_sv_TrackSize  = newtree->Branch("sv_TrackSize", &sv_TrackSize_, "sv_TrackSize/I");
 TBranch* branch_sv_LXYSig     = newtree->Branch("sv_LXYSig", &sv_LXYSig_, "sv_LXYSig/F");
 TBranch* branch_sv_LXYZSig    = newtree->Branch("sv_LXYZSig", &sv_LXYZSig_, "sv_LXYZSig/F");
 TBranch* branch_sv_LXY        = newtree->Branch("sv_LXY", &sv_LXY_, "sv_LXY/F");
 TBranch* branch_sv_LXYZ       = newtree->Branch("sv_LXYZ", &sv_LXYZ_, "sv_LXYZ/F");
 TBranch* branch_sv_dir_x      = newtree->Branch("sv_dir_x", &sv_dir_x_, "sv_dir_x/F");
 TBranch* branch_sv_dir_y      = newtree->Branch("sv_dir_y", &sv_dir_y_, "sv_dir_y/F");
 TBranch* branch_sv_dir_z      = newtree->Branch("sv_dir_z", &sv_dir_z_, "sv_dir_z/F");
 TBranch* branch_sv_mass       = newtree->Branch("sv_mass", &sv_mass_, "sv_mass/F");
 TBranch* branch_sv_eta        = newtree->Branch("sv_eta", &sv_eta_, "sv_eta/F");
 TBranch* branch_sv_phi        = newtree->Branch("sv_phi", &sv_phi_, "sv_phi/F");
 TBranch* branch_sv_pt         = newtree->Branch("sv_pt", &sv_pt_, "sv_pt/F");
 TBranch* branch_sv_p          = newtree->Branch("sv_p", &sv_p_, "sv_p/F");
 TBranch* branch_sv_px         = newtree->Branch("sv_px", &sv_px_, "sv_px/F");
 TBranch* branch_sv_py         = newtree->Branch("sv_py", &sv_py_, "sv_py/F");
 TBranch* branch_sv_pz         = newtree->Branch("sv_pz", &sv_pz_, "sv_pz/F");
 TBranch* branch_sv_energy     = newtree->Branch("sv_energy", &sv_energy_, "sv_energy/F");
 TBranch* branch_sv_Beta       = newtree->Branch("sv_Beta", &sv_Beta_, "sv_Beta/F");
 TBranch* branch_sv_Gamma      = newtree->Branch("sv_Gamma", &sv_Gamma_, "sv_Gamma/F");
 TBranch* branch_sv_CTau0      = newtree->Branch("sv_CTau0", &sv_CTau0_, "sv_CTau0/F");
 TBranch* branch_sv_NDof       = newtree->Branch("sv_NDof", &sv_NDof_, "sv_NDof/F");
 TBranch* branch_sv_Chi2       = newtree->Branch("sv_Chi2", &sv_Chi2_, "sv_Chi2/F");
 TBranch* branch_sv_Angle3D    = newtree->Branch("sv_Angle3D", &sv_Angle3D_, "sv_Angle3D/F");
 TBranch* branch_sv_Angle2D    = newtree->Branch("sv_Angle2D", &sv_Angle2D_, "sv_Angle2D/F");
 TBranch* branch_sv_tracks_Sumcharge = newtree->Branch("sv_tracks_Sumcharge", &sv_tracks_Sumcharge_, "sv_tracks_Sumcharge/I");
 TBranch* branch_sv_tracks_Sumpt     = newtree->Branch("sv_tracks_Sumpt", &sv_tracks_Sumpt_, "sv_tracks_Sumpt/F");
 TBranch* branch_sv_match            = newtree->Branch("sv_match", &sv_match_, "sv_match/F");

 TBranch* branch_sv_Xpos     = newtree->Branch("sv_Xpos" , &sv_Xpos_,"sv_Xpos/F");
 TBranch* branch_sv_Ypos     = newtree->Branch("sv_Ypos" , &sv_Ypos_,"sv_Ypos/F");
 TBranch* branch_sv_Zpos     = newtree->Branch("sv_Zpos" , &sv_Zpos_,"sv_Zpos/F");
 TBranch* branch_sv_xError   = newtree->Branch("sv_xError" , &sv_xError_,"sv_xError/F");
 TBranch* branch_sv_yError   = newtree->Branch("sv_yError" , &sv_yError_,"sv_yError/F");
 TBranch* branch_sv_zError   = newtree->Branch("sv_zError" , &sv_zError_,"sv_zError/F");
 TBranch* branch_sv_lx       = newtree->Branch("sv_lx" , &sv_lx_,"sv_lx/F");
 TBranch* branch_sv_ly       = newtree->Branch("sv_ly" , &sv_ly_,"sv_ly/F");
 TBranch* branch_sv_lz       = newtree->Branch("sv_lz" , &sv_lz_,"sv_lz/F");
 TBranch* branch_sv_hasMuon  = newtree->Branch("sv_hasMuon", &sv_hasMuon_, "sv_hasMuon/I");

 TBranch* branch_sv_track_sumdxySig = newtree->Branch("sv_track_sumdxySig" , &sv_track_sumdxySig_, "sv_track_sumdxySig/F");

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

 //================================= Jets Variables ==============================================//
 Float_t alljet_charge_,         alljet_et_,    alljet_pt_, alljet_eta_, alljet_NEmEnFraction_,
         alljet_phi_,            alljet_theta_, alljet_en_, alljet_chEmEn_, alljet_NHadEnFraction_,
         alljet_chMuEn_,         alljet_chMuEnFraction_,    alljet_numberOfDaughters_,
         alljet_muonEnergyFraction_, alljet_muonMultiplicity_,
         alljet_neutralEmEnergy_,    alljet_neutralHadronEnergy_, alljet_NHadMultiplicity_,
         alljet_NMultiplicity_,      alljet_chHadEn_,       alljet_muonEnergy_, alljet_chargedMultiplicity_ ,
   alljet_chargedEmEnergyFraction, alljet_chargedHadronEnergyFraction,alljet_notSmeard_pt,alljet_raw_pt,HT_alljet;

 Float_t alljet_CsvV2, alljet_DeepCsv_udsg, alljet_DeepCsv_b, alljet_DeepCsv_c, alljet_DeepCsv_bb, alljet_HadronFlavor;
 Int_t   alljets_size, alljet_bjet_L,alljet_bjet_M,alljet_bjet_T;

 TBranch* branch_alljet_charge_ = newtree->Branch("alljet_charge", &alljet_charge_, "alljet_charge/F");
 TBranch* branch_alljet_et_     = newtree->Branch("alljet_et", &alljet_et_, "alljet_et/F");
 TBranch* branch_alljet_pt_     = newtree->Branch("alljet_pt", &alljet_pt_, "alljet_pt/F");
 TBranch* branch_alljet_eta_    = newtree->Branch("alljet_eta", &alljet_eta_, "alljet_eta/F");
 TBranch* branch_alljet_phi_    = newtree->Branch("alljet_phi", &alljet_phi_, "alljet_phi/F");
 TBranch* branch_alljet_theta_  = newtree->Branch("alljet_theta", &alljet_theta_, "alljet_theta/F");
 TBranch* branch_alljet_en_     = newtree->Branch("alljet_en", &alljet_en_, "alljet_en/F");
 TBranch* branch_alljet_chEmEn_                     = newtree->Branch("alljet_chEmEn", &alljet_chEmEn_, "alljet_chEmEn/F");
 TBranch* branch_alljet_NEmEnFraction_              = newtree->Branch("alljet_NEmEnFraction", &alljet_NEmEnFraction_, "alljet_NEmEnFraction/F");
 TBranch* branch_alljet_chHadEn_                    = newtree->Branch("alljet_chHadEn", &alljet_chHadEn_, "alljet_chHadEn/F");
 TBranch* branch_alljet_NHadEnFraction_             = newtree->Branch("alljet_NHadEnFraction", &alljet_NHadEnFraction_, "alljet_NHadEnFraction/F");
 TBranch* branch_alljet_chMuEn_                     = newtree->Branch("alljet_chMuEn", &alljet_chMuEn_, "alljet_chMuEn/F");
 TBranch* branch_alljet_chMuEnFraction_             = newtree->Branch("alljet_chMuEnFraction", &alljet_chMuEnFraction_, "alljet_chMuEnFraction/F");
 TBranch* branch_alljet_numberOfDaughters_          = newtree->Branch("alljet_numberOfDaughters", &alljet_numberOfDaughters_, "alljet_numberOfDaughters/F");
 TBranch* branch_alljet_muonEnergy_                 = newtree->Branch("alljet_muonEnergy", &alljet_muonEnergy_, "alljet_muonEnergy/F");
 TBranch* branch_alljet_muonEnergyFraction_         = newtree->Branch("alljet_muonEnergyFraction", &alljet_muonEnergyFraction_, "alljet_muonEnergyFraction/F");
 TBranch* branch_alljet_muonMultiplicity_           = newtree->Branch("alljet_muonMultiplicity", &alljet_muonMultiplicity_, "alljet_muonMultiplicity/F");
 TBranch* branch_alljet_neutralEmEnergy_            = newtree->Branch("alljet_neutralEmEnergy", &alljet_neutralEmEnergy_, "alljet_neutralEmEnergy/F");
 TBranch* branch_alljet_neutralHadronEnergy_        = newtree->Branch("alljet_neutralHadronEnergy", &alljet_neutralHadronEnergy_, "alljet_neutralHadronEnergy/F");
 TBranch* branch_alljet_NHadMultiplicity_           = newtree->Branch("alljet_NHadMultiplicity", &alljet_NHadMultiplicity_, "alljet_NHadMultiplicity/F");
 TBranch* branch_alljet_NMultiplicity_              = newtree->Branch("alljet_NMultiplicity", &alljet_NMultiplicity_, "alljet_NMultiplicity/F");
 TBranch* branch_alljet_chargedMultiplicity_        = newtree->Branch("alljet_chargedMultiplicity", &alljet_chargedMultiplicity_, "alljet_chargedMultiplicity/F");
 TBranch* branch_alljet_chargedEmEnergyFraction     = newtree->Branch("alljet_chargedEmEnergyFraction", &alljet_chargedEmEnergyFraction, "alljet_chargedEmEnergyFraction/F");
 TBranch* branch_alljet_chargedHadronEnergyFraction = newtree->Branch("alljet_chargedHadronEnergyFraction", &alljet_chargedHadronEnergyFraction, "alljet_chargedHadronEnergyFraction/F");

 TBranch* branch_alljets_size =  newtree->Branch("alljets_size", &alljets_size, "alljets_size/I");

 TBranch* branch_alljet_bjet_L =  newtree->Branch("alljet_bjet_L", &alljet_bjet_L, "alljet_bjet_L/I");
 TBranch* branch_alljet_bjet_M =  newtree->Branch("alljet_bjet_M", &alljet_bjet_M, "alljet_bjet_M/I");
 TBranch* branch_alljet_bjet_T =  newtree->Branch("alljet_bjet_T", &alljet_bjet_T, "alljet_bjet_T/I");

 TBranch* branchalltjet_notSmeard_pt                = newtree->Branch("alljet_notSmeard_pt" ,&alljet_notSmeard_pt,"alljet_notSmeard_pt/F");
 TBranch* branch_alljet_raw_pt                      = newtree->Branch("alljet_raw_pt",&alljet_notSmeard_pt,"alljet_notSmeard_pt/F");

 TBranch* branch_alljet_CsvV2        =  newtree->Branch("alljet_CsvV2", &alljet_CsvV2, "alljet_CsvV2/F");
 TBranch* branch_alljet_DeepCsv_udsg =  newtree->Branch("alljet_DeepCsv_udsg", &alljet_DeepCsv_udsg, "alljet_DeepCsv_udsg/F");
 TBranch* branch_alljet_DeepCsv_b    =  newtree->Branch("alljet_DeepCsv_b", &alljet_DeepCsv_b, "alljet_DeepCsv_b/F");
 TBranch* branch_alljet_DeepCsv_c    =  newtree->Branch("alljet_DeepCsv_c", &alljet_DeepCsv_c,"alljet_DeepCsv_c/F");
 TBranch* branch_alljet_DeepCsv_bb   =  newtree->Branch("alljet_DeepCsv_bb", &alljet_DeepCsv_bb,"alljet_DeepCsv_bb/F");
 TBranch* branch_alljet_HadronFlavor =  newtree->Branch("alljet_HadronFlavor", &alljet_HadronFlavor,"alljet_HadronFlavor/F");
 TBranch* branch_HT_alljet           =  newtree->Branch("HT_alljet" ,&HT_alljet,"HT_alljet/F");

 //================================= Second Jet Variables ==============================================//                                                                                                 
 Float_t firstjet_charge_,          firstjet_et_,    firstjet_pt_, firstjet_eta_, firstjet_NEmEnFraction_,
         firstjet_phi_,             firstjet_theta_, firstjet_en_, firstjet_chEmEn_, firstjet_NHadEnFraction_,
         firstjet_chMuEn_,          firstjet_chMuEnFraction_,      firstjet_numberOfDaughters_,
         firstjet_muonEnergyFraction_, firstjet_muonMultiplicity_,
         firstjet_neutralEmEnergy_,    firstjet_neutralHadronEnergy_,           firstjet_NHadMultiplicity_,
         firstjet_NMultiplicity_,      firstjet_chHadEn_,          firstjet_muonEnergy_, firstjet_chargedMultiplicity_ ,
   firstjet_chargedEmEnergyFraction, firstjet_chargedHadronEnergyFraction, firstjet_notSmeard_pt, firstjet_raw_pt;

 Float_t firstjet_CsvV2, firstjet_DeepCsv_udsg, firstjet_DeepCsv_b, firstjet_DeepCsv_c, firstjet_DeepCsv_bb, firstjet_HadronFlavor;
 Int_t   firstjet_bjet_L,firstjet_bjet_M,firstjet_bjet_T,firstjets_size;

 TBranch* branch_firstjet_charge_   = newtree->Branch("firstjet_charge", &firstjet_charge_, "firstjet_charge/F");
 TBranch* branch_firstjet_et_       = newtree->Branch("firstjet_et", &firstjet_et_, "firstjet_et/F");
 TBranch* branch_firstjet_pt_       = newtree->Branch("firstjet_pt", &firstjet_pt_, "firstjet_pt/F");
 TBranch* branch_firstjet_eta_      = newtree->Branch("firstjet_eta", &firstjet_eta_, "firstjet_eta/F");
 TBranch* branch_firstjet_phi_      = newtree->Branch("firstjet_phi", &firstjet_phi_, "firstjet_phi/F");
 TBranch* branch_firstjet_theta_    = newtree->Branch("firstjet_theta", &firstjet_theta_, "firstjet_theta/F");
 TBranch* branch_firstjet_en_       = newtree->Branch("firstjet_en", &firstjet_en_, "firstjet_en/F");
 TBranch* branch_firstjet_chEmEn_                     = newtree->Branch("firstjet_chEmEn", &firstjet_chEmEn_, "firstjet_chEmEn/F");
 TBranch* branch_firstjet_NEmEnFraction_              = newtree->Branch("firstjet_NEmEnFraction", &firstjet_NEmEnFraction_, "firstjet_NEmEnFraction/F");
 TBranch* branch_firstjet_chHadEn_                    = newtree->Branch("firstjet_chHadEn", &firstjet_chHadEn_, "firstjet_chHadEn/F");
 TBranch* branch_firstjet_NHadEnFraction_             = newtree->Branch("firstjet_NHadEnFraction", &firstjet_NHadEnFraction_, "firstjet_NHadEnFraction/F");
 TBranch* branch_firstjet_chMuEn_                     = newtree->Branch("firstjet_chMuEn", &firstjet_chMuEn_, "firstjet_chMuEn/F");
 TBranch* branch_firstjet_chMuEnFraction_             = newtree->Branch("firstjet_chMuEnFraction", &firstjet_chMuEnFraction_, "firstjet_chMuEnFraction/F");
 TBranch* branch_firstjet_numberOfDaughters_          = newtree->Branch("firstjet_numberOfDaughters", &firstjet_numberOfDaughters_, "firstjet_numberOfDaughters/F");
 TBranch* branch_firstjet_muonEnergy_                 = newtree->Branch("firstjet_muonEnergy", &firstjet_muonEnergy_, "firstjet_muonEnergy/F");
 TBranch* branch_firstjet_muonEnergyFraction_         = newtree->Branch("firstjet_muonEnergyFraction", &firstjet_muonEnergyFraction_, "firstjet_muonEnergyFraction/F");
 TBranch* branch_firstjet_muonMultiplicity_           = newtree->Branch("firstjet_muonMultiplicity", &firstjet_muonMultiplicity_, "firstjet_muonMultiplicity/F");
 TBranch* branch_firstjet_neutralEmEnergy_            = newtree->Branch("firstjet_neutralEmEnergy", &firstjet_neutralEmEnergy_, "firstjet_neutralEmEnergy/F");
 TBranch* branch_firstjet_neutralHadronEnergy_        = newtree->Branch("firstjet_neutralHadronEnergy", &firstjet_neutralHadronEnergy_, "firstjet_neutralHadronEnergy/F");
 TBranch* branch_firstjet_NHadMultiplicity_           = newtree->Branch("firstjet_NHadMultiplicity", &firstjet_NHadMultiplicity_, "firstjet_NHadMultiplicity/F");
 TBranch* branch_firstjet_NMultiplicity_              = newtree->Branch("firstjet_NMultiplicity", &firstjet_NMultiplicity_, "firstjet_NMultiplicity/F");
 TBranch* branch_firstjet_chargedMultiplicity_        = newtree->Branch("firstjet_chargedMultiplicity", &firstjet_chargedMultiplicity_, "firstjet_chargedMultiplicity/F");
 TBranch* branch_firstjet_chargedEmEnergyFraction     = newtree->Branch("firstjet_chargedEmEnergyFraction", &firstjet_chargedEmEnergyFraction, "firstjet_chargedEmEnergyFraction/F");
 TBranch* branch_firstjet_chargedHadronEnergyFraction = newtree->Branch("firstjet_chargedHadronEnergyFraction", &firstjet_chargedHadronEnergyFraction,"firstjet_chargedHadronEnergyFraction/F");

 TBranch* branch_firstjets_size  =  newtree->Branch("firstjets_size",  &firstjets_size,  "firstjets_size/I");
 TBranch* branch_firstjet_bjet_L =  newtree->Branch("firstjet_bjet_L", &firstjet_bjet_L, "firstjet_bjet_L/I");
 TBranch* branch_firstjet_bjet_M =  newtree->Branch("firstjet_bjet_M", &firstjet_bjet_M, "firstjet_bjet_M/I");
 TBranch* branch_firstjet_bjet_T =  newtree->Branch("firstjet_bjet_T", &firstjet_bjet_T, "firstjet_bjet_T/I");

 TBranch* branch_firstjet_notSmeard_pt                = newtree->Branch("firstjet_notSmeard_pt" ,&firstjet_notSmeard_pt,"firstjet_notSmeard_pt/F");
 TBranch* branch_firstjet_raw_pt                      = newtree->Branch("firstjet_raw_pt",&firstjet_notSmeard_pt,"firstjet_notSmeard_pt/F");


 TBranch* branch_firstjet_CsvV2        =  newtree->Branch("firstjet_CsvV2", &firstjet_CsvV2, "firstjet_CsvV2/F");
 TBranch* branch_firstjet_DeepCsv_udsg =  newtree->Branch("firstjet_DeepCsv_udsg", &firstjet_DeepCsv_udsg, "firstjet_DeepCsv_udsg/F");
 TBranch* branch_firstjet_DeepCsv_b    =  newtree->Branch("firstjet_DeepCsv_b", &firstjet_DeepCsv_b, "firstjet_DeepCsv_b/F");
 TBranch* branch_firstjet_DeepCsv_c    =  newtree->Branch("firstjet_DeepCsv_c", &firstjet_DeepCsv_c,"firstjet_DeepCsv_c/F");
 TBranch* branch_firstjet_DeepCsv_bb   =  newtree->Branch("firstjet_DeepCsv_bb", &firstjet_DeepCsv_bb,"firstjet_DeepCsv_bb/F");
 TBranch* branch_firstjet_HadronFlavor =  newtree->Branch("firstjet_HadronFlavor", &firstjet_HadronFlavor,"firstjet_HadronFlavor/F");


 //================================= Thrid Jet Variables ==============================================//                                                                                                   
 Float_t mujet_charge_,             mujet_et_,             mujet_pt_, mujet_eta_,    mujet_NEmEnFraction_,
         mujet_phi_,                mujet_theta_,          mujet_en_, mujet_chEmEn_, mujet_NHadEnFraction_,
         mujet_chMuEn_,             mujet_chMuEnFraction_, mujet_numberOfDaughters_,
         mujet_muonEnergyFraction_, mujet_muonMultiplicity_,
         mujet_neutralEmEnergy_,    mujet_neutralHadronEnergy_, mujet_NHadMultiplicity_,
         mujet_NMultiplicity_,      mujet_chHadEn_,         mujet_muonEnergy_, mujet_chargedMultiplicity_ ,
   mujet_chargedEmEnergyFraction, mujet_chargedHadronEnergyFraction, mujet_notSmeard_pt, mujet_raw_pt;

 Float_t mujet_CsvV2, mujet_DeepCsv_udsg, mujet_DeepCsv_b, mujet_DeepCsv_c, mujet_DeepCsv_bb, mujet_HadronFlavor, mujet_RSecMu, mujet_M,hnl_mass;

 TBranch* branch_mujet_charge_                     = newtree->Branch("mujet_charge", &mujet_charge_, "mujet_charge/F");
 TBranch* branch_mujet_et_                         = newtree->Branch("mujet_et",     &mujet_et_, "mujet_et/F");
 TBranch* branch_mujet_pt_                         = newtree->Branch("mujet_pt",     &mujet_pt_, "mujet_pt/F");
 TBranch* branch_mujet_eta_                        = newtree->Branch("mujet_eta",    &mujet_eta_, "mujet_eta/F");
 TBranch* branch_mujet_phi_                        = newtree->Branch("mujet_phi",    &mujet_phi_, "mujet_phi/F");
 TBranch* branch_mujet_theta_                      = newtree->Branch("mujet_theta",  &mujet_theta_, "mujet_theta/F");
 TBranch* branch_mujet_en_                         = newtree->Branch("mujet_en",     &mujet_en_, "mujet_en/F");
 TBranch* branch_mujet_chEmEn_                     = newtree->Branch("mujet_chEmEn", &mujet_chEmEn_, "mujet_chEmEn/F");
 TBranch* branch_mujet_NEmEnFraction_              = newtree->Branch("mujet_NEmEnFraction", &mujet_NEmEnFraction_, "mujet_NEmEnFraction/F");
 TBranch* branch_mujet_chHadEn_                    = newtree->Branch("mujet_chHadEn",       &mujet_chHadEn_, "mujet_chHadEn/F");
 TBranch* branch_mujet_NHadEnFraction_             = newtree->Branch("mujet_NHadEnFraction",&mujet_NHadEnFraction_, "mujet_NHadEnFraction/F");
 TBranch* branch_mujet_chMuEn_                     = newtree->Branch("mujet_chMuEn",        &mujet_chMuEn_, "mujet_chMuEn/F");
 TBranch* branch_mujet_chMuEnFraction_             = newtree->Branch("mujet_chMuEnFraction",&mujet_chMuEnFraction_, "mujet_chMuEnFraction/F");
 TBranch* branch_mujet_numberOfDaughters_          = newtree->Branch("mujet_numberOfDaughters", &mujet_numberOfDaughters_, "mujet_numberOfDaughters/F");
 TBranch* branch_mujet_muonEnergy_                 = newtree->Branch("mujet_muonEnergy",        &mujet_muonEnergy_, "mujet_muonEnergy/F");
 TBranch* branch_mujet_muonEnergyFraction_         = newtree->Branch("mujet_muonEnergyFraction",&mujet_muonEnergyFraction_, "mujet_muonEnergyFraction/F");
 TBranch* branch_mujet_muonMultiplicity_           = newtree->Branch("mujet_muonMultiplicity",  &mujet_muonMultiplicity_, "mujet_muonMultiplicity/F");
 TBranch* branch_mujet_neutralEmEnergy_            = newtree->Branch("mujet_neutralEmEnergy",   &mujet_neutralEmEnergy_, "mujet_neutralEmEnergy/F");
 TBranch* branch_mujet_neutralHadronEnergy_        = newtree->Branch("mujet_neutralHadronEnergy",&mujet_neutralHadronEnergy_, "mujet_neutralHadronEnergy/F");
 TBranch* branch_mujet_NHadMultiplicity_           = newtree->Branch("mujet_NHadMultiplicity",   &mujet_NHadMultiplicity_, "mujet_NHadMultiplicity/F");
 TBranch* branch_mujet_NMultiplicity_              = newtree->Branch("mujet_NMultiplicity",      &mujet_NMultiplicity_, "mujet_NMultiplicity/F");
 TBranch* branch_mujet_chargedMultiplicity_        = newtree->Branch("mujet_chargedMultiplicity",         &mujet_chargedMultiplicity_, "mujet_chargedMultiplicity/F");
 TBranch* branch_mujet_chargedEmEnergyFraction     = newtree->Branch("mujet_chargedEmEnergyFraction",     &mujet_chargedEmEnergyFraction, "mujet_chargedEmEnergyFraction/F");
 TBranch* branch_mujet_chargedHadronEnergyFraction = newtree->Branch("mujet_chargedHadronEnergyFraction", &mujet_chargedHadronEnergyFraction,"mujet_chargedHadronEnergyFraction/F");

 TBranch* branch_mujet_notSmeard_pt                = newtree->Branch("mujet_notSmeard_pt" ,&mujet_notSmeard_pt,"mujet_notSmeard_pt/F");
 TBranch* branch_mujet_raw_pt                      = newtree->Branch("mujet_raw_pt",&mujet_notSmeard_pt,"mujet_notSmeard_pt/F");


 TBranch* branch_mujet_CsvV2        =  newtree->Branch("mujet_CsvV2", &mujet_CsvV2, "mujet_CsvV2/F");
 TBranch* branch_mujet_DeepCsv_udsg =  newtree->Branch("mujet_DeepCsv_udsg", &mujet_DeepCsv_udsg, "mujet_DeepCsv_udsg/F");
 TBranch* branch_mujet_DeepCsv_b    =  newtree->Branch("mujet_DeepCsv_b", &mujet_DeepCsv_b, "mujet_DeepCsv_b/F");
 TBranch* branch_mujet_DeepCsv_c    =  newtree->Branch("mujet_DeepCsv_c", &mujet_DeepCsv_c,"mujet_DeepCsv_c/F");
 TBranch* branch_mujet_DeepCsv_bb   =  newtree->Branch("mujet_DeepCsv_bb", &mujet_DeepCsv_bb,"mujet_DeepCsv_bb/F");
 TBranch* branch_mujet_HadronFlavor =  newtree->Branch("mujet_HadronFlavor", &mujet_HadronFlavor,"mujet_HadronFlavor/F");
 TBranch* branch_mujet_RSecMu       =  newtree->Branch("mujet_RSecMu", &mujet_RSecMu,"mujet_RSecMu/F");

 TBranch* branch_mujet_M  = newtree->Branch("mujet_M", &mujet_M,"mujet_M/F");
 TBranch* branch_hnl_mass = newtree->Branch("hnl_mass", &hnl_mass, "hnl_mass/F");

 //================================= Prompt Jet Variables ==============================================//  
 Float_t prompt_mujet_charge_,             prompt_mujet_et_,             prompt_mujet_pt_, prompt_mujet_eta_,    prompt_mujet_NEmEnFraction_,
   prompt_mujet_phi_,                prompt_mujet_theta_,          prompt_mujet_en_, prompt_mujet_chEmEn_, prompt_mujet_NHadEnFraction_,
   prompt_mujet_chMuEn_,             prompt_mujet_chMuEnFraction_, prompt_mujet_numberOfDaughters_,
   prompt_mujet_muonEnergyFraction_, prompt_mujet_muonMultiplicity_,
   prompt_mujet_neutralEmEnergy_,    prompt_mujet_neutralHadronEnergy_, prompt_mujet_NHadMultiplicity_,
   prompt_mujet_NMultiplicity_,      prompt_mujet_chHadEn_,         prompt_mujet_muonEnergy_, prompt_mujet_chargedMultiplicity_ ,
   prompt_mujet_chargedEmEnergyFraction, prompt_mujet_chargedHadronEnergyFraction, prompt_mujet_notSmeard_pt, prompt_mujet_raw_pt;

 Float_t prompt_mujet_CsvV2, prompt_mujet_DeepCsv_udsg, prompt_mujet_DeepCsv_b, prompt_mujet_DeepCsv_c, prompt_mujet_DeepCsv_bb, prompt_mujet_HadronFlavor, prompt_mujet_RSecMu, prompt_mujet_M;

 TBranch* branch_prompt_mujet_charge_                     = newtree->Branch("prompt_mujet_charge", &prompt_mujet_charge_, "prompt_mujet_charge/F");
 TBranch* branch_prompt_mujet_et_                         = newtree->Branch("prompt_mujet_et",     &prompt_mujet_et_, "prompt_mujet_et/F");
 TBranch* branch_prompt_mujet_pt_                         = newtree->Branch("prompt_mujet_pt",     &prompt_mujet_pt_, "prompt_mujet_pt/F");
 TBranch* branch_prompt_mujet_eta_                        = newtree->Branch("prompt_mujet_eta",    &prompt_mujet_eta_, "prompt_mujet_eta/F");
 TBranch* branch_prompt_mujet_phi_                        = newtree->Branch("prompt_mujet_phi",    &prompt_mujet_phi_, "prompt_mujet_phi/F");
 TBranch* branch_prompt_mujet_theta_                      = newtree->Branch("prompt_mujet_theta",  &prompt_mujet_theta_, "prompt_mujet_theta/F");
 TBranch* branch_prompt_mujet_en_                         = newtree->Branch("prompt_mujet_en",     &prompt_mujet_en_, "prompt_mujet_en/F");
 TBranch* branch_prompt_mujet_chEmEn_                     = newtree->Branch("prompt_mujet_chEmEn", &prompt_mujet_chEmEn_, "prompt_mujet_chEmEn/F");
 TBranch* branch_prompt_mujet_NEmEnFraction_              = newtree->Branch("prompt_mujet_NEmEnFraction", &prompt_mujet_NEmEnFraction_, "prompt_mujet_NEmEnFraction/F");
 TBranch* branch_prompt_mujet_chHadEn_                    = newtree->Branch("prompt_mujet_chHadEn",       &prompt_mujet_chHadEn_, "prompt_mujet_chHadEn/F");
 TBranch* branch_prompt_mujet_NHadEnFraction_             = newtree->Branch("prompt_mujet_NHadEnFraction",&prompt_mujet_NHadEnFraction_, "prompt_mujet_NHadEnFraction/F");
 TBranch* branch_prompt_mujet_chMuEn_                     = newtree->Branch("prompt_mujet_chMuEn",        &prompt_mujet_chMuEn_, "prompt_mujet_chMuEn/F");
 TBranch* branch_prompt_mujet_chMuEnFraction_             = newtree->Branch("prompt_mujet_chMuEnFraction",&prompt_mujet_chMuEnFraction_, "prompt_mujet_chMuEnFraction/F");
 TBranch* branch_prompt_mujet_numberOfDaughters_          = newtree->Branch("prompt_mujet_numberOfDaughters", &prompt_mujet_numberOfDaughters_, "prompt_mujet_numberOfDaughters/F");
 TBranch* branch_prompt_mujet_muonEnergy_                 = newtree->Branch("prompt_mujet_muonEnergy",        &prompt_mujet_muonEnergy_, "prompt_mujet_muonEnergy/F");
 TBranch* branch_prompt_mujet_muonEnergyFraction_         = newtree->Branch("prompt_mujet_muonEnergyFraction",&prompt_mujet_muonEnergyFraction_, "prompt_mujet_muonEnergyFraction/F");
 TBranch* branch_prompt_mujet_muonMultiplicity_           = newtree->Branch("prompt_mujet_muonMultiplicity",  &prompt_mujet_muonMultiplicity_, "prompt_mujet_muonMultiplicity/F");
 TBranch* branch_prompt_mujet_neutralEmEnergy_            = newtree->Branch("prompt_mujet_neutralEmEnergy",   &prompt_mujet_neutralEmEnergy_, "prompt_mujet_neutralEmEnergy/F");
 TBranch* branch_prompt_mujet_neutralHadronEnergy_        = newtree->Branch("prompt_mujet_neutralHadronEnergy",&prompt_mujet_neutralHadronEnergy_, "prompt_mujet_neutralHadronEnergy/F");
 TBranch* branch_prompt_mujet_NHadMultiplicity_           = newtree->Branch("prompt_mujet_NHadMultiplicity",   &prompt_mujet_NHadMultiplicity_, "prompt_mujet_NHadMultiplicity/F");
 TBranch* branch_prompt_mujet_NMultiplicity_              = newtree->Branch("prompt_mujet_NMultiplicity",      &prompt_mujet_NMultiplicity_, "prompt_mujet_NMultiplicity/F");
 TBranch* branch_prompt_mujet_chargedMultiplicity_        = newtree->Branch("prompt_mujet_chargedMultiplicity",         &prompt_mujet_chargedMultiplicity_, "prompt_mujet_chargedMultiplicity/F");
 TBranch* branch_prompt_mujet_chargedEmEnergyFraction     = newtree->Branch("prompt_mujet_chargedEmEnergyFraction",     &prompt_mujet_chargedEmEnergyFraction, "prompt_mujet_chargedEmEnergyFraction/F");
 TBranch* branch_prompt_mujet_chargedHadronEnergyFraction = newtree->Branch("prompt_mujet_chargedHadronEnergyFraction", &prompt_mujet_chargedHadronEnergyFraction,"prompt_mujet_chargedHadronEnergyFraction/F");

 TBranch* branch_prompt_mujet_notSmeard_pt                = newtree->Branch("prompt_mujet_notSmeard_pt" ,&prompt_mujet_notSmeard_pt,"prompt_mujet_notSmeard_pt/F");
 TBranch* branch_prompt_mujet_raw_pt                      = newtree->Branch("prompt_mujet_raw_pt",&prompt_mujet_notSmeard_pt,"prompt_mujet_notSmeard_pt/F");
 TBranch* branch_prompt_mujet_CsvV2        =  newtree->Branch("prompt_mujet_CsvV2", &prompt_mujet_CsvV2, "prompt_mujet_CsvV2/F");
 TBranch* branch_prompt_mujet_DeepCsv_udsg =  newtree->Branch("prompt_mujet_DeepCsv_udsg", &prompt_mujet_DeepCsv_udsg, "prompt_mujet_DeepCsv_udsg/F");
 TBranch* branch_prompt_mujet_DeepCsv_b    =  newtree->Branch("prompt_mujet_DeepCsv_b", &prompt_mujet_DeepCsv_b, "prompt_mujet_DeepCsv_b/F");
 TBranch* branch_prompt_mujet_DeepCsv_c    =  newtree->Branch("prompt_mujet_DeepCsv_c", &prompt_mujet_DeepCsv_c,"prompt_mujet_DeepCsv_c/F");
 TBranch* branch_prompt_mujet_DeepCsv_bb   =  newtree->Branch("prompt_mujet_DeepCsv_bb", &prompt_mujet_DeepCsv_bb,"prompt_mujet_DeepCsv_bb/F");
 TBranch* branch_prompt_mujet_HadronFlavor =  newtree->Branch("prompt_mujet_HadronFlavor", &prompt_mujet_HadronFlavor,"prompt_mujet_HadronFlavor/F");
 TBranch* branch_prompt_mujet_RSecMu       =  newtree->Branch("prompt_mujet_RSecMu", &prompt_mujet_RSecMu,"prompt_mujet_RSecMu/F");
 TBranch* branch_prompt_mujet_M=  newtree->Branch("prompt_mujet_M", &prompt_mujet_M,"prompt_mujet_M/F");



 //================================= second jet Variables ==============================================//    

 Float_t secondjet_charge_,          secondjet_et_,    secondjet_pt_, secondjet_eta_, secondjet_NEmEnFraction_,
         secondjet_phi_,             secondjet_theta_, secondjet_en_, secondjet_chEmEn_, secondjet_NHadEnFraction_,
         secondjet_chMuEn_,          secondjet_chMuEnFraction_,      secondjet_numberOfDaughters_,
         secondjet_muonEnergyFraction_, secondjet_muonMultiplicity_,
         secondjet_neutralEmEnergy_,    secondjet_neutralHadronEnergy_,           secondjet_NHadMultiplicity_,
         secondjet_NMultiplicity_,      secondjet_chHadEn_,          secondjet_muonEnergy_, secondjet_chargedMultiplicity_ ,
   secondjet_chargedEmEnergyFraction, secondjet_chargedHadronEnergyFraction, secondjet_notSmeard_pt, secondjet_raw_pt;

 Float_t secondjet_CsvV2, secondjet_DeepCsv_udsg, secondjet_DeepCsv_b, secondjet_DeepCsv_c, secondjet_DeepCsv_bb, secondjet_HadronFlavor;
 Int_t   secondjet_bjet_L,secondjet_bjet_M,secondjet_bjet_T,secondjets_size;

 TBranch* branch_secondjet_charge_   = newtree->Branch("secondjet_charge", &secondjet_charge_, "secondjet_charge/F");
 TBranch* branch_secondjet_et_       = newtree->Branch("secondjet_et", &secondjet_et_, "secondjet_et/F");
 TBranch* branch_secondjet_pt_       = newtree->Branch("secondjet_pt", &secondjet_pt_, "secondjet_pt/F");
 TBranch* branch_secondjet_eta_      = newtree->Branch("secondjet_eta", &secondjet_eta_, "secondjet_eta/F");
 TBranch* branch_secondjet_phi_      = newtree->Branch("secondjet_phi", &secondjet_phi_, "secondjet_phi/F");
 TBranch* branch_secondjet_theta_    = newtree->Branch("secondjet_theta", &secondjet_theta_, "secondjet_theta/F");
 TBranch* branch_secondjet_en_       = newtree->Branch("secondjet_en", &secondjet_en_, "secondjet_en/F");
 TBranch* branch_secondjet_chEmEn_                     = newtree->Branch("secondjet_chEmEn", &secondjet_chEmEn_, "secondjet_chEmEn/F");
 TBranch* branch_secondjet_NEmEnFraction_              = newtree->Branch("secondjet_NEmEnFraction", &secondjet_NEmEnFraction_, "secondjet_NEmEnFraction/F");
 TBranch* branch_secondjet_chHadEn_                    = newtree->Branch("secondjet_chHadEn", &secondjet_chHadEn_, "secondjet_chHadEn/F");
 TBranch* branch_secondjet_NHadEnFraction_             = newtree->Branch("secondjet_NHadEnFraction", &secondjet_NHadEnFraction_, "secondjet_NHadEnFraction/F");
 TBranch* branch_secondjet_chMuEn_                     = newtree->Branch("secondjet_chMuEn", &secondjet_chMuEn_, "secondjet_chMuEn/F");
 TBranch* branch_secondjet_chMuEnFraction_             = newtree->Branch("secondjet_chMuEnFraction", &secondjet_chMuEnFraction_, "secondjet_chMuEnFraction/F");
 TBranch* branch_secondjet_numberOfDaughters_          = newtree->Branch("secondjet_numberOfDaughters", &secondjet_numberOfDaughters_, "secondjet_numberOfDaughters/F");
 TBranch* branch_secondjet_muonEnergy_                 = newtree->Branch("secondjet_muonEnergy", &secondjet_muonEnergy_, "secondjet_muonEnergy/F");
 TBranch* branch_secondjet_muonEnergyFraction_         = newtree->Branch("secondjet_muonEnergyFraction", &secondjet_muonEnergyFraction_, "secondjet_muonEnergyFraction/F");
 TBranch* branch_secondjet_muonMultiplicity_           = newtree->Branch("secondjet_muonMultiplicity", &secondjet_muonMultiplicity_, "secondjet_muonMultiplicity/F");
 TBranch* branch_secondjet_neutralEmEnergy_            = newtree->Branch("secondjet_neutralEmEnergy", &secondjet_neutralEmEnergy_, "secondjet_neutralEmEnergy/F");
 TBranch* branch_secondjet_neutralHadronEnergy_        = newtree->Branch("secondjet_neutralHadronEnergy", &secondjet_neutralHadronEnergy_, "secondjet_neutralHadronEnergy/F");
 TBranch* branch_secondjet_NHadMultiplicity_           = newtree->Branch("secondjet_NHadMultiplicity", &secondjet_NHadMultiplicity_, "secondjet_NHadMultiplicity/F");
 TBranch* branch_secondjet_NMultiplicity_              = newtree->Branch("secondjet_NMultiplicity", &secondjet_NMultiplicity_, "secondjet_NMultiplicity/F");
 TBranch* branch_secondjet_chargedMultiplicity_        = newtree->Branch("secondjet_chargedMultiplicity", &secondjet_chargedMultiplicity_, "secondjet_chargedMultiplicity/F");
 TBranch* branch_secondjet_chargedEmEnergyFraction     = newtree->Branch("secondjet_chargedEmEnergyFraction", &secondjet_chargedEmEnergyFraction, "secondjet_chargedEmEnergyFraction/F");
 TBranch* branch_secondjet_chargedHadronEnergyFraction = newtree->Branch("secondjet_chargedHadronEnergyFraction", &secondjet_chargedHadronEnergyFraction,"secondjet_chargedHadronEnergyFraction/F");

 TBranch* branch_secondjets_size  =  newtree->Branch("secondjets_size",  &secondjets_size, "secondjets_size/I");
 TBranch* branch_secondjet_bjet_L =  newtree->Branch("secondjet_bjet_L", &secondjet_bjet_L, "secondjet_bjet_L/I");
 TBranch* branch_secondjet_bjet_M =  newtree->Branch("secondjet_bjet_M", &secondjet_bjet_M, "secondjet_bjet_M/I");
 TBranch* branch_secondjet_bjet_T =  newtree->Branch("secondjet_bjet_T", &secondjet_bjet_T, "secondjet_bjet_T/I");

 TBranch* branch_secondjet_notSmeard_pt                = newtree->Branch("secondjet_notSmeard_pt" ,&secondjet_notSmeard_pt,"secondjet_notSmeard_pt/F");
 TBranch* branch_secondjet_raw_pt                      = newtree->Branch("secondjet_raw_pt",&secondjet_notSmeard_pt,"secondjet_notSmeard_pt/F");

 TBranch* branch_secondjet_CsvV2        =  newtree->Branch("secondjet_CsvV2", &secondjet_CsvV2, "secondjet_CsvV2/F");
 TBranch* branch_secondjet_DeepCsv_udsg =  newtree->Branch("secondjet_DeepCsv_udsg", &secondjet_DeepCsv_udsg, "secondjet_DeepCsv_udsg/F");
 TBranch* branch_secondjet_DeepCsv_b    =  newtree->Branch("secondjet_DeepCsv_b", &secondjet_DeepCsv_b, "secondjet_DeepCsv_b/F");
 TBranch* branch_secondjet_DeepCsv_c    =  newtree->Branch("secondjet_DeepCsv_c", &secondjet_DeepCsv_c,"secondjet_DeepCsv_c/F");
 TBranch* branch_secondjet_DeepCsv_bb   =  newtree->Branch("secondjet_DeepCsv_bb", &secondjet_DeepCsv_bb,"secondjet_DeepCsv_bb/F");
 TBranch* branch_secondjet_HadronFlavor =  newtree->Branch("secondjet_HadronFlavor", &secondjet_HadronFlavor,"secondjet_HadronFlavor/F");

  
//============================ Missing Energy Variables =========================================//
 Float_t pfMet_et_, pfMet_pt_, pfMet_phi_, pfMet_en_, pfMet_sumEt_, caloMet_phi_,caloMet_pt_, Mass_METL,tmass_METL;

 TBranch* branch_pfMet_et_    = newtree->Branch("pfMet_et", &pfMet_et_, "pfMet_et/F");
 TBranch* branch_pfMet_pt_    = newtree->Branch("pfMet_pt", &pfMet_pt_, "pfMet_pt/F");
 TBranch* branch_pfMet_phi_   = newtree->Branch("pfMet_phi", &pfMet_phi_, "pfMet_phi/F");
 TBranch* branch_pfMet_en_    = newtree->Branch("pfMet_en", &pfMet_en_, "pfMet_en/F");
 TBranch* branch_pfMet_sumEt_ = newtree->Branch("pfMet_sumEt", &pfMet_sumEt_, "pfMet_sumEt/F");
 TBranch* branch_caloMet_pt_  = newtree->Branch("caloMet_pt", &caloMet_pt_, "caloMet_pt/F");
 TBranch* branch_caloMet_phi_ = newtree->Branch("caloMet_phi", &caloMet_phi_, "caloMet_phi/F");
 TBranch* branch_Mass_METL    = newtree->Branch("Mass_METL", &Mass_METL, "Mass_METL/F"); 

 TBranch* branch_tmass_METL    = newtree->Branch("tmass_METL", &tmass_METL, "tmass_METL/F");

 //=========================== systematic ================================================================//

 Float_t mujetPt_JECUp, mujetPt_JECDown, mujetSmearedPt_unUp, mujetSmearedPt_unDown, sys_prefire_weight, sys_prefire_weightup, sys_prefire_weightdown;

 TBranch* branch_mujetPt_JECUp           = newtree->Branch("mujetPt_JECUp",&mujetPt_JECUp,"mujetPt_JECUp/F");
 TBranch* branch_mujetPt_JECDown         = newtree->Branch("mujetPt_JECDown",&mujetPt_JECDown,"mujetPt_JECDown/F");
 TBranch* branch_mujetSmearedPt_unUp     = newtree->Branch("mujetSmearedPt_unUp",&mujetSmearedPt_unUp,"mujetSmearedPt_unUp/F");
 TBranch* branch_mujetSmearedPt_unDown   = newtree->Branch("mujetSmearedPt_unDown",&mujetSmearedPt_unDown,"mujetSmearedPt_unDown/F");
 TBranch* branch_sys_prefire_weight      = newtree->Branch("sys_prefire_weight",&sys_prefire_weight,"sys_prefire_weight/F");
 TBranch* branch_sys_prefire_weightup    = newtree->Branch("sys_prefire_weightup",&sys_prefire_weightup,"sys_prefire_weightup/F");
 TBranch* branch_sys_prefire_weightdown  = newtree->Branch("sys_prefire_weightdown",&sys_prefire_weightdown,"sys_prefire_weightdown/F");

//======================= Start the running over input branches ==========================================//
 for (int i=0;i<oldtree->GetEntries(); i++) {

 //for (int i=0;i<100000; i++) {
    if (i%10000==0) cout<<i<<endl;
    oldtree->GetEntry(i);

    unsigned npt_ = -1;
    for(unsigned i=0; i<npT->size(); ++i){
      h_nTrueInteractions50->Fill(npT->at(i));
      h_nTrueInteractions100->Fill(npT->at(i));
      npt_ = i;
    }

    unsigned wt = -1;
    for(unsigned j=0; j<prefire_weight->size(); ++j){
      wt = j;
    }


    //if (passIsoMu24All==0) continue;  // cut on the trigger!
    if (passIsoMu24 == 0 && passIsoMuTk24 == 0) continue;
    ///if ( passMu3_PFJet40  == 0 && passMu8_TrkIsoVVL  == 0 && passMu17_TrkIsoVVL == 0 && passIsoMu24 == 0 && passIsoMuTk24 == 0) continue;
    
    //cout<<"========================== This the new event ==========================="<<endl;

    Float_t   minPt_prompt = -1000;
    unsigned  FirstMuon = -1;
    int count_prompt = 0;
    for(unsigned i=0; i<mu_isLooseMuon->size(); i++){
      if(mu_isLooseMuon->size() == 1) continue;
      if(mu_isLooseMuon->at(i)  == 0 ) continue;
      if (mu_isTightMuon->at(i)==0.  || deltaBeta->at(i)>isoCut
	  || mu_ptTunePMuonBestTrack->at(i) < 25 
	  || abs(mu_etaTunePMuonBestTrack->at(i)) > 2.4 
	  || mu_absdxyTunePMuonBestTrack->at(i) > 0.005 
	  || mu_absdzTunePMuonBestTrack->at(i) > 0.1) continue;  
      count_prompt++;
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
    sort(mu_ptTunePMuonBestTrack->begin(), mu_ptTunePMuonBestTrack->end() , greater<int>()); 

    for(unsigned i=0; i<mu_ptTunePMuonBestTrack->size(); i++){
      //if(i != 1) continue;
      if(i == FirstMuon || FirstMuon == -1) continue;
      if(mu_isLooseMuon->size() == 1) continue;
      if(mu_isLooseMuon->at(i)  == 0 ) continue;
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
      for(unsigned i=0; i<sv_pt->size(); i++){
	if(sv_hasMuon->at(i) == 0) continue;
	if( sv_pt->at(i) < mu_ptTunePMuonBestTrack->at(SecondMuon)) continue;
	if(sv_pt->at(i) > minPt_sv){
	  minPt_sv = sv_pt->at(i);
	  SecondVertex=i;
	}
      }
    }

/*    
    unsigned  SecondVertex = -1;
    Float_t   min_diff = 1000;
    Float_t minRvtx = 999;
    if(SecondMuon != -1){
      float mu_pT  = mu_ptTunePMuonBestTrack->at(SecondMuon);
      float mu_eta = mu_etaTunePMuonBestTrack->at(SecondMuon);
      float mu_phi = mu_phiTunePMuonBestTrack->at(SecondMuon);
      sort(sv_pt->begin(), sv_pt->end() , greater<int>());
      for(unsigned i=0; i<sv_pt->size(); i++){
	if(sv_hasMuon->at(i) == 0) continue;
	if( sv_pt->at(i) < mu_pT) continue;
	for(unsigned j=0; j<sv_tracks_pt->at(i).size(); j++){
	  float trk_pT = sv_tracks_pt->at(i).at(j) ;
	  float trk_eta = sv_tracks_eta->at(i).at(j) ;
          float trk_phi = sv_tracks_phi->at(i).at(j) ; 
	  double DeltapT = abs(mu_pT-trk_pT);
          //float R = sqrt(((mu_eta-trk_eta)*(mu_eta-trk_eta))+((mu_phi-trk_phi)*(mu_phi-trk_phi)));
	  if(DeltapT < min_diff && DeltapT < 0.1 ){
	  //if( R < minRvtx && R < 0.1) {
	       min_diff = DeltapT;
	      //minRvtx = R;
	    SecondVertex=i;
	  }
	}
      }      
    }

*/
    //first track
    unsigned  sv_FirstTrack = -1;
    Float_t  firstTrack_minPt = -1000;
    Float_t sumdxySig = 0;
    if(SecondVertex != -1 ){
      for(unsigned j=0; j<sv_tracks_pt->at(SecondVertex).size(); j++){
        sumdxySig = +abs(sv_tracks_dxySig->at(SecondVertex).at(j));
        if(sv_tracks_pt->at(SecondVertex).at(j) > firstTrack_minPt){
          firstTrack_minPt =sv_tracks_pt->at(SecondVertex).at(j);
          sv_FirstTrack = j;
        }
      }
    }
    //second track
    unsigned  sv_SecondTrack = -1;
    Float_t  secondTrack_minPt = -1000;
    if(SecondVertex != -1 ){
      for(unsigned j=0; j<sv_tracks_pt->at(SecondVertex).size(); j++){
        if(j == sv_FirstTrack) continue;
        if(sv_tracks_pt->at(SecondVertex).at(j) > secondTrack_minPt){
          secondTrack_minPt =sv_tracks_pt->at(SecondVertex).at(j);
          sv_SecondTrack = j;
        }
      }
    }
    //leading jet                                                                                                                                                                                          
    unsigned  alljet = -1;
    Float_t pt_1stjet = -1000;
    int alljet_count = 0;
    int bjet1= 0;
    int bjet2= 0;
    int bjet3= 0;
    Float_t sum_HT = 0; 
    for(unsigned i=0; i<jet_pt->size(); i++){
      if(jet_pt->at(i) > 20 && fabs(jet_eta->at(i)) <= 3.0) sum_HT += jetSmearedPt->at(i); 
      if(SecondMuon == -1) continue;
      if( (fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.90
	   && jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
          ||
	  ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0
	    && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 )
          ||
	  ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
	    && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2)
          ||
	  (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10)) {
	alljet_count++;
	if(jet_pt->at(i) > pt_1stjet) {
	  pt_1stjet = jet_pt->at(i);
	  alljet = i;
	}
	if(jet_DeepCsv_b->at(i) > 0.2219 ){ bjet1++;}
	if(jet_DeepCsv_b->at(i) > 0.6324 ){ bjet2++;}
	if(jet_DeepCsv_b->at(i) > 0.8958 ){ bjet3++;}
      }
    }
      
    unsigned  mujet = -1;
    Float_t minR = 1000;
    for(unsigned i=0; i<jet_pt->size(); i++){
      if(SecondMuon == -1) continue;
      float mu_eta2 = mu_etaTunePMuonBestTrack->at(SecondMuon);
      float mu_phi2 = mu_phiTunePMuonBestTrack->at(SecondMuon);
      float mu_eta = mu_etaTunePMuonBestTrack->at(FirstMuon);
      float mu_phi = mu_phiTunePMuonBestTrack->at(FirstMuon);
      float R1 = sqrt(((mu_eta-jet_eta->at(i))*(mu_eta-jet_eta->at(i)))+((mu_phi-jet_phi->at(i))*(mu_phi-jet_phi->at(i))));
      float R2 = sqrt(((mu_eta2-jet_eta->at(i))*(mu_eta2-jet_eta->at(i)))+((mu_phi2-jet_phi->at(i))*(mu_phi2-jet_phi->at(i))));
      if(R1 < 0.4 || R2 > 0.7 || jet_pt->at(i) < 20) continue;
      if( (fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.90
           && jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
          ||
          ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0
            && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 )
          ||
          ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
            && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2)
          ||
          (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10)) {
        if(R2 < minR) {
          minR = R2;
          mujet = i;
        }
      }
    }
    
    unsigned  prompt_mujet = -1;
    Float_t prompt_minR = 1000;
    for(unsigned i=0; i<jet_pt->size(); i++){
      if(FirstMuon == -1) continue;
      float mu_eta = mu_etaTunePMuonBestTrack->at(FirstMuon);
      float mu_phi = mu_phiTunePMuonBestTrack->at(FirstMuon);
      float R = sqrt(((mu_eta-jet_eta->at(i))*(mu_eta-jet_eta->at(i)))+((mu_phi-jet_phi->at(i))*(mu_phi-jet_phi->at(i))));
      if(R > 0.4 || jet_pt->at(i) < 20) continue;
      if( (fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.90
           && jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
          ||
          ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0
            && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 )
          ||
          ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
            && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2)
          ||
          (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10)) {
        if(R < prompt_minR) {
          prompt_minR = R;
          prompt_mujet = i;
        }
      }
    }


    //sub-leading jet
    unsigned  firstjet = -1;
    Float_t pt_firstjet = -1000;
    int firstjet_count = 0;
    int bjet11= 0;
    int bjet22= 0;
    int bjet33= 0;
    for(unsigned i=0; i<jet_pt->size(); i++){
      if(SecondMuon == -1) continue;
      float mu_eta = mu_etaTunePMuonBestTrack->at(FirstMuon);
      float mu_phi = mu_phiTunePMuonBestTrack->at(FirstMuon);
      float R = sqrt(((mu_eta-jet_eta->at(i))*(mu_eta-jet_eta->at(i)))+((mu_phi-jet_phi->at(i))*(mu_phi-jet_phi->at(i))));
      if(R < 0.4 || jet_pt->at(i) < 20) continue;
      if( (fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.90
           && jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
          ||
          ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0
            && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 )
          ||
          ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
            && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2)
          ||
          (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10)) {
	firstjet_count++;
        if(jet_pt->at(i) > pt_firstjet) {
          pt_firstjet = jet_pt->at(i);
          firstjet = i;
        }
	if(jet_DeepCsv_b->at(i) > 0.2219 ){ bjet11++;}
        if(jet_DeepCsv_b->at(i) > 0.6324 ){ bjet22++;}
        if(jet_DeepCsv_b->at(i) > 0.8958 ){ bjet33++;}
      }
    }
    
    unsigned  secondjet = -1;
    Float_t pt_secondjet = -1000;
    int secondjet_count = 0;
    int bjet111= 0;
    int bjet222= 0;
    int bjet333= 0;
    for(unsigned i=0; i<jet_pt->size(); i++){
      if(SecondMuon == -1) continue;
      float mu1_eta = mu_etaTunePMuonBestTrack->at(FirstMuon);
      float mu1_phi = mu_phiTunePMuonBestTrack->at(FirstMuon);
      float mu2_eta = mu_etaTunePMuonBestTrack->at(SecondMuon);
      float mu2_phi = mu_phiTunePMuonBestTrack->at(SecondMuon);
      float R1 = sqrt(((mu1_eta-jet_eta->at(i))*(mu1_eta-jet_eta->at(i)))+((mu1_phi-jet_phi->at(i))*(mu1_phi-jet_phi->at(i))));
      float R2 = sqrt(((mu2_eta-jet_eta->at(i))*(mu2_eta-jet_eta->at(i)))+((mu2_phi-jet_phi->at(i))*(mu2_phi-jet_phi->at(i))));
      if(R1 < 0.4 || R2 < 0.7 || jet_pt->at(i) < 20) continue;
      if(i == mujet) continue;
      if( (fabs(jet_eta->at(i)) <= 2.7 && jet_neutralHadronEnergyFraction->at(i) <=  0.90
           && jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_chargedMultiplicity->at(i)+jet_neutralMultiplicity->at(i) >= 1 )
          ||
          ( fabs(jet_eta->at(i)) <= 2.4 && jet_chargedHadronEnergyFraction->at(i) >= 0
            && jet_chargedMultiplicity->at(i) >= 0 && jet_chargedEmEnergyFraction->at(i) <= 0.99 )
          ||
          ( fabs(jet_eta->at(i)) <= 3.0 && jet_neutralHadronEnergyFraction->at(i) <= 0.98
            && jet_neutralEmEnergyFraction->at(i) >= 0.01 && jet_neutralMultiplicity->at(i) >= 2)
          ||
          (jet_neutralEmEnergyFraction->at(i) <= 0.90 && jet_neutralMultiplicity->at(i) >= 10)) {
	secondjet_count++;
        if(jet_pt->at(i) > pt_secondjet) {
          pt_secondjet = jet_pt->at(i);
          secondjet = i;
        }
        if(jet_DeepCsv_b->at(i) > 0.2219 ){ bjet111++;}
        if(jet_DeepCsv_b->at(i) > 0.6324 ){ bjet222++;}
        if(jet_DeepCsv_b->at(i) > 0.8958 ){ bjet333++;}
      }
    }

   
    if(FirstMuon != -1 && SecondMuon != -1 && SecondVertex != -1 && mujet != -1){
    //if(FirstMuon != -1 && SecondMuon != -1 && SecondVertex != -1 ){
                         //GenCtau = (lep2_gen_MomCTau0 > -1) ? lep2_gen_MomCTau0 : 0;
      //pile up weight
      TrueNumInteractions = (isMC) ? npT->at(npt_) : 0;
      mc_gen_weight = (isMC) ? gen_weight  : 0;
      mc_lhe_weight = (isMC) ? lhe_weight  : 0;
      mc_lhe_ctau   = (isMC) ? lhe_ctau    : 0;

      sys_prefire_weight     = (isMC) ? prefire_weight->at(wt)     :0;
      sys_prefire_weightup   = (isMC) ? prefire_weightup->at(wt)   :0;
      sys_prefire_weightdown = (isMC) ? prefire_weightdown->at(wt) :0;


      //pv info     
      pvX_        =  pvX ;
      pvY_        =  pvY ;
      pvZ_        =  pvZ ;
      pvXErr_     =  pvXErr ;
      pvYErr_     =  pvYErr ;
      pvZErr_     =  pvZErr ;
      pvLxy_      =  pvLxy  ;
      pvLxyz_     =  pvLxyz ;
      pvLxySig_   =  pvLxySigma  ;
      pvLxyzSig_  =  pvLxyzSigma ;
      pvChi2_     =  pvChi2 ;
      pvSumPtSq_  =  pvSumPtSq ;
      //numberPV_   =  numberPV;
      //pvNTrack_   =  pvNTrack ;

  
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
   mu_promptAbsdzError     = mu_absdzErrorTunePMuonBestTrack->at(FirstMuon);
   mu_promptAbsdzSig       = (mu_promptAbsdz/mu_promptAbsdzError) ;//mu_absdzSigTunePMuonBestTrack->at(FirstMuon);  
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
   mu_promptsize                   =  count_prompt;

   mu_promptMotherID  = (isMC) ? mu_DecayChain->at(FirstMuon): 0;

   mu_prompt3dIP     = mu_3dIP->at(FirstMuon);
   mu_prompt3dIPSig  = mu_3dIPSig->at(FirstMuon);
   mu_prompt2dIP     = mu_2dIP->at(FirstMuon);
   mu_prompt2dIPSig  = mu_2dIPSig->at(FirstMuon);

   double A_eff_prompt;

   if     (abs(mu_promptEta)<=.8)                      { A_eff_prompt=0.0735; }
   else if(abs(mu_promptEta)>0.8 && abs(mu_promptEta)<=1.3) { A_eff_prompt=0.0619; }
   else if(abs(mu_promptEta)>1.3 && abs(mu_promptEta)<=2.0) { A_eff_prompt=0.0465; }
   else if(abs(mu_promptEta)>2.0 && abs(mu_promptEta)<=2.2) { A_eff_prompt=0.0433; }
   else if(abs(mu_promptEta)>2.2 && abs(mu_promptEta)<=2.4) { A_eff_prompt=0.0577; }
   float charged_prompt =  mu_pfSumChargedHadronPt->at(FirstMuon);
   float neutral_prompt =  mu_pfSumNeutralHadronEt->at(FirstMuon);
   float sumPhotonEt_prompt = mu_PFSumPhotonEt->at(FirstMuon);
   float pileup_prompt =  mu_pfSumPUPt->at(FirstMuon);

   RelIso_prompt = (charged_prompt + std::max(0.0, neutral_prompt+sumPhotonEt_prompt-  mu_promptRhoIso*A_eff_prompt))/mu_promptPt;


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

    mu_secondInnerTrackFraction   = mu_InnerTrackValidFraction->at(SecondMuon);
    mu_secondSegmentCompatibility = mu_segmentCompatibilityMuonBestTrack->at(SecondMuon);
    mu_secondTrkKink              = mu_trkKinkMuonBestTrack->at(SecondMuon);
    mu_secondChi2LocalPosition    = mu_chi2LocalPositionMuonBestTrack->at(SecondMuon);
    mu_secondGlobalMuon           = mu_isGlobalMuon->at(SecondMuon);

    mu_secondRPCTofDirection        = mu_RPCTofDirection->at(SecondMuon);
    mu_secondRPCTofNDof             = mu_RPCTofNDof->at(SecondMuon);
    mu_secondRPCTofTimeAtIpInOut    = mu_RPCTofTimeAtIpInOut->at(SecondMuon);
    mu_secondRPCTofTimeAtIpInOutErr = mu_RPCTofTimeAtIpInOutErr->at(SecondMuon);
    mu_secondRPCTofTimeAtIpOutIn    = mu_RPCTofTimeAtIpOutIn->at(SecondMuon);
    mu_secondRPCTofTimeAtIpOutInErr = mu_RPCTofTimeAtIpOutInErr->at(SecondMuon);

    double A_eff_second;

    if     (abs(mu_secondEta)<=.8)                      { A_eff_second=0.0735; }
    else if(abs(mu_secondEta)>0.8 && abs(mu_secondEta)<=1.3) { A_eff_second=0.0619; }
    else if(abs(mu_secondEta)>1.3 && abs(mu_secondEta)<=2.0) { A_eff_second=0.0465; }
    else if(abs(mu_secondEta)>2.0 && abs(mu_secondEta)<=2.2) { A_eff_second=0.0433; }
    else if(abs(mu_secondEta)>2.2 && abs(mu_secondEta)<=2.4) { A_eff_second=0.0577; }

    float R = sqrt(((mu_secondEta-mu_promptEta)*(mu_secondEta-mu_promptEta))+((mu_secondPhi-mu_promptPhi)*(mu_secondPhi-mu_promptPhi)));
    float charged =  mu_pfSumChargedHadronPt->at(SecondMuon);
    float neutral =  mu_pfSumNeutralHadronEt->at(SecondMuon);
    float sumPhotonEt = mu_PFSumPhotonEt->at(SecondMuon);
    float pileup =  mu_pfSumPUPt->at(SecondMuon); 

    double deltaBetaR3 = (charged + std::max(0.0, neutral+sumPhotonEt-0.5*pileup))/mu_secondPt;
    RelIso_second = (charged + std::max(0.0, neutral+sumPhotonEt- mu_secondRhoIso*A_eff_second))/mu_secondPt;

    mu_secondMotherID  = (isMC) ? mu_DecayChain->at(SecondMuon): 0;

    mu_second3dIP     = mu_3dIP->at(SecondMuon);
    mu_second3dIPSig  = mu_3dIPSig->at(SecondMuon);
    mu_second2dIP     = mu_2dIP->at(SecondMuon);
    mu_second2dIPSig  = mu_2dIPSig->at(SecondMuon);

    TLorentzVector Mu1;
    TLorentzVector Mu2; 
    TLorentzVector Jet1;
    TLorentzVector Jet2;
    //Mu1.SetPtEtaPhiE(mu_promptPt,mu_promptEta,mu_promptPhi,mu_promptE);
    //Mu2.SetPtEtaPhiE(mu_secondPt,mu_secondEta,mu_secondPhi,mu_secondE);

    float px1 = mu_pxTunePMuonBestTrack->at(FirstMuon);
    float px2 = mu_pxTunePMuonBestTrack->at(SecondMuon);
    float py1 = mu_pyTunePMuonBestTrack->at(FirstMuon);
    float py2 = mu_pyTunePMuonBestTrack->at(SecondMuon);
    float pz1 = mu_pzTunePMuonBestTrack->at(FirstMuon);
    float pz2 = mu_pzTunePMuonBestTrack->at(SecondMuon);
    float mT  = sqrt(pow(mu_promptPt + mu_secondPt, 2) - pow(px1 + px2, 2) - pow(py1 + py2, 2));

    Mu1.SetPxPyPzE(px1,py1,pz1,mu_promptE);
    Mu2.SetPxPyPzE(px2,py2,pz2,mu_secondE);

    float DiMuMass = (Mu1 + Mu2).M();
    float DeltaR = Mu1.DeltaR(Mu2);
    float DeltaPhi = abs(Mu1.DeltaPhi(Mu2));

    mu_DeltaBetaR3 = deltaBetaR3;
    mu_DiMuMass    = DiMuMass;
    mu_Size        = count;
    mu_DeltaPhi    = DeltaPhi ;
    mu_DeltaR_vec  = DeltaR ;
    mu_DeltaR      = R ;

    mu_mT          = mT ;
    mu_nbLoose     = nb_loose;

   //secondary vertex info   

      sv_TrackSize_  =  (SecondVertex != -1) ? sv_munTracks->at(SecondVertex) : -999;
      sv_LXYSig_     =  (SecondVertex != -1) ? sv_LxySig->at(SecondVertex)    : -999;
      sv_LXYZSig_    =  (SecondVertex != -1) ? sv_LxyzSig->at(SecondVertex)   : -999;
      sv_LXY_        =  (SecondVertex != -1) ? sv_Lxy->at(SecondVertex)       : -999;
      sv_LXYZ_       =  (SecondVertex != -1) ? sv_Lxyz->at(SecondVertex)      : -999;
      sv_mass_       =  (SecondVertex != -1) ? sv_mass->at(SecondVertex)      : -999;
      sv_eta_        =  (SecondVertex != -1) ? sv_eta->at(SecondVertex)       : -999;
      sv_phi_        =  (SecondVertex != -1) ? sv_phi->at(SecondVertex)       : -999;
      sv_pt_         =  (SecondVertex != -1) ? sv_pt->at(SecondVertex)        : -999;
      sv_p_          =  (SecondVertex != -1) ? sv_p->at(SecondVertex)         : -999;
      sv_px_         =  (SecondVertex != -1) ? sv_px->at(SecondVertex)        : -999;
      sv_py_         =  (SecondVertex != -1) ? sv_py->at(SecondVertex)        : -999;
      sv_pz_         =  (SecondVertex != -1) ? sv_pz->at(SecondVertex)        : -999;
      sv_energy_     =  (SecondVertex != -1) ? sv_energy->at(SecondVertex)    : -999;
      sv_Beta_       =  (SecondVertex != -1) ? sv_Beta->at(SecondVertex)      : -999;
      sv_Gamma_      =  (SecondVertex != -1) ? sv_Gamma->at(SecondVertex)     : -999;
      sv_CTau0_      =  (SecondVertex != -1) ? sv_CTau0->at(SecondVertex)     : -999;
      sv_NDof_       =  (SecondVertex != -1) ? sv_NDof->at(SecondVertex)      : -999;
      sv_Chi2_       =  (SecondVertex != -1) ? sv_Chi2->at(SecondVertex)      : -999;
      sv_Angle3D_    =  (SecondVertex != -1) ?   sv_Angle3D->at(SecondVertex) : -999;
      sv_Angle2D_    =  (SecondVertex != -1) ? sv_Angle2D->at(SecondVertex)   : -999;
      sv_tracks_Sumcharge_ = (SecondVertex != -1) ? sv_tracks_Sumcharge->at(SecondVertex): -999;
      sv_tracks_Sumpt_     = (SecondVertex != -1) ? sv_tracks_Sumpt->at(SecondVertex)    : -999;
      sv_match_            = (SecondVertex != -1) ? sv_match_dxyz->at(SecondVertex)      : -999;
      sv_track_sumdxySig_  = (SecondVertex != -1) ? sumdxySig                            : -999;
      sv_hasMuon_          = (SecondVertex != -1) ? sv_hasMuon->at(SecondVertex)         : -999;

      TLorentzVector vtx;

      //vtx.SetPtEtaPhiE(sv_pt_,sv_eta_,sv_phi_,sv_energy_);
      vtx.SetPxPyPzE(sv_px_, sv_py_, sv_pz_,sv_energy_);
      //float vtxmumass =  (SecondVertex != -1) ? (ele1 + vtx).M() : -999 ;

      float vtxmumass =  (Mu1 + vtx).M();

      //ele_Size         = count;
      vtxmu_mass       = vtxmumass;

      sv_Xpos_     = (SecondVertex != -1) ? sv_X->at(SecondVertex)  : -999;
      sv_Ypos_     = (SecondVertex != -1) ? sv_Y->at(SecondVertex)  : -999;
      sv_Zpos_     = (SecondVertex != -1) ? sv_Z->at(SecondVertex)  : -999;
      sv_xError_   = (SecondVertex != -1) ? sv_xErr->at(SecondVertex)  : -999;
      sv_yError_   = (SecondVertex != -1) ? sv_yErr->at(SecondVertex)  : -999;
      sv_zError_   = (SecondVertex != -1) ? sv_zErr->at(SecondVertex)  : -999;
      sv_lx_ = (SecondVertex != -1) ? sv_lx->at(SecondVertex)  : -999;
      sv_ly_ = (SecondVertex != -1) ? sv_ly->at(SecondVertex)  : -999;
      sv_lx_ = (SecondVertex != -1) ? sv_lz->at(SecondVertex)  : -999;

      //first track in sv info
      firstTrack_eta    =  (SecondVertex != -1) ? sv_tracks_eta->at(SecondVertex).at(sv_FirstTrack)    : -999;
      firstTrack_phi    =  (SecondVertex != -1) ? sv_tracks_phi->at(SecondVertex).at(sv_FirstTrack)    : -999;
      firstTrack_pt     =  (SecondVertex != -1) ? sv_tracks_pt->at(SecondVertex).at(sv_FirstTrack)     : -999;
      firstTrack_dxySig =  (SecondVertex != -1) ? sv_tracks_dxySig->at(SecondVertex).at(sv_FirstTrack) : -999;
      firstTrack_dxy    =  (SecondVertex != -1) ? sv_tracks_dxy->at(SecondVertex).at(sv_FirstTrack)    : -999;
      firstTrack_dxyz   =  (SecondVertex != -1) ? sv_tracks_dxyz->at(SecondVertex).at(sv_FirstTrack)   : -999;
      firstTrack_charge =  (SecondVertex != -1) ? sv_tracks_charge->at(SecondVertex).at(sv_FirstTrack) : -999;
      firstTrack_en     =  (SecondVertex != -1) ? sv_tracks_p->at(SecondVertex).at(sv_FirstTrack)      : -999;

      //second track in sv info
      secondTrack_eta    =  (SecondVertex != -1) ? sv_tracks_eta->at(SecondVertex).at(sv_SecondTrack)   : -999;
      secondTrack_phi    =  (SecondVertex != -1) ? sv_tracks_phi->at(SecondVertex).at(sv_SecondTrack)   : -999;
      secondTrack_pt     =  (SecondVertex != -1) ? sv_tracks_pt->at(SecondVertex).at(sv_SecondTrack)    : -999;
      secondTrack_dxySig =  (SecondVertex != -1) ? sv_tracks_dxySig->at(SecondVertex).at(sv_SecondTrack): -999;
      secondTrack_dxy    =  (SecondVertex != -1) ? sv_tracks_dxy->at(SecondVertex).at(sv_SecondTrack)   : -999;
      secondTrack_dxyz   =  (SecondVertex != -1) ? sv_tracks_dxyz->at(SecondVertex).at(sv_SecondTrack)  : -999;
      secondTrack_charge =  (SecondVertex != -1) ? sv_tracks_charge->at(SecondVertex).at(sv_SecondTrack): -999;
      secondTrack_en     =  (SecondVertex != -1) ? sv_tracks_p->at(SecondVertex).at(sv_SecondTrack)     : -999;


      sv_tracks_charge_n.clear();
      sv_tracks_eta_n.clear();
      sv_tracks_phi_n.clear();
      sv_tracks_pt_n.clear();
      sv_tracks_p_n.clear();
      sv_tracks_dxySig_n.clear();
      sv_tracks_dxy_n.clear();
      sv_tracks_dxyz_n.clear();
      sv_mass_n.clear(); 
      sv_eta_n.clear();
      sv_phi_n.clear(); 
      sv_pt_n.clear(); 
      sv_p_n.clear(); 
      sv_px_n.clear(); 
      sv_py_n.clear(); 
      sv_pz_n.clear(); 
      sv_energy_n.clear(); 
      sv_Beta_n.clear(); 
      sv_Gamma_n.clear(); 
      sv_CTau0_n.clear(); 
      sv_NDof_n.clear(); 
      sv_Chi2_n.clear(); 
      sv_Angle3D_n.clear(); 
      sv_Angle2D_n.clear();

      for(unsigned j=0; j<sv_tracks_pt->at(SecondVertex).size(); j++){
	sv_mass_n   .push_back( sv_mass->at(SecondVertex)     );
	sv_eta_n    .push_back( sv_eta->at(SecondVertex)      );
	sv_phi_n    .push_back( sv_phi->at(SecondVertex)      );
	sv_pt_n     .push_back( sv_pt->at(SecondVertex)       );
	sv_p_n      .push_back( sv_p->at(SecondVertex)        );
	sv_px_n     .push_back( sv_px->at(SecondVertex)       );
	sv_py_n     .push_back( sv_py->at(SecondVertex)       );
	sv_pz_n     .push_back( sv_pz->at(SecondVertex)       );
	sv_energy_n .push_back( sv_energy->at(SecondVertex)   );
	sv_Beta_n   .push_back( sv_Beta->at(SecondVertex)     );
	sv_Gamma_n  .push_back( sv_Gamma->at(SecondVertex)    );
	sv_CTau0_n  .push_back( sv_CTau0->at(SecondVertex)    );
	sv_NDof_n   .push_back( sv_NDof->at(SecondVertex)     );
	sv_Chi2_n   .push_back( sv_Chi2->at(SecondVertex)     );
	sv_Angle3D_n.push_back( sv_Angle3D->at(SecondVertex)  );
	sv_Angle2D_n.push_back( sv_Angle2D->at(SecondVertex)  );
	sv_tracks_charge_n.push_back(sv_tracks_charge->at(SecondVertex).at(j));
	sv_tracks_eta_n.push_back(sv_tracks_eta->at(SecondVertex).at(j));
	sv_tracks_phi_n.push_back(sv_tracks_phi->at(SecondVertex).at(j));
	sv_tracks_pt_n.push_back(sv_tracks_pt->at(SecondVertex).at(j));
	sv_tracks_p_n.push_back(sv_tracks_p->at(SecondVertex).at(j));
	sv_tracks_dxySig_n.push_back(sv_tracks_dxySig->at(SecondVertex).at(j));
	sv_tracks_dxy_n.push_back( sv_tracks_dxy->at(SecondVertex).at(j));	
        sv_tracks_dxyz_n.push_back( sv_tracks_dxyz->at(SecondVertex).at(j));

      }


      //leading jet whatever info
    
      alljet_charge_                     = (alljet != -1) ? jet_charge->at(alljet)                      : -999;
      alljet_et_                         = (alljet != -1) ? jet_et->at(alljet)                          : -999;
      alljet_pt_                         = (alljet != -1) ? jetSmearedPt->at(alljet)                    : -999;
      alljet_eta_                        = (alljet != -1) ? jet_eta->at(alljet)                         : -999;
      alljet_phi_                        = (alljet != -1) ? jet_phi->at(alljet)                         : -999;
      alljet_theta_                      = (alljet != -1) ? jet_theta->at(alljet)                       : -999;
      alljet_en_                         = (alljet != -1) ? jet_en->at(alljet)                          : -999;
      alljet_chEmEn_                     = (alljet != -1) ? jet_chargedEmEnergy->at(alljet)             : -999;
      alljet_NEmEnFraction_              = (alljet != -1) ? jet_neutralEmEnergyFraction->at(alljet)     : -999;
      alljet_chHadEn_                    = (alljet != -1) ? jet_chargedHadronEnergy->at(alljet)         : -999;
      alljet_NHadEnFraction_             = (alljet != -1) ? jet_neutralHadronEnergyFraction->at(alljet) : -999;
      alljet_chMuEn_                     = (alljet != -1) ? jet_chargedMuEnergy->at(alljet)             : -999;
      alljet_chMuEnFraction_             = (alljet != -1) ? jet_chargedMuEnergyFraction->at(alljet)     : -999;
      alljet_numberOfDaughters_          = (alljet != -1) ? jet_numberOfDaughters->at(alljet)           : -999;
      alljet_muonEnergy_                 = (alljet != -1) ? jet_muonEnergy->at(alljet)                  : -999;
      alljet_muonEnergyFraction_         = (alljet != -1) ? jet_muonEnergyFraction->at(alljet)          : -999;
      alljet_muonMultiplicity_           = (alljet != -1) ? jet_muonMultiplicity->at(alljet)            : -999;
      alljet_neutralEmEnergy_            = (alljet != -1) ? jet_neutralEmEnergy->at(alljet)             : -999;
      alljet_neutralHadronEnergy_        = (alljet != -1) ? jet_neutralHadronEnergy->at(alljet)         : -999;
      alljet_NHadMultiplicity_           = (alljet != -1) ? jet_neutralHadronMultiplicity->at(alljet)   : -999;
      alljet_NMultiplicity_              = (alljet != -1) ? jet_neutralMultiplicity->at(alljet)         : -999;
      alljet_chargedMultiplicity_        = (alljet != -1) ? jet_chargedMultiplicity->at(alljet)         : -999;
      alljet_chargedEmEnergyFraction     = (alljet != -1) ? jet_chargedEmEnergyFraction->at(alljet)     : -999;
      alljet_chargedHadronEnergyFraction = (alljet != -1) ? jet_chargedHadronEnergyFraction->at(alljet) : -999;

      alljet_notSmeard_pt                = (alljet != -1) ? jet_pt->at(alljet) : -999;
      alljet_raw_pt                      = (alljet != -1) ? jet_ptuncorrected->at(alljet) : -999;

      alljet_CsvV2          = (alljet != -1) ? jet_CsvV2->at(alljet)        : -999;
      alljet_DeepCsv_udsg   = (alljet != -1) ? jet_DeepCsv_udsg->at(alljet) : -999;
      alljet_DeepCsv_b      = (alljet != -1) ? jet_DeepCsv_b->at(alljet)    : -999;
      alljet_DeepCsv_c      = (alljet != -1) ? jet_DeepCsv_c->at(alljet)    : -999;
      alljet_DeepCsv_bb     = (alljet != -1) ? jet_DeepCsv_bb->at(alljet)   : -999;
      alljet_HadronFlavor   = (alljet != -1) ? jet_HadronFlavor->at(alljet) : -999;

      alljets_size = alljet_count;
      alljet_bjet_L = (alljet != -1) ? bjet1 : -999;
      alljet_bjet_M = (alljet != -1) ? bjet2 : -999;
      alljet_bjet_T = (alljet != -1) ? bjet3 : -999;

      HT_alljet =  sum_HT;

      //leading jet whithout first jet info  

      firstjet_charge_                     = (firstjet != -1) ? jet_charge->at(firstjet)                      : -999;
      firstjet_et_                         = (firstjet != -1) ? jet_et->at(firstjet)                          : -999;
      firstjet_pt_                         = (firstjet != -1) ? jetSmearedPt->at(firstjet)                    : -999;
      firstjet_eta_                        = (firstjet != -1) ? jet_eta->at(firstjet)                         : -999;
      firstjet_phi_                        = (firstjet != -1) ? jet_phi->at(firstjet)                         : -999;
      firstjet_theta_                      = (firstjet != -1) ? jet_theta->at(firstjet)                       : -999;
      firstjet_en_                         = (firstjet != -1) ? jet_en->at(firstjet)                          : -999;
      firstjet_chEmEn_                     = (firstjet != -1) ? jet_chargedEmEnergy->at(firstjet)             : -999;
      firstjet_NEmEnFraction_              = (firstjet != -1) ? jet_neutralEmEnergyFraction->at(firstjet)     : -999;
      firstjet_chHadEn_                    = (firstjet != -1) ? jet_chargedHadronEnergy->at(firstjet)         : -999;
      firstjet_NHadEnFraction_             = (firstjet != -1) ? jet_neutralHadronEnergyFraction->at(firstjet) : -999;
      firstjet_chMuEn_                     = (firstjet != -1) ? jet_chargedMuEnergy->at(firstjet)             : -999;
      firstjet_chMuEnFraction_             = (firstjet != -1) ? jet_chargedMuEnergyFraction->at(firstjet)     : -999;
      firstjet_numberOfDaughters_          = (firstjet != -1) ? jet_numberOfDaughters->at(firstjet)           : -999;
      firstjet_muonEnergy_                 = (firstjet != -1) ? jet_muonEnergy->at(firstjet)                  : -999;
      firstjet_muonEnergyFraction_         = (firstjet != -1) ? jet_muonEnergyFraction->at(firstjet)          : -999;
      firstjet_muonMultiplicity_           = (firstjet != -1) ? jet_muonMultiplicity->at(firstjet)            : -999;
      firstjet_neutralEmEnergy_            = (firstjet != -1) ? jet_neutralEmEnergy->at(firstjet)             : -999;
      firstjet_neutralHadronEnergy_        = (firstjet != -1) ? jet_neutralHadronEnergy->at(firstjet)         : -999;
      firstjet_NHadMultiplicity_           = (firstjet != -1) ? jet_neutralHadronMultiplicity->at(firstjet)   : -999;
      firstjet_NMultiplicity_              = (firstjet != -1) ? jet_neutralMultiplicity->at(firstjet)         : -999;
      firstjet_chargedMultiplicity_        = (firstjet != -1) ? jet_chargedMultiplicity->at(firstjet)         : -999;
      firstjet_chargedEmEnergyFraction     = (firstjet != -1) ? jet_chargedEmEnergyFraction->at(firstjet)     : -999;
      firstjet_chargedHadronEnergyFraction = (firstjet != -1) ? jet_chargedHadronEnergyFraction->at(firstjet) : -999;

      firstjet_CsvV2          = (firstjet != -1) ? jet_CsvV2->at(firstjet)        : -999;
      firstjet_DeepCsv_udsg   = (firstjet != -1) ? jet_DeepCsv_udsg->at(firstjet) : -999;
      firstjet_DeepCsv_b      = (firstjet != -1) ? jet_DeepCsv_b->at(firstjet)    : -999;
      firstjet_DeepCsv_c      = (firstjet != -1) ? jet_DeepCsv_c->at(firstjet)    : -999;
      firstjet_DeepCsv_bb     = (firstjet != -1) ? jet_DeepCsv_bb->at(firstjet)   : -999;
      firstjet_HadronFlavor   = (firstjet != -1) ? jet_HadronFlavor->at(firstjet) : -999;

      firstjet_notSmeard_pt                = (firstjet != -1) ? jet_pt->at(firstjet) : -999;
      firstjet_raw_pt                      = (firstjet != -1) ? jet_ptuncorrected->at(firstjet) : -999;

      firstjets_size  = (firstjet != -1) ? firstjet_count : -999;
      firstjet_bjet_L = (firstjet != -1) ? bjet11         : -999;
      firstjet_bjet_M = (firstjet != -1) ? bjet22         : -999;
      firstjet_bjet_T = (firstjet != -1) ? bjet33         : -999;

      //jet with in ele info
      mujet_charge_                     = (mujet != -1) ? jet_charge->at(mujet)                      : -999;
      mujet_et_                         = (mujet != -1) ? jet_et->at(mujet)                          : -999;
      mujet_pt_                         = (mujet != -1) ? jetSmearedPt->at(mujet)                    : -999;
      mujet_eta_                        = (mujet != -1) ? jet_eta->at(mujet)                         : -999;
      mujet_phi_                        = (mujet != -1) ? jet_phi->at(mujet)                         : -999;
      mujet_theta_                      = (mujet != -1) ? jet_theta->at(mujet)                       : -999;
      mujet_en_                         = (mujet != -1) ? jet_en->at(mujet)                          : -999;
      mujet_chEmEn_                     = (mujet != -1) ? jet_chargedEmEnergy->at(mujet)             : -999;
      mujet_NEmEnFraction_              = (mujet != -1) ? jet_neutralEmEnergyFraction->at(mujet)     : -999;
      mujet_chHadEn_                    = (mujet != -1) ? jet_chargedHadronEnergy->at(mujet)         : -999;
      mujet_NHadEnFraction_             = (mujet != -1) ? jet_neutralHadronEnergyFraction->at(mujet) : -999;
      mujet_chMuEn_                     = (mujet != -1) ? jet_chargedMuEnergy->at(mujet)             : -999;
      mujet_chMuEnFraction_             = (mujet != -1) ? jet_chargedMuEnergyFraction->at(mujet)     : -999;
      mujet_numberOfDaughters_          = (mujet != -1) ? jet_numberOfDaughters->at(mujet)           : -999;
      mujet_muonEnergy_                 = (mujet != -1) ? jet_muonEnergy->at(mujet)                  : -999;
      mujet_muonEnergyFraction_         = (mujet != -1) ? jet_muonEnergyFraction->at(mujet)          : -999;
      mujet_muonMultiplicity_           = (mujet != -1) ? jet_muonMultiplicity->at(mujet)            : -999;
      mujet_neutralEmEnergy_            = (mujet != -1) ? jet_neutralEmEnergy->at(mujet)             : -999;
      mujet_neutralHadronEnergy_        = (mujet != -1) ? jet_neutralHadronEnergy->at(mujet)         : -999;
      mujet_NHadMultiplicity_           = (mujet != -1) ? jet_neutralHadronMultiplicity->at(mujet)   : -999;
      mujet_NMultiplicity_              = (mujet != -1) ? jet_neutralMultiplicity->at(mujet)         : -999;
      mujet_chargedMultiplicity_        = (mujet != -1) ? jet_chargedMultiplicity->at(mujet)         : -999;
      mujet_chargedEmEnergyFraction     = (mujet != -1) ? jet_chargedEmEnergyFraction->at(mujet)     : -999;
      mujet_chargedHadronEnergyFraction = (mujet != -1) ? jet_chargedHadronEnergyFraction->at(mujet) : -999;

      mujet_notSmeard_pt                = (mujet != -1) ? jet_pt->at(mujet) : -999;
      mujet_raw_pt                      = (mujet != -1) ? jet_ptuncorrected->at(mujet) : -999;

      mujet_CsvV2          = (mujet != -1) ? jet_CsvV2->at(mujet)        : -999;
      mujet_DeepCsv_udsg   = (mujet != -1) ? jet_DeepCsv_udsg->at(mujet) : -999;
      mujet_DeepCsv_b      = (mujet != -1) ? jet_DeepCsv_b->at(mujet)    : -999;
      mujet_DeepCsv_c      = (mujet != -1) ? jet_DeepCsv_c->at(mujet)    : -999;
      mujet_DeepCsv_bb     = (mujet != -1) ? jet_DeepCsv_bb->at(mujet)   : -999;
      mujet_HadronFlavor   = (mujet != -1) ? jet_HadronFlavor->at(mujet) : -999;
      mujet_RSecMu        = (mujet != -1) ? minR : -999;

      TLorentzVector mu_jet;
      mu_jet.SetPtEtaPhiE(mujet_pt_,mujet_eta_,mujet_phi_,mujet_en_);
      float mujetmass =  (mujet != -1) ? (Mu1 + mu_jet).M() : -999 ;
      float hnl_m = (mujet != -1) ?  (mu_jet).M() : -999;
      hnl_mass = hnl_m ;
      mujet_M = mujetmass;

      mujetPt_JECUp         =  (mujet != -1) ? jetPt_JECUp->at(mujet)          : -999;
      mujetPt_JECDown       =  (mujet != -1) ? jetPt_JECDown->at(mujet)        : -999;
      mujetSmearedPt_unUp   =  (mujet != -1) ? jetSmearedPt_unUp->at(mujet)    : -999;
      mujetSmearedPt_unDown =  (mujet != -1) ? jetSmearedPt_unDown->at(mujet)  : -999;

      // jets closer to prompt mu 
      prompt_mujet_charge_                     = (prompt_mujet != -1) ? jet_charge->at(prompt_mujet)                      : -999;
      prompt_mujet_et_                         = (prompt_mujet != -1) ? jet_et->at(prompt_mujet)                          : -999;
      prompt_mujet_pt_                         = (prompt_mujet != -1) ? jetSmearedPt->at(prompt_mujet)                    : -999;
      prompt_mujet_eta_                        = (prompt_mujet != -1) ? jet_eta->at(prompt_mujet)                         : -999;
      prompt_mujet_phi_                        = (prompt_mujet != -1) ? jet_phi->at(prompt_mujet)                         : -999;
      prompt_mujet_theta_                      = (prompt_mujet != -1) ? jet_theta->at(prompt_mujet)                       : -999;
      prompt_mujet_en_                         = (prompt_mujet != -1) ? jet_en->at(prompt_mujet)                          : -999;
      prompt_mujet_chEmEn_                     = (prompt_mujet != -1) ? jet_chargedEmEnergy->at(prompt_mujet)             : -999;
      prompt_mujet_NEmEnFraction_              = (prompt_mujet != -1) ? jet_neutralEmEnergyFraction->at(prompt_mujet)     : -999;
      prompt_mujet_chHadEn_                    = (prompt_mujet != -1) ? jet_chargedHadronEnergy->at(prompt_mujet)         : -999;
      prompt_mujet_NHadEnFraction_             = (prompt_mujet != -1) ? jet_neutralHadronEnergyFraction->at(prompt_mujet) : -999;
      prompt_mujet_chMuEn_                     = (prompt_mujet != -1) ? jet_chargedMuEnergy->at(prompt_mujet)             : -999;
      prompt_mujet_chMuEnFraction_             = (prompt_mujet != -1) ? jet_chargedMuEnergyFraction->at(prompt_mujet)     : -999;
      prompt_mujet_numberOfDaughters_          = (prompt_mujet != -1) ? jet_numberOfDaughters->at(prompt_mujet)           : -999;
      prompt_mujet_muonEnergy_                 = (prompt_mujet != -1) ? jet_muonEnergy->at(prompt_mujet)                  : -999;
      prompt_mujet_muonEnergyFraction_         = (prompt_mujet != -1) ? jet_muonEnergyFraction->at(prompt_mujet)          : -999;
      prompt_mujet_muonMultiplicity_           = (prompt_mujet != -1) ? jet_muonMultiplicity->at(prompt_mujet)            : -999;
      prompt_mujet_neutralEmEnergy_            = (prompt_mujet != -1) ? jet_neutralEmEnergy->at(prompt_mujet)             : -999;
      prompt_mujet_neutralHadronEnergy_        = (prompt_mujet != -1) ? jet_neutralHadronEnergy->at(prompt_mujet)         : -999;
      prompt_mujet_NHadMultiplicity_           = (prompt_mujet != -1) ? jet_neutralHadronMultiplicity->at(prompt_mujet)   : -999;
      prompt_mujet_NMultiplicity_              = (prompt_mujet != -1) ? jet_neutralMultiplicity->at(prompt_mujet)         : -999;
      prompt_mujet_chargedMultiplicity_        = (prompt_mujet != -1) ? jet_chargedMultiplicity->at(prompt_mujet)         : -999;
      prompt_mujet_chargedEmEnergyFraction     = (prompt_mujet != -1) ? jet_chargedEmEnergyFraction->at(prompt_mujet)     : -999;
      prompt_mujet_chargedHadronEnergyFraction = (prompt_mujet != -1) ? jet_chargedHadronEnergyFraction->at(prompt_mujet) : -999;
      prompt_mujet_notSmeard_pt                = (prompt_mujet != -1) ? jet_pt->at(prompt_mujet) : -999;
      prompt_mujet_raw_pt                      = (prompt_mujet != -1) ? jet_ptuncorrected->at(prompt_mujet) : -999;

      prompt_mujet_CsvV2          = (prompt_mujet != -1) ? jet_CsvV2->at(prompt_mujet)        : -999;
      prompt_mujet_DeepCsv_udsg   = (prompt_mujet != -1) ? jet_DeepCsv_udsg->at(prompt_mujet) : -999;
      prompt_mujet_DeepCsv_b      = (prompt_mujet != -1) ? jet_DeepCsv_b->at(prompt_mujet)    : -999;
      prompt_mujet_DeepCsv_c      = (prompt_mujet != -1) ? jet_DeepCsv_c->at(prompt_mujet)    : -999;
      prompt_mujet_DeepCsv_bb     = (prompt_mujet != -1) ? jet_DeepCsv_bb->at(prompt_mujet)   : -999;
      prompt_mujet_HadronFlavor   = (prompt_mujet != -1) ? jet_HadronFlavor->at(prompt_mujet) : -999;
      prompt_mujet_RSecMu        = (prompt_mujet != -1) ? minR : -999;

      //leading jet whithout first and second ele info 

      secondjet_charge_                     = (secondjet != -1) ? jet_charge->at(secondjet)                      : -999;
      secondjet_et_                         = (secondjet != -1) ? jet_et->at(secondjet)                          : -999;
      secondjet_pt_                         = (secondjet != -1) ? jetSmearedPt->at(secondjet)                    : -999;
      secondjet_eta_                        = (secondjet != -1) ? jet_eta->at(secondjet)                         : -999;
      secondjet_phi_                        = (secondjet != -1) ? jet_phi->at(secondjet)                         : -999;
      secondjet_theta_                      = (secondjet != -1) ? jet_theta->at(secondjet)                       : -999;
      secondjet_en_                         = (secondjet != -1) ? jet_en->at(secondjet)                          : -999;
      secondjet_chEmEn_                     = (secondjet != -1) ? jet_chargedEmEnergy->at(secondjet)             : -999;
      secondjet_NEmEnFraction_              = (secondjet != -1) ? jet_neutralEmEnergyFraction->at(secondjet)     : -999;
      secondjet_chHadEn_                    = (secondjet != -1) ? jet_chargedHadronEnergy->at(secondjet)         : -999;
      secondjet_NHadEnFraction_             = (secondjet != -1) ? jet_neutralHadronEnergyFraction->at(secondjet) : -999;
      secondjet_chMuEn_                     = (secondjet != -1) ? jet_chargedMuEnergy->at(secondjet)             : -999;
      secondjet_chMuEnFraction_             = (secondjet != -1) ? jet_chargedMuEnergyFraction->at(secondjet)     : -999;
      secondjet_numberOfDaughters_          = (secondjet != -1) ? jet_numberOfDaughters->at(secondjet)           : -999;
      secondjet_muonEnergy_                 = (secondjet != -1) ? jet_muonEnergy->at(secondjet)                  : -999;
      secondjet_muonEnergyFraction_         = (secondjet != -1) ? jet_muonEnergyFraction->at(secondjet)          : -999;
      secondjet_muonMultiplicity_           = (secondjet != -1) ? jet_muonMultiplicity->at(secondjet)            : -999;
      secondjet_neutralEmEnergy_            = (secondjet != -1) ? jet_neutralEmEnergy->at(secondjet)             : -999;
      secondjet_neutralHadronEnergy_        = (secondjet != -1) ? jet_neutralHadronEnergy->at(secondjet)         : -999;
      secondjet_NHadMultiplicity_           = (secondjet != -1) ? jet_neutralHadronMultiplicity->at(secondjet)   : -999;
      secondjet_NMultiplicity_              = (secondjet != -1) ? jet_neutralMultiplicity->at(secondjet)         : -999;
      secondjet_chargedMultiplicity_        = (secondjet != -1) ? jet_chargedMultiplicity->at(secondjet)         : -999;
      secondjet_chargedEmEnergyFraction     = (secondjet != -1) ? jet_chargedEmEnergyFraction->at(secondjet)     : -999;
      secondjet_chargedHadronEnergyFraction = (secondjet != -1) ? jet_chargedHadronEnergyFraction->at(secondjet) : -999;

      secondjet_notSmeard_pt                = (secondjet != -1) ? jet_pt->at(secondjet) : -999;
      secondjet_raw_pt                      = (secondjet != -1) ? jet_ptuncorrected->at(secondjet) : -999;

      //cout<<"jet not smearing pt = "<<secondjet_notSmeard_pt<<" jet raw pt = "<<secondjet_raw_pt<<" eta = "<<secondjet_eta_ <<" phi " << secondjet_phi_ <<" Rho = " <<ele_second_rhoIso<<endl; 

      secondjet_CsvV2          = (secondjet != -1) ? jet_CsvV2->at(secondjet)        : -999;
      secondjet_DeepCsv_udsg   = (secondjet != -1) ? jet_DeepCsv_udsg->at(secondjet) : -999;
      secondjet_DeepCsv_b      = (secondjet != -1) ? jet_DeepCsv_b->at(secondjet)    : -999;
      secondjet_DeepCsv_c      = (secondjet != -1) ? jet_DeepCsv_c->at(secondjet)    : -999;
      secondjet_DeepCsv_bb     = (secondjet != -1) ? jet_DeepCsv_bb->at(secondjet)   : -999;
      secondjet_HadronFlavor   = (secondjet != -1) ? jet_HadronFlavor->at(secondjet) : -999;

      secondjets_size  = (secondjet != -1) ? secondjet_count : -999;
      secondjet_bjet_L = bjet111;
      secondjet_bjet_M = bjet222;
      secondjet_bjet_T = bjet333;
    //missing energy Info
      pfMet_et_    =  pfMet_et;
      pfMet_pt_    =  pfMet_pt;
      pfMet_phi_   =  pfMet_phi;
      pfMet_en_    =  pfMet_en;
      pfMet_sumEt_ =  pfMet_sumEt;
      caloMet_pt_  =  caloMet_pt;
      caloMet_phi_ =  caloMet_phi;
      
      TLorentzVector MET;
      MET.SetPxPyPzE(pfMet_px,pfMet_py,pfMet_pz,pfMet_en_);

      //tmass_METL = sqrt(pow(mu_promptEt + pfMet_et, 2) - (pow(mu_promptPt + pfMet_pt, 2)));

      tmass_METL = sqrt(pow(mu_promptPt + pfMet_pt, 2) - pow(px1 + pfMet_px, 2) - pow(py1 + pfMet_py, 2));
      float MetLMass = (Mu1 + MET).M();

      Mass_METL = MetLMass;
      newtree->Fill();
    }
    
  }

  h_nTrueInteractions50->Write();
  h_nTrueInteractions100->Write();

  //newfile->Write(); 
  //newtree->Print();
  newtree->AutoSave();

  delete newfile;
  return 0;
}

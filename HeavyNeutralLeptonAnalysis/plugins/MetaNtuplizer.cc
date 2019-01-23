// -*- C++ -*-
//
// Package:    URNtuples/MetaNtuplizer.cc
// Class:      MetaNtuplizer.cc
// 
/**\class MetaNtuplizer.cc MetaNtuplizer.cc.cc URNtuples/MetaNtuplizer.cc/plugins/MetaNtuplizer.cc.cc

 Description: Computes and stores in the rootfile the Meta Information needed

*/
//
// Original Author:  Mauro Verzetti
//         Created:  Thu, 20 Nov 2014 11:34:09 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DQMServices/Core/interface/DQMStore.h>
#include <DQMServices/Core/interface/MonitorElement.h>
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "TTree.h"
#include "TObjString.h"
#include "TH1F.h"

#include <map>
#include <string>
#include <iostream>
#include <sstream> 

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

using namespace std;
//
// class declaration
//

class MetaNtuplizer : public edm::EDAnalyzer {
public:
  explicit MetaNtuplizer(const edm::ParameterSet&);
  ~MetaNtuplizer() {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::InputTag weights_src_;
  edm::EDGetTokenT< LHEEventProduct > weights_srcToken_;
  edm::EDGetTokenT< LHERunInfoProduct > header_token_;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pu_token_;

  TTree *meta_tree_;
  std::map<std::string, std::string> to_json_;
  bool string_dumped_, isMC_, hasLhe_, useWeighted_, triedWeighted_;
  TH1F *pu_distro_;
  unsigned int lumi_;
  unsigned int run_;
  unsigned long long processed_ = 0;
  long long processedWeighted_ = 0; 
  vector<double> sumw_;
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
MetaNtuplizer::MetaNtuplizer(const edm::ParameterSet& iConfig):
  weights_src_(iConfig.getParameter<edm::InputTag>("weightsSrc") ),
  weights_srcToken_(consumes<LHEEventProduct>(weights_src_)),	
  header_token_(consumes<LHERunInfoProduct,edm::InRun>(weights_src_)),
  pu_token_(consumes< std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puSrc"))),
  string_dumped_(false),
  isMC_(iConfig.getParameter<bool>("isMC")),
	hasLhe_(iConfig.getParameter<bool>("hasLHE")),
	sumw_()
{
  useWeighted_ = true;
  triedWeighted_ = false;
  //dump direct information
  to_json_.insert(std::make_pair<std::string, std::string>("tuple_commit", iConfig.getParameter<std::string>("commit"))); 
  to_json_.insert(std::make_pair<std::string, std::string>("tuple_user", iConfig.getParameter<std::string>("user"))); 
  to_json_.insert(std::make_pair<std::string, std::string>("tuple_cmsswVersion", iConfig.getParameter<std::string>("cmsswVersion"))); 
  to_json_.insert(std::make_pair<std::string, std::string>("tuple_date", iConfig.getParameter<std::string>("date"))); 
  to_json_.insert(std::make_pair<std::string, std::string>("tuple_globalTag", iConfig.getParameter<std::string>("globalTag")));
  to_json_.insert(std::make_pair<std::string, std::string>("tuple_args", iConfig.getParameter<std::string>("args")));

  edm::Service<TFileService> fs;
  meta_tree_ = fs->make<TTree>( "meta"  , "File Meta Information");
  meta_tree_->Branch("run", &run_);
  meta_tree_->Branch("lumi", &lumi_);
  meta_tree_->Branch("processed", &processed_);
  meta_tree_->Branch("processedWeighted", &processedWeighted_);
  meta_tree_->Branch("sum_weigts", &sumw_);

  pu_distro_   = fs->make<TH1F>("PUDistribution", "PUDistribution", 100, 0, 100);
}

//
// member functions
//

// ------------ method called once each job just after ending the event loop  ------------
void MetaNtuplizer::beginJob()
{

}

void MetaNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup&)
{
	processed_++;
	double weight = 1.;

	int npu = -1;
	if(isMC_) {
	  edm::Handle< std::vector<PileupSummaryInfo> > pu_info;
	  iEvent.getByToken(pu_token_, pu_info);

	  for(const auto& PVI : *pu_info) {
	      int BX = PVI.getBunchCrossing();
	      if(BX == 0) {
		npu = PVI.getTrueNumInteractions();
		pu_distro_->Fill(npu);
		break;
	      }
	  }
	}

	if(hasLhe_) //lheinfo.isValid() && lheinfo->weights().size() > 0)
	{
		edm::Handle<LHEEventProduct> lheinfo;
		iEvent.getByToken(weights_srcToken_, lheinfo);
		if(!lheinfo.isValid()) {
			throw cms::Exception("RuntimeError") << "The Handle of LHEInfo I was trying to access is not valid!" << std::endl;
		}
		if(lheinfo->weights().size() == 0) 
			throw cms::Exception("RuntimeError") << "The LHEInfo I got works but has not weights!" << std::endl;

		//std::cout << weight << std::endl;
		size_t nws = lheinfo->weights().size();
		if(!sumw_.size()) {
			sumw_.reserve(nws);
			for(size_t i=0; i<nws; ++i) {
				sumw_.push_back(0.);
			}
		} else {
			if(nws != sumw_.size())
				throw cms::Exception("RuntimeError") << "I set up for " << sumw_.size() << 
					" LHE weights, but this event has " << nws << "!" << std::endl;
		}

		weight = lheinfo->weights()[0].wgt;
		for(size_t i=0; i<nws; ++i) {
			sumw_[i] += lheinfo->weights()[i].wgt;
		}
	}
	processedWeighted_ += (weight < 0. ? -1. : 1.);

}

void MetaNtuplizer::endJob() 
{
  edm::Service<TFileService> fs;

  std::stringstream stream;
  stream << "{" << std::endl;
  for(auto entry = to_json_.begin(); entry != to_json_.end(); ++entry)
    {
      stream << "   \"" << entry->first << "\" : \"" << entry->second << "\"," << std::endl;
    }
  stream << "}" << std::endl;
  
  fs->make<TObjString>(stream.str().c_str());
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MetaNtuplizer::endLuminosityBlock(edm::LuminosityBlock const& block, edm::EventSetup const&)
{
  lumi_ = block.luminosityBlock();
  run_ = block.run();

  meta_tree_->Fill();
  processed_ = 0;
  processedWeighted_ = 0;
	for(size_t i=0; i<sumw_.size(); ++i) sumw_[i] = 0;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MetaNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_FWK_MODULE(MetaNtuplizer);

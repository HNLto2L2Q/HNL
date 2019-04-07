#include "HNL/DisplacedAdaptiveVertexFinder/plugins/DisplacedInclusiveVertexFinder.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"


typedef TemplatedInclusiveVertexFinder<reco::TrackCollection,reco::Vertex> DisplacedInclusiveVertexFinder;
typedef TemplatedInclusiveVertexFinder<edm::View<reco::Candidate>,reco::VertexCompositePtrCandidate > InclusiveCandidateVertexFinder;


DEFINE_FWK_MODULE(DisplacedInclusiveVertexFinder);
DEFINE_FWK_MODULE(InclusiveCandidateVertexFinder);

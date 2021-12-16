// -*- C++ -*-
//
// Package:    DarkSUSY/fastJetNtupler
// Class:      fastJetNtupler
//
/**\class fastJetNtupler fastJetNtupler.cc DarkSUSY/fastJetNtupler/plugins/fastJetNtupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Balashangar Kailasapathy
//         Created:  Thu, 16 Dec 2021 03:57:05 GMT
//
//


// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//======Additional header files==================
#include "FWCore/Framework/interface/EventSetup.h" 
#include "FWCore/Framework/interface/ESHandle.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" //to work with reco::GenParticle
#include "DataFormats/PatCandidates/interface/Jet.h" //to work with pat::Jets
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //to work with particles inside jets
#include <vector> 
#include <string> 
#include <map>
#include <iostream>
#include <fstream>
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/ServiceRegistry/interface/Service.h" // to use TFileService
#include "CommonTools/UtilAlgos/interface/TFileService.h" // to use TFileService

//Muons
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"      


//Missing Energy
//==============
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/CorrMETData.h"


#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"


//Header file for Fastjet Analysis
//====================================

#include "fastjet/ClusterSequence.hh"
#include "fastjet/config.h"
#include "fastjet/SISConePlugin.hh"

//#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DarkSUSY/fastJetNtupler/interface/Myheaderfile1.h"       // My header file
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"

void print_jets (const std::vector<fastjet::PseudoJet> &);
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class fastJetNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit fastJetNtupler(const edm::ParameterSet&);
      ~fastJetNtupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  	edm::EDGetTokenT<std::vector<pat::Muon>			> patmuonToken;
	edm::EDGetTokenT<std::vector<pat::MET>			> patMetToken;
		
	edm::Service<TFileService> fs;
	
	TFile *f = new TFile("myfile.root","RECREATE");
	TTree *T = new TTree("T","A ROOT tree ");
	
	float muonpt;
	int incJetSize;
	
	float jetrap;
	float jetphi;
	float jetpt;
	float invmass;
	
	float metpt;
	float metpx;
	float metpy;
	float metphi; 
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
fastJetNtupler::fastJetNtupler(const edm::ParameterSet& iConfig)
 :
  patmuonToken(consumes<std::vector<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>	("muonTag"))),
				patMetToken(consumes<std::vector<pat::MET> >(iConfig.getUntrackedParameter<edm::InputTag>	("metTag")))

{
   //now do what ever initialization is needed

}


fastJetNtupler::~fastJetNtupler()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
	f->Write();
	f->Close();

}


//
// member functions
//

// ------------ method called for each event  ------------
void
fastJetNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   	edm::Handle<std::vector<pat::Muon>> patMuon;
	iEvent.getByToken(patmuonToken, patMuon);
	
	edm::Handle<std::vector<pat::MET>> patMet;
	iEvent.getByToken(patMetToken, patMet); 
	
	std::vector<fastjet::PseudoJet> input_particles;
	std::vector<fastjet::PseudoJet> selected_input_particles;
	
//Fastjet------------------------------------	
	
		for (std::vector<pat::Muon>::const_iterator itMuon=patMuon->begin(); itMuon!=patMuon->end(); ++itMuon) 
		{
			muonpt = itMuon->pt();
			T->Branch("muonpt", &muonpt,"muonpt/F");
			T->Fill();
			
					
			input_particles.push_back(fastjet::PseudoJet(itMuon->px(),itMuon->py(),itMuon->pz(),itMuon->energy()));
		}
	 
		std::cout <<  " Number of particles before applying cuts (ie, before using selector) : " <<input_particles.size() << std::endl;
	
		fastjet::Selector particle_selector = fastjet::SelectorAbsRapRange(1.0,2.5);
		selected_input_particles = particle_selector(input_particles);   
		std::cout <<  " Number of particles after applying cuts : " <<selected_input_particles.size() << std::endl;
		
		// ------Jet definition
		double R = 0.7;
		fastjet::JetDefinition jet_def(fastjet::kt_algorithm,R);
		std::cout<<"Jet definition used here: "<<jet_def.description()<<std::endl;
	
	
		fastjet::ClusterSequence clust_seq(input_particles, jet_def);  //Jet clustering for particles with above jet definition

		//Jet Selector
		fastjet::Selector jet_selector = fastjet::SelectorNHardest(5) * fastjet::SelectorAbsRapMax(2.0);
		
	
		double ptmin = 5.0;
		std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin)); //get the resulting jets ordered in pt
		
		
		std::cout<< "Number  of Inclusive jets = "<<inclusive_jets.size()<<std::endl;
		incJetSize = inclusive_jets.size();

		T->Branch("incJetSize", &incJetSize,"incJetSize/i");
		T->Fill();
		
			
		int pdg_id = 13;   	//  - pdg_id        the PDG id of the particle
   		int vertex_no = 1;	//  - vertex_number the id of the vertex it originates from
		
		for (unsigned int i = 0; i < inclusive_jets.size(); i++)
		{
			Myheaderfile1* infojet = new Myheaderfile1(pdg_id,vertex_no);	
			inclusive_jets[i].set_user_info(infojet);
			const int & pdgid = inclusive_jets[i].user_info<Myheaderfile1>().pdg_id();
		
			//std::cout<<"This is the pgdID: "<<pdgid<<std::endl;
		}
		
		double rapmin, rapmax;

		std::cout << "Selected particles: " << particle_selector.description() << std::endl;
		particle_selector.get_rapidity_extent(rapmin, rapmax);
		std::cout << "  with a total rapidity range of [" << rapmin << ", " << rapmax << "]" << std::endl;

		std::cout << "Selected jets: " << jet_selector.description() << std::endl;
		jet_selector.get_rapidity_extent(rapmin, rapmax);
		std::cout << "  with a total rapidity range of [" << rapmin << ", " << rapmax << "]" << std::endl;
		
		
		printf("\n ");
		printf("Details for each jet\n");
		printf("%5s %15s %15s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt" ,"pt_perp", "invMass" , "n Constituents");   // label the columns
		
		for (unsigned int i = 0; i < inclusive_jets.size(); i++)
		{
			std::vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();
			
			
			printf("%5u %15.8f %15.8f %15.8f %15.8f %15.8f %15.8u\n",i, inclusive_jets[i].rap(), inclusive_jets[i].phi(), inclusive_jets[i].pt(),inclusive_jets[i].perp(), inclusive_jets[i].m(), (unsigned int) constituents.size());
			
			jetrap  = inclusive_jets[i].rap();
			jetphi = inclusive_jets[i].phi();
			jetpt = inclusive_jets[i].pt();
			invmass = inclusive_jets[i].m();
			
			T->Branch("jetrap", &jetrap,"jetrap/F");
			T->Branch("jetphi", &jetphi,"jetphi/F");
			T->Branch("jetpt", &jetpt,"jetpt/F");
			T->Branch("invmass", &invmass,"invmass/F");
			
			T->Fill();
					
		}

		//Exclusive Jets
		if (incJetSize >2)      
		{
			std::vector<fastjet::PseudoJet> exclusive_jets = clust_seq.exclusive_jets(1);
			std::cout<< "Number of Exclusive jets = "<<exclusive_jets.size()<<std::endl;
			int ExcJetSize = exclusive_jets.size();
			
		}
		
		
	const pat::MET &met = patMet->front();
	
	metpt=met.pt();
	metpx=met.px();
	metpy=met.py();
	metphi=met.phi();
	
	T->Branch("metpt", &metpt,"metpt/F");
	T->Branch("metpx", &metpx,"metpx/F");
	T->Branch("metpy", &metpx,"metpy/F");
	T->Branch("metphi", &metpx,"metphi/F");
	
	T->Fill();
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
fastJetNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
fastJetNtupler::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
fastJetNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(fastJetNtupler);

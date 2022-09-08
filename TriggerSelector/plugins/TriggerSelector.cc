// -*- C++ -*-
//
// Package:    DarkSUSY/TriggerSelector
// Class:      TriggerSelector
// 
/**\class TriggerSelector TriggerSelector.cc DarkSUSY/TriggerSelector/plugins/TriggerSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Balashangar Kailasapathy
//         Created:  Thu, 28 Jul 2022 11:57:32 GMT

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include <string>

// user include files

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" 

#include "TH1.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"


//
// class declaration
//


class TriggerSelector : public edm::stream::EDFilter<> {
   public:
      explicit TriggerSelector(const edm::ParameterSet&);
      ~TriggerSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;


	edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
			  
	  // TFileService
	edm::Service<TFileService> fs;

	std::vector<std::string> triggerList;
	std::map<TString, size_t> triggerIdxList; 
	
	unsigned int *passtrig_;
	std::map<std::string, int> counts;  
	int num;
	
	//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    // ----------member data ---------------------------
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
TriggerSelector::TriggerSelector(const edm::ParameterSet& iConfig)
{
 //now do what ever initialization is needed

  edm::InputTag triggerTag("TriggerResults", "", "HLT");
 // edm::InputTag trigObjTag("selectedPatTrigger");

 trigResultsToken = consumes<edm::TriggerResults>(triggerTag);
 // trigObjCollToken = consumes<pat::TriggerObjectStandAloneCollection>(trigObjTag);

  triggerList = iConfig.getParameter<std::vector<std::string> >("TriggerList");
  passtrig_ = new unsigned int();
  
  
}


TriggerSelector::~TriggerSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)


	std::vector<std::pair<std::string, int > > vec;
	for (auto i = counts.begin();i!=counts.end();i++){
	vec.push_back(make_pair(i->first, i->second));
		}
  //sort(vec.begin(),vec.end(),sortByCount);

  std::cout <<" Results "<<std::endl;
	for (auto i = vec.begin(); i!= vec.end();i++){
    std::cout <<" ++ "<<std::left<<std::setw(80)<<i->first<<std::right<<std::setw(20)<<i->second<<std::endl;
	} 
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool TriggerSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  // Retrieve TriggerResults 
  
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(trigResultsToken, triggerBits);
  
  
  if(!triggerBits.isValid())
  {
    throw cms::Exception("TriggerResults collection not valid!"); 
  } 
  
  if (triggerBits.isValid())
  {
  std::cout << "TriggerResults found, number of HLT paths: " << triggerBits->size()<<std::endl;
 
  const edm::TriggerNames& trns = iEvent.triggerNames(*triggerBits);
    std::vector<std::string>  names = trns.triggerNames();
    edm::TriggerResultsByName resultsByNameHLT = iEvent.triggerResultsByName(*triggerBits);
	
	
	
	
  for (auto inm = names.begin(); inm!= names.end(); inm++)
  {
	  
     //std::cout <<*inm<<std::endl;
	 //std::cout<<"Is it fired?  --  " <<resultsByNameHLT.accept(*inm)<<std::endl;
      
	 if (resultsByNameHLT.accept(*inm)) //&&	
		//inm->find("HLT") == 0 &&
		//inm->find("ZeroBias") ==  std::string::npos &&
		//inm->find("Physics") ==  std::string::npos &&
		//inm->find("Random") ==  std::string::npos)
	  {
		if (counts.find(*inm) == counts.end())
		{
		counts[*inm]=0;
		}
		num = counts[*inm]++;
		std::cout<<"this is the count   :"<<num<<std::endl;
		}
  }
}
  
  
  
  /*
  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerBits);

	for(unsigned int itrg = 0; itrg < triggerBits->size(); ++itrg) {
		TString trigname = triggerNames.triggerName(itrg);
	bool accept = triggerBits->accept(itrg);
    if(accept) passtrig_[itrg] = 1;
    else       passtrig_[itrg] = 0;
	
	std::cout<<itrg<<"  "<<trigname<<std::setw(30) <<*passtrig_<<std::endl;
  }


  // Find trigger indexes (only once)
  if(triggerIdxList.size()==0)
  { 
     for(auto&& t : triggerList)
     { 
         for(size_t i=0, n=triggerBits->size(); i<n; ++i)
         {
		 	if(TString(triggerNames.triggerName(i)).Contains((t+"_v").c_str()))
			 { 
				triggerIdxList[(t+"_v").c_str()] = i;
				std::cout<<(i)<<std::endl;
	         break; 
             } 
         }
     }
  }

 */  


   using namespace edm;
   
   
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void TriggerSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void TriggerSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TriggerSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TriggerSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TriggerSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TriggerSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TriggerSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  //Please change this to state exactly what you do use, even if it is no parameters
  
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerSelector);

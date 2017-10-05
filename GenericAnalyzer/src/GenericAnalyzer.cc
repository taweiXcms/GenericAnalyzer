// Original Author:  Ta-Wei Wang,,,
// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "TH1F.h"
#include "TLorentzVector.h"

#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677
#define KSHORT_MASS 0.497614
#define KSTAR_MASS  0.89594
#define PHI_MASS    1.019455
#define JPSI_MASS   3.096916
#define PSI2S_MASS  3.686109

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
//#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
//
// class declaration
//

class GenericAnalyzer : public edm::EDAnalyzer {
	public:
	explicit GenericAnalyzer(const edm::ParameterSet&);
	~GenericAnalyzer();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	void ShowMyOffspring(const reco::Candidate* p, int layer);

	private:
	virtual void beginJob() ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void endRun(edm::Run const&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
	
	// ----------member data ---------------------------
//    edm::InputTag genLabel_;
    edm::EDGetTokenT< reco::GenParticleCollection > genLabel_;
    edm::InputTag muonLabel_;
    edm::InputTag trackLabel_;
    edm::InputTag jetLabel_;

    bool doHepMC_; 
    bool doGenParticle_; 
    bool doGenJet_; 
    bool doRecoTrk_; 
    bool doRecoMuon_; 

	//My output histogram
	TH1F *h_mass1,*h_mass2;
	TH1F *Bpt, *Beta;
	TH1F *Bupt, *Bueta;
	TH1F *Bdpt, *Bdeta;
	TH1F *Bspt, *Bseta;
	TH1F *genmu_eta;
	TH1F *gensig_eta;
	TH1F *genemb_eta;
};

GenericAnalyzer::GenericAnalyzer(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
//    genLabel_           = iConfig.getParameter<edm::InputTag>("GenLabel");
    genLabel_           = consumes< reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>("GenLabel"));
    trackLabel_         = iConfig.getParameter<edm::InputTag>("TrackLabel");
    muonLabel_          = iConfig.getParameter<edm::InputTag>("MuonLabel");
    jetLabel_          = iConfig.getParameter<edm::InputTag>("JetLabel");

    doHepMC_ = iConfig.getUntrackedParameter<bool>("doHepMC",true);
    doGenParticle_ = iConfig.getUntrackedParameter<bool>("doGenParticle",true);
    doGenJet_ = iConfig.getUntrackedParameter<bool>("doGenJet",true);
    doRecoTrk_ = iConfig.getUntrackedParameter<bool>("doRecoTrk",true);
    doRecoMuon_ = iConfig.getUntrackedParameter<bool>("doRecoMuon",true);
	edm::Service<TFileService> fs;

	//My output histogram
	// 1 => track pairs, 2 => muon pairs
	h_mass1 = fs->make<TH1F>("h_mass1" , "h_mass1" , 320 , 0. , 10. );
	h_mass2 = fs->make<TH1F>("h_mass2" , "h_mass2" , 3200 , 0. , 160. );
	Bpt = fs->make<TH1F>("Bpt" , "Bpt" , 200 , 0. , 20. );
	Beta = fs->make<TH1F>("Beta" , "Beta" , 200 , -6. , 6. );
	Bupt = fs->make<TH1F>("Bupt" , "Bpt" , 200 , 0. , 20. );
	Bueta = fs->make<TH1F>("Bueta" , "Beta" , 200 , -6. , 6. );
	Bdpt = fs->make<TH1F>("Bdpt" , "Bpt" , 200 , 0. , 20. );
	Bdeta = fs->make<TH1F>("Bdeta" , "Beta" , 200 , -6. , 6. );
	Bspt = fs->make<TH1F>("Bspt" , "Bpt" , 200 , 0. , 20. );
	Bseta = fs->make<TH1F>("Bseta" , "Beta" , 200 , -6. , 6. );
	genmu_eta = fs->make<TH1F>("genmu_eta" , "genmu_eta" , 200 , -6. , 6. );
	gensig_eta = fs->make<TH1F>("gensig_eta" , "gensig_eta" , 200 , -6. , 6. );
	genemb_eta = fs->make<TH1F>("genemb_eta" , "genemb_eta" , 200 , -10. , 10. );
}

GenericAnalyzer::~GenericAnalyzer()
{
}

// ------------ method called for each event  ------------
void
GenericAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace std;
	using namespace edm;
	
	cout<<"**********"<<endl;

	if(doHepMC_){
		Handle<HepMCProduct> mc;
	    iEvent.getByLabel("generator",mc);
		const HepMC::GenEvent* evt; 
	    evt = mc->GetEvent();
		for(HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it){
	    	if((*it)->status() == 1 || (*it)->status() == 2){
	        //int pdg_id = (*it)->pdg_id();
	        //float eta = (*it)->momentum().eta();
			}
		}
	}

	if(doGenParticle_){
//	    edm::Handle<reco::GenParticleCollection> parts;
//	    iEvent.getByLabel(genLabel_,parts);
        edm::Handle<reco::GenParticleCollection> parts;
        iEvent.getByToken(genLabel_, parts);
		bool isSig = 0;
	    bool GetGammaSignal = false;
	    bool GetBpSignal = false;
	    bool GetB0Signal = false;
	    bool GetBsSignal = false;
	    bool GetD0Signal = false;
	    bool GetDsSignal = false;
		////My analysis
	    for(unsigned int i = 0; i < parts->size(); ++i){
	    	const reco::GenParticle& p = (*parts)[i];
	        //if(p.status() != 1 && p.status() != 2) continue;
	        //float pt = p.pt();
	        int pdg = p.pdgId();

			//Printing out gen info
			//if(abs(pdg) == 531 || abs(pdg) == 521 || abs(pdg) == 511 || (abs(pdg) == 11 && p.mother()->pdgId() == 310) || (abs(pdg) == 310 && abs(p.mother()->pdgId()) == 311) ){
			//if(abs(pdg) == 531 || abs(pdg) == 521 || abs(pdg) == 511 ){
			//if(abs(pdg) == 421 || abs(pdg) == 22){
			if(abs(pdg) == 431){
			//if(abs(int(pdg/100) % 100) == 5){
			//if(1){
				if(p.numberOfMothers() != 0) cout<<"======my mother Id/pt: "<<p.mother()->pdgId()<<" / "<<p.mother()->pt()<<endl;
				int ndau = p.numberOfDaughters();
				cout<<"I am(pdg/status/pt/#daughters): "<<pdg<<" / "<<p.status()<<" / "<<p.pt()<<" / "<<ndau<<endl;
				if(ndau != 0){
					for(int i = 0; i < ndau; i++){
						const reco::Candidate* MyDau = (p.daughter(i));
	                    int ndaudau = MyDau->numberOfDaughters();
	                    cout<< "   daughter: " << MyDau->pdgId() << " / " << MyDau->status() <<" / "<<MyDau->pt()<<" / "<<ndaudau<<endl;
						ShowMyOffspring(MyDau, 2);
					}
				}
			}
			//Saving Gamma info
//			if(abs(pdg) == 22){	
			if(abs(pdg) == 22 && p.pt()>2.){	
					GetGammaSignal = true;
			}

			//Saving D0 info
			if(abs(pdg) == 421){	
				if(p.numberOfDaughters()>1)
				if(abs(p.daughter(0)->pdgId()) == 321) 
				if(abs(p.daughter(1)->pdgId()) == 211) 
				{
					Bpt->Fill(p.pt());
					Beta->Fill(p.eta());
					Bupt->Fill(p.pt());
					Bueta->Fill(p.eta());
					GetD0Signal = true;
				}
			}

			//Saving Ds info
			if(abs(pdg) == 431){	
				if(p.numberOfDaughters()>1)
				if(abs(p.daughter(0)->pdgId()) == 333) 
				if(abs(p.daughter(1)->pdgId()) == 211) 
				if(p.daughter(0)->numberOfDaughters()>1)
				if(abs(p.daughter(0)->daughter(0)->pdgId())==321)
				if(abs(p.daughter(0)->daughter(1)->pdgId())==321)
				{
					Bpt->Fill(p.pt());
					Beta->Fill(p.eta());
					Bupt->Fill(p.pt());
					Bueta->Fill(p.eta());
					GetDsSignal = true;
				}
			}

			//Saving B+ info
			if(abs(pdg) == 521){	
				if(p.numberOfDaughters()>1)
				if(abs(p.daughter(0)->pdgId()) == 443) 
				if(abs(p.daughter(1)->pdgId()) == 321) 
				if(p.daughter(0)->numberOfDaughters()>1)
				if(abs(p.daughter(0)->daughter(0)->pdgId())==13)
				if(abs(p.daughter(0)->daughter(1)->pdgId())==13)
				{
					Bpt->Fill(p.pt());
					Beta->Fill(p.eta());
					Bupt->Fill(p.pt());
					Bueta->Fill(p.eta());
					GetBpSignal = true;
				}
			}

			//Saving B0 info
            if(abs(pdg) == 511){
				if(p.numberOfDaughters()>1)
				if(abs(p.daughter(0)->pdgId()) == 443) 
				if(p.daughter(0)->numberOfDaughters()>1)
				if(abs(p.daughter(0)->daughter(0)->pdgId())==13)
				if(abs(p.daughter(0)->daughter(1)->pdgId())==13)
				if(abs(p.daughter(1)->pdgId()) == 313) 
				if(p.daughter(1)->numberOfDaughters()>1)
				if(abs(p.daughter(1)->daughter(0)->pdgId())==211 || abs(p.daughter(1)->daughter(0)->pdgId())==321)
				if(abs(p.daughter(1)->daughter(1)->pdgId())==211 || abs(p.daughter(1)->daughter(1)->pdgId())==321)
				{
					Bpt->Fill(p.pt());
					Beta->Fill(p.eta());
					Bdpt->Fill(p.pt());
					Bdeta->Fill(p.eta());
					GetB0Signal = true;
				}
			}

			//Saving Bs info
            if(abs(pdg) == 531){
				if(p.numberOfDaughters()>1)
				if(abs(p.daughter(0)->pdgId()) == 443) 
				if(p.daughter(0)->numberOfDaughters()>1)
				if(abs(p.daughter(0)->daughter(0)->pdgId())==13)
				if(abs(p.daughter(0)->daughter(1)->pdgId())==13)
				if(abs(p.daughter(1)->pdgId()) == 333) 
				if(p.daughter(1)->numberOfDaughters()>1)
				if(abs(p.daughter(1)->daughter(0)->pdgId())==321)
				if(abs(p.daughter(1)->daughter(1)->pdgId())==321)
				{
					Bpt->Fill(p.pt());
					Beta->Fill(p.eta());
					Bspt->Fill(p.pt());
					Bseta->Fill(p.eta());
					GetBsSignal = true;
				}
			}
			//Saving muon eta
			if(abs(pdg) == 13){	
				if(abs(p.mother()->pdgId()) == 443) if(abs(p.mother()->mother()->pdgId()) == 521)//B+
				//if(abs(p.mother()->pdgId()) == 443) if(abs(p.mother()->mother()->pdgId()) == 511)//B0
				genmu_eta->Fill(p.eta());
			}

			////filling only signal particles
			//B+
			if(abs(pdg) == 521 || abs(pdg) == 13 || abs(pdg) == 443){
				if(abs(pdg) == 521 && abs(p.daughter(0)->pdgId()) == 443) isSig = 1;		
				if(abs(pdg) == 443 && abs(p.mother()->pdgId()) == 521) isSig = 1;
				if(abs(pdg) == 13 && abs(p.mother()->pdgId()) == 443) if(abs(p.mother()->mother()->pdgId()) == 521) isSig = 1;
			}
			if(p.pt() > 0){
				if(isSig == 1)	{
					gensig_eta->Fill(p.eta());
					//cout << pdg << endl;
				}
				else if(p.status() == 1) genemb_eta->Fill(p.eta());
			}
			isSig = 0;
		}
		if(GetGammaSignal) cout<<"GetGammaSignal"<<endl;
		else cout<<"NoGammaSignal"<<endl;
		if(GetD0Signal) cout<<"GetD0Signal"<<endl;
		else cout<<"NoD0Signal"<<endl;
		if(GetDsSignal) cout<<"GetDsSignal"<<endl;
		else cout<<"NoDsSignal"<<endl;
		if(GetBpSignal) cout<<"GetBpSignal"<<endl;
		else cout<<"NoBpSignal"<<endl;
		if(GetB0Signal) cout<<"GetB0Signal"<<endl;
		else cout<<"NoB0Signal"<<endl;
		if(GetBsSignal) cout<<"GetBsSignal"<<endl;
		else cout<<"NoBsSignal"<<endl;
	}

	if(doGenJet_){
	    edm::Handle<vector<reco::GenJet> > jets;
	    iEvent.getByLabel(jetLabel_,jets);
	    for(unsigned int j = 0 ; j < jets->size(); ++j){
	       //const reco::GenJet& jet = (*jets)[j];
	       //float jet_p = jet.p();
		}
	}

	if(doRecoTrk_){
		Handle<reco::TrackCollection> tracks;
		iEvent.getByLabel(trackLabel_, tracks);
		for( reco::TrackCollection::const_iterator it1 = tracks->begin();
			it1 != tracks->end(); it1++ ) {
			for( reco::TrackCollection::const_iterator it2 = tracks->begin();
				it2 != tracks->end(); it2++ ) {
				if (it1->charge()<0 || it2->charge()>0) continue;
				TLorentzVector v1(it1->px(),it1->py(),it1->pz(),sqrt(it1->p()*it1->p()+PION_MASS*PION_MASS));
				TLorentzVector v2(it2->px(),it2->py(),it2->pz(),sqrt(it2->p()*it2->p()+PION_MASS*PION_MASS));
				TLorentzVector v = v1 + v2;
				h_mass1->Fill(v.Mag());
			}
		}
	}
	
	if(doRecoMuon_){
		Handle<reco::MuonCollection> muons;
		iEvent.getByLabel(muonLabel_, muons);
		for( reco::MuonCollection::const_iterator it1 = muons->begin();
			it1 != muons->end(); it1++ ) {
			//std::cout<<"px:"<<it1->px()<<std::endl;
			for( reco::MuonCollection::const_iterator it2 = muons->begin();
				it2 != muons->end(); it2++ ) {
				if (it1->charge()<0 || it2->charge()>0) continue;
				TLorentzVector v1(it1->px(),it1->py(),it1->pz(),sqrt(it1->p()*it1->p()+MUON_MASS*MUON_MASS));
				TLorentzVector v2(it2->px(),it2->py(),it2->pz(),sqrt(it2->p()*it2->p()+MUON_MASS*MUON_MASS));
				TLorentzVector v = v1 + v2;
				h_mass2->Fill(v.Mag());
			}
		}
	}

	//Do SIM track. Developing
	/*
	if(0){
	    Handle<vector<reco::Track> > etracks;
	    iEvent.getByLabel(trackLabel_, etracks);
	    Handle<TrackingParticleCollection>  TPCollectionHfake;		
	    Handle<edm::View<reco::Track> >  trackCollection;
	    iEvent.getByLabel(trackLabel_, trackCollection);
	    ESHandle<TrackAssociatorBase> theAssociator;
	    const TrackAssociatorByHits *theAssociatorByHits;
	    reco::RecoToSimCollection recSimColl;
	    iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
	    theAssociatorByHits = (const TrackAssociatorByHits*) theAssociator.product();
	    recSimColl= theAssociatorByHits->associateRecoToSim(trackCollection,TPCollectionHfake,&iEvent);
	    for(unsigned it=0; it<etracks->size(); ++it){
		    reco::TrackRef trackRef=reco::TrackRef(etracks,it);
	    	reco::RecoToSimCollection::const_iterator matchedSim = recSimColl.find(edm::RefToBase<reco::Track>(trackRef));
	        if(matchedSim == recSimColl.end()){
	      	}
			else{
	        	const TrackingParticle* tparticle = matchedSim->val[0].first.get();
				cout<<"Sim tk pdgID: " << tparticle->pdgId() << endl;
	    	    if (tparticle->parentVertex().isNonnull() && !tparticle->parentVertex()->sourceTracks().empty()) {
					cout<<"Sim tk Mother pdgID: " << tparticle->parentVertex()->sourceTracks()[0]->pdgId() << endl;
		    	} 
	
	         	if (!tparticle->genParticle().empty()) {
	 				const HepMC::GenParticle *gp;
		            const HepMC::GenVertex *vtx = gp->production_vertex();
	        	    const HepMC::GenParticle *genMom;
	    	     	if (vtx != 0 && vtx->particles_in_size() > 0) {
	        	    	genMom = *vtx->particles_in_const_begin();
					}
			        if (genMom) {
		         		cout<<"Sim tk Mother2 pdgID: " << genMom->pdg_id() << endl;
		     	 	}
	         	}
			}
		}
	}
	*/


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
GenericAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenericAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
GenericAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenericAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenericAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenericAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenericAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void
GenericAnalyzer::ShowMyOffspring(const reco::Candidate* p, int layer){
    int ndau = p->numberOfDaughters();
    if(ndau != 0){
        for(int i = 0; i < ndau; i++){
			const reco::Candidate* MyDau = (p->daughter(i));
    		int ndaudau = MyDau->numberOfDaughters();
			for(int l = 0; l < layer; l++){ std::cout<<"   ";}
		    std::cout<< "daughter: " << MyDau->pdgId() << " / " << MyDau->status() <<" / "<<MyDau->pt()<<" / "<<ndaudau<<std::endl;
			ShowMyOffspring(MyDau, layer+1);
		}
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenericAnalyzer);

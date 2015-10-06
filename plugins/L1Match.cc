/*
 *  L1Match.cc
 *  
 *
 *  Created by Florian Scheuch on 09.04.14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 *
 *	
 *	Ghost definition?: Ghost is deltaR > 0.1 away from generating particle
 *	TODO: Have a look at the events with ghosts
 *
 */

// system include files
#include <memory>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/PhysObjectMatcher.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"

//Include Muon Propagator
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"


class L1Match : public edm::EDAnalyzer {
	public:
		explicit L1Match(const edm::ParameterSet&);
		~L1Match();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual bool equalsGenParticles(reco::GenParticle gp1, reco::GenParticle gp2);
	
	TH1D* relMuonNumber;
	TH1D* numberOfMatches;
	TH1D* numberOfMatchesGhostEvent;
	TH1D* numberOfMatchesNonGhostEvent;

	TH1D* delRallEvents;
	TH1D* delRghostEvents;
	TH1D* ptAllEvents;
	TH1D* ptGhostEvents;
	TH1D* etaAllEvents;
	TH1D* etaGhostEvents;
	TH1D* phiAllEvents;
	TH1D* phiGhostEvents;
	TH1D* pdgIdAllEvents;
	TH1D* pdgIdGhostEvents;
	TH1D* chargeAllEvents;
	TH1D* chargeGhostEvents;
	TH1D* deltaPhiAllEvents;
	TH1D* deltaPhiGhostEvents;
	TH1D* deltaThetaAllEvents;
	TH1D* deltaThetaGhostEvents;
	TH1D* relPtInAllEvents;
	TH1D* relPtInGhostEvents;
	TH1D* delVzAllEvents;
	TH1D* delVzNonGhostEvents;
	TH1D* delVzGhostEvents;

	TH1D* delRnonGhostEvents;
	TH1D* ptNonGhostEvents;
	TH1D* etaNonGhostEvents;
	TH1D* phiNonGhostEvents;
	TH1D* pdgIdNonGhostEvents;
	TH1D* chargeNonGhostEvents;
	TH1D* deltaPhiNonGhostEvents;
	TH1D* deltaThetaNonGhostEvents;
	TH1D* relPtInNonGhostEvents;

	TH2D* corr1overPtAndDeltaRinAllEvents;
	TH2D* corr1overPtAndDeltaRinNonGhostEvents;
	TH2D* corr1overPtAndDeltaRinGhostEvents;
	TH2D* corrPtAndDeltaRinAllEvents;
	TH2D* corrPtAndDeltaRinNonGhostEvents;
	TH2D* corrPtAndDeltaRinGhostEvents;
	TH2D* corrRelPtAndDeltaRinAllEvents;
	TH2D* corrRelPtAndDeltaRinNonGhostEvents;
	TH2D* corrRelPtAndDeltaRinGhostEvents;
	TH2D* corrDeltaPhiAndDeltaThetaInAllEvents;
	TH2D* corrDeltaPhiAndDeltaThetaInNonGhostEvents;
	TH2D* corrDeltaPhiAndDeltaThetaInGhostEvents;

	int distinguished;
	int notDistinguished;
	int noMuons;
	int over15Below20;
	int noEvents;
	const double GRANULARITY = .087;
};

L1Match::L1Match(const edm::ParameterSet& iConfig) {
	//now do what ever initialization is needed
}

L1Match::~L1Match() {
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}

void L1Match::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	noEvents++;

//	edm::Handle<edm::View<reco::Muon>> muonHandle;
//	iEvent.getByLabel("muons", muonHandle);

//	edm::Handle<reco::MuonCollection> muonCollectionHandle;
//	iEvent.getByLabel("muons", muonCollectionHandle);
	
	// Get GenParticleCollection
	edm::Handle<reco::GenParticleCollection> genParticleHandle;
 	iEvent.getByLabel("genParticles", genParticleHandle);
	
// 	edm::Handle<edm::View<l1extra::L1MuonParticle>> l1MuonHandle;
//	iEvent.getByLabel("l1extraParticles", l1MuonHandle);

	edm::Handle<l1extra::L1MuonParticleCollection> l1MuonCollectionHandle;
	iEvent.getByLabel("l1extraParticles", l1MuonCollectionHandle);

	
	std::cout << l1MuonCollectionHandle -> size() << std::endl;
	if(l1MuonCollectionHandle -> size() > 1)
		std::cout << l1MuonCollectionHandle -> size() << std::endl;

	//Get association with genParticle matches
//	edm::Handle<reco::GenParticleMatch> genParticleMatchHandle;
//	iEvent.getByLabel("muonMatch",genParticleMatchHandle);
//	
//	edm::Handle<edm::ValueMap<reco::CandidatePtr>> rawmatches;
//	iEvent.getByLabel("muonL1Match", rawmatches);
//
//	edm::Handle<edm::ValueMap<int>> qualities;
//	iEvent.getByLabel(edm::InputTag("muonL1Match", "quality"), qualities);

	//std::cout << "Content of Handle: " << rawmatches->size() << std::endl;

	//std::vector<pat::TriggerObjectStandAlone> match = muonCollectionHandle->at(0).triggerMatchesByFilter("l1");

	//std::cout << muonHandle->size() << " MuonHandleSize" << std::endl;
	//std::cout << l1MuonCollectionHandle->size() << " L1MuonCollectionHandleSize" << std::endl;
	//std::cout << genParticleMatchHandle->size() << " GenParticleMatchHandle" << std::endl;


	std::vector<reco::GenParticle> processed_particles; // Container for GenParticles that already have been associated
//	std::vector<edm::RefToBase<l1extra::L1MuonParticle>> processed_particles_ref;

//	std::vector<std::vector<l1extra::L1MuonParticle>> matchedVector;	

//	std::vector<int> numberOfRefs;

//	bool ghostEvent = false;

//	for(unsigned int idx = 0; idx < muonHandle->size(); idx++){
//		edm::Ptr<reco::Muon> muPtr = muonHandle->ptrAt(idx);
//		reco::CandidatePtr rawmatchptr = (*rawmatches)[muPtr];
//		if(rawmatchptr.isNonnull()){
//			//const l1extra::L1MuonParticle &l1p = dynamic_cast<const l1extra::L1MuonParticle &>(*rawmatchptr);
//			//const L1MuGMTExtendedCand &gmtCand = l1p.gmtMuonCand();
//			//bool rpc = gmtCand.isRPC();
//			//int quality = (*qualities)[muPtr];
//			//std::cout << rpc << std::endl;
//			//std::cout << l1p.pt() << " pt" << std::endl;
//			//std::cout << quality << " Quality" << std::endl;
//		}
//	}

	// Search through all reco mu
//	for(unsigned int idx = 0; idx < l1MuonHandle->size(); idx++){		
//		edm::RefToBase<l1extra::L1MuonParticle> muonRef = l1MuonHandle->refAt(idx);
//		l1extra::L1MuonParticle amuon(*(muonRef.get()));
//		
//		// Check if muon in barrel region
//		if(amuon.eta()>.8 || amuon.eta()<-.8) {
//			continue;
//		}
//		std::cout << amuon.eta() << std::endl;
//		// Just look for global muons
//		//if(!amuon.isGlobalMuon()) continue;
//		
//		// Get associated gen particle
//		reco::GenParticleRef genparticleref = (*genParticleMatchHandle)[muonRef];
//
//		// Check if genParticle is present
//		if (!genparticleref.isNonnull() || !genparticleref.isAvailable()){
//			std::cout << "Ref is null" << std::endl;
//			continue;
//		}
//		
//		const reco::GenParticle genparticle = *(genparticleref.get());
//
//		// Check if genparticle and reco mu are from the same vertex
//		//if(fabs(genparticle.vz() - amuon.vz()) > 0.01) continue;
//		//if(fabs(genparticle.vy() - amuon.vy()) > 0.01) continue;
//		//if(fabs(genparticle.vx() - amuon.vx()) > 0.01) continue;
//		
//		// Fill data
//		double del_R_reco_gen = sqrt(
//									 (genparticle.phi()-amuon.phi()) * (genparticle.phi()-amuon.phi()) +
//									 (genparticle.eta()-amuon.eta()) * (genparticle.eta()-amuon.eta())
//									 ); // Distance of reco mu and matching gen mu
//		del_R_reco_gen = fmod(del_R_reco_gen, TMath::Pi()); //Normalize from 0 to Pi
//		
//		delRallEvents -> Fill(del_R_reco_gen);
//		delVzAllEvents -> Fill(amuon.vz() - genparticle.vz());
//		ptAllEvents -> Fill(amuon.pt());
//		etaAllEvents -> Fill(amuon.eta());
//		phiAllEvents -> Fill(amuon.phi());
//		pdgIdAllEvents -> Fill(amuon.pdgId());
//		chargeAllEvents -> Fill(amuon.charge());
//		
//		double del_phi_all = sqrt((genparticle.phi()-amuon.phi()) * (genparticle.phi()-amuon.phi())); // Distance of reco mu and matching gen mu
//		del_phi_all = fmod(del_phi_all, TMath::Pi()); //Normalize from 0 to Pi
//		deltaPhiAllEvents -> Fill(del_phi_all);
//		
//		double del_theta_all = sqrt((genparticle.theta()-amuon.theta()) * (genparticle.theta()-amuon.theta())); // Distance of reco mu and matching gen mu
//		del_theta_all = fmod(del_theta_all, TMath::Pi()); //Normalize from 0 to Pi
//		deltaThetaAllEvents -> Fill(del_theta_all);
//
//		relPtInAllEvents -> Fill(amuon.pt() / genparticle.pt());
//
//		corr1overPtAndDeltaRinAllEvents -> Fill(1./amuon.pt(), del_R_reco_gen);
//		corrPtAndDeltaRinAllEvents -> Fill(amuon.pt(), del_R_reco_gen);
//		corrRelPtAndDeltaRinAllEvents -> Fill(amuon.pt() / genparticle.pt(), del_R_reco_gen);
//		corrDeltaPhiAndDeltaThetaInAllEvents -> Fill(del_phi_all, del_theta_all);
//
//		bool ghost = false;
//
//		// Go through handled genParticles, if genParticle was handled before -> this reco muon or the first muon is a muon ghost
//		for(unsigned int jdx = 0; jdx < processed_particles.size(); jdx++){
//			reco::GenParticle gp = processed_particles.at(jdx);
//			if(equalsGenParticles(gp, genparticle)){
//				numberOfRefs[jdx]++;
//				matchedVector[jdx].push_back(amuon);
//				ghost = true;
//				ghostEvent = true;
//				std::cout << "Ghost Event!" << std::endl;
//			}
//		}
//		if(!ghost){
//			// Fill the processed particles in the vector
//			std::vector<l1extra::L1MuonParticle> v;
//			v.push_back(amuon);
//			matchedVector.push_back(v);
//			processed_particles.push_back(genparticle);
//			processed_particles_ref.push_back(muonRef);
//			numberOfRefs.push_back(1); //Count 1 reference
//		}
//	}
//
//	//noMuons += matchedVector.size();
//	relMuonNumber -> Fill(l1MuonHandle->size() - processed_particles.size());
//	//std::cout << "Number of Refs: " << numberOfRefs.size() << std::endl;
//	for(unsigned int idx = 0; idx < numberOfRefs.size(); idx++){
//		if(numberOfRefs.at(idx) == 4){
//			std::cout << "BX:" << iEvent.id() << ": " << numberOfRefs.at(idx) << std::endl;
//		}
//		//std::cout << "Number of Refs (" << idx << "): " << numberOfRefs.at(idx) << std::endl;
//		numberOfMatches->Fill(numberOfRefs.at(idx));
//	}
//	if(ghostEvent){
//		for(unsigned int idx = 0; idx < numberOfRefs.size(); idx++){
//			numberOfMatchesGhostEvent->Fill(numberOfRefs.at(idx));
//		}
//	} else {
//		for(unsigned int idx = 0; idx < numberOfRefs.size(); idx++){
//			numberOfMatchesNonGhostEvent->Fill(numberOfRefs.at(idx));
//		}
//	}
//
//	noMuons += genParticleHandle->size();
//
//
//	// Go around again and fill in ghosts and non ghost events
//	for(unsigned int idx = 0; idx < l1MuonHandle->size(); idx++){		
//		edm::RefToBase<l1extra::L1MuonParticle> muonRef = l1MuonHandle->refAt(idx);
//		l1extra::L1MuonParticle amuon(*(muonRef.get()));
//		
//		// Check if muon in barrel region
//		if(amuon.eta()>.8 && amuon.eta()<-.8) continue;
//		
//		// Just look for global muons
//		//if(!amuon.isGlobalMuon()) continue;
//		
//		// Get associated gen particle
//		reco::GenParticleRef genparticleref = (*genParticleMatchHandle)[muonRef];
//
//		// Check if genParticle is present
//		if (!genparticleref.isNonnull() || !genparticleref.isAvailable()){
//			std::cout << "Ref is null" << std::endl;
//			continue;
//		}
//		
//		const reco::GenParticle genparticle = *(genparticleref.get());
//		
//		double del_R_reco_gen = sqrt(
//									 (genparticle.phi()-amuon.phi()) * (genparticle.phi()-amuon.phi()) +
//									 (genparticle.eta()-amuon.eta()) * (genparticle.eta()-amuon.eta())
//									 ); // Distance of reco mu and matching gen mu
//		del_R_reco_gen = fmod(del_R_reco_gen, TMath::Pi()); //Normalize from 0 to Pi
//
//		double del_phi_all = sqrt((genparticle.phi()-amuon.phi()) * (genparticle.phi()-amuon.phi())); // Distance of reco mu and matching gen mu
//		del_phi_all = fmod(del_phi_all, TMath::Pi()); //Normalize from 0 to Pi
//		
//		double del_theta_all = sqrt((genparticle.theta()-amuon.theta()) * (genparticle.theta()-amuon.theta())); // Distance of reco mu and matching gen mu
//		del_theta_all = fmod(del_theta_all, TMath::Pi()); //Normalize from 0 to Pi
//
//		if(!ghostEvent){
//			delRnonGhostEvents -> Fill(del_R_reco_gen);
//			delVzNonGhostEvents -> Fill(amuon.vz() - genparticle.vz());
//			ptNonGhostEvents -> Fill(amuon.pt());
//			etaNonGhostEvents -> Fill(amuon.eta());
//			phiNonGhostEvents -> Fill(amuon.phi());
//			pdgIdNonGhostEvents -> Fill(amuon.pdgId());
//			chargeNonGhostEvents -> Fill(amuon.charge());
//			deltaPhiNonGhostEvents -> Fill(del_phi_all);
//			deltaThetaNonGhostEvents -> Fill(del_theta_all);
//			relPtInNonGhostEvents -> Fill(amuon.pt() / genparticle.pt());
//					
//			corr1overPtAndDeltaRinNonGhostEvents -> Fill(1./amuon.pt(), del_R_reco_gen);
//			corrPtAndDeltaRinNonGhostEvents -> Fill(amuon.pt(), del_R_reco_gen);
//			corrRelPtAndDeltaRinNonGhostEvents -> Fill(amuon.pt() / genparticle.pt(), del_R_reco_gen);
//			corrDeltaPhiAndDeltaThetaInNonGhostEvents -> Fill(del_phi_all, del_theta_all);
//		} else {
//			delRghostEvents -> Fill(del_R_reco_gen);
//			delVzGhostEvents -> Fill(amuon.vz() - genparticle.vz());
//			ptGhostEvents -> Fill(amuon.pt());
//			etaGhostEvents -> Fill(amuon.eta());
//			phiGhostEvents -> Fill(amuon.phi());
//			pdgIdGhostEvents -> Fill(amuon.pdgId());
//			chargeGhostEvents -> Fill(amuon.charge());
//			deltaPhiGhostEvents -> Fill(del_phi_all);
//			deltaThetaGhostEvents -> Fill(del_theta_all);
//			relPtInGhostEvents -> Fill(amuon.pt() / genparticle.pt());
//				
//			corr1overPtAndDeltaRinGhostEvents -> Fill(1./amuon.pt(), del_R_reco_gen);
//			corrPtAndDeltaRinGhostEvents -> Fill(amuon.pt(), del_R_reco_gen);
//			corrRelPtAndDeltaRinGhostEvents -> Fill(amuon.pt() / genparticle.pt(), del_R_reco_gen);
//			corrDeltaPhiAndDeltaThetaInGhostEvents -> Fill(del_phi_all, del_theta_all);
//		}
//	}
	bool trig23 = false;
	bool trig15 = false;
	int counting = 0;
	for(unsigned int idx = 0; idx < l1MuonCollectionHandle->size(); idx++){
		if(l1MuonCollectionHandle->at(idx).pt() > 20){
			trig23 = true;
			break;
		}
		if(l1MuonCollectionHandle->at(idx).pt() > 10){
			counting++;
		}
		if(l1MuonCollectionHandle->at(idx).pt() > 15){
			trig15 = true;
		}
	}
	
	if(counting >= 2 && !trig23 && trig15){
		over15Below20++;
	}
		
		
//		//bool trig18 = false;
//		//bool trig10 = false;
//		for(unsigned int jdx = 0; jdx < matchedVector[idx].size(); jdx++){
//			if(matchedVector[idx][jdx].pt() > 25) {
//				trig23 = true;
//				break;
//			}
//		//	if(matchedVector[idx][jdx].pt() > 10 && trig10){
//		//		trig18 = true;
//		//		break;
//		//	}
////			if(matchedVector[idx][jdx].pt() > 10){
////				counting++;
////			}
//		}
//		if(counting >= 2 && !trig23)
//			over15Below20++;
	
	
//	for(unsigned int idx = 0; idx < matchedVector.size(); idx++){
		for(unsigned int jdx = 0; jdx < l1MuonCollectionHandle->size(); jdx++){
			for(unsigned int kdx = 0; kdx < l1MuonCollectionHandle->size(); kdx++){
				if(jdx == kdx) continue;
				l1extra::L1MuonParticle one =  l1MuonCollectionHandle->at(jdx);
				l1extra::L1MuonParticle two =  l1MuonCollectionHandle->at(kdx);
				if(((int)(one.phi()/GRANULARITY) == (int)(two.phi()/GRANULARITY)) && ((int)(one.eta()/GRANULARITY) == (int)(two.eta()/GRANULARITY))){
					notDistinguished++;
				} else {
					distinguished++;
				}
			}
		}
	
//	}


}

void L1Match::beginJob() {
	const double PI  = 3.141592653589793238462;
	
	distinguished = 0;
	notDistinguished = 0;
	noMuons = 0;
	over15Below20 = 0;
	noEvents = 0;

	edm::Service<TFileService> fs;
	relMuonNumber = fs -> make<TH1D>("NumberOfL1 muons - matches", "NumberOfL1 muons - matches", 20, -0.5, 19.5);
	numberOfMatches = fs -> make<TH1D>("NumberOf L1 muon matches", "NumberOf L1 muon matches", 20, -.5, 19.5);
	numberOfMatchesGhostEvent = fs -> make<TH1D>("NumberOf L1 muon matches in ghost events", "NumberOf L1 muon matches in ghost events", 20, -.5, 19.5);
	numberOfMatchesNonGhostEvent = fs -> make<TH1D>("NumberOf L1 muon matches in non ghost events", "NumberOf L1 muon matches in non ghost events", 20, -.5, 19.5);

	//TFileDirectory dir_number_of_ghosts = fs->mkdir("Number_of_ghosts_sum");
	delRallEvents = fs -> make<TH1D>("Delta R between L1Mu and GEN mu in all events", "Delta R between L1Mu and GEN mu in all events", 1000, 0, 2*PI);
	delRnonGhostEvents = fs -> make<TH1D>("Delta R between L1Mu and GEN mu in nonGhost events", "Delta R between L1Mu and GEN mu in nonGhost events", 1000, 0, 2*PI);
	delRghostEvents = fs -> make<TH1D>("Delta R between L1Mu and GEN mu in ghost events", "Delta R between L1Mu and GEN mu in ghost events", 1000, 0, 2*PI);

	delVzAllEvents = fs -> make<TH1D>("Delta Vz between L1Mu and GEN mu in all events", "Delta Vz between L1Mu and GEN mu in all events", 1000, -10, 10);
	delVzNonGhostEvents = fs -> make<TH1D>("Delta Vz between L1Mu and GEN mu in non ghost events", "Delta Vz between L1Mu and GEN mu in non ghost events", 1000, -10, 10);
	delVzGhostEvents = fs -> make<TH1D>("Delta Vz between L1Mu and GEN mu in ghost events", "Delta Vz between L1Mu and GEN mu in ghost events", 1000, -10, 10);
	
	ptAllEvents = fs -> make<TH1D>("pt of L1Mu in all events", "pt of L1Mu in all events", 1000, 0, 1000);
	ptNonGhostEvents = fs -> make<TH1D>("pt of L1Mu in non ghost events", "pt of L1Mu in non ghost events", 1000, 0, 1000);
	ptGhostEvents = fs -> make<TH1D>("pt of L1Mu in ghost events", "pt of L1Mu in ghost events", 1000, 0, 1000);
	
	etaAllEvents  = fs -> make<TH1D>("eta of L1Mu in all events", "eta of L1Mu in all events", 1000, 0, 10);
	etaNonGhostEvents = fs -> make<TH1D>("eta of L1Mu in non ghost events", "eta of L1Mu in non ghost events", 1000, 0, 10);
	etaGhostEvents = fs -> make<TH1D>("eta of L1Mu in ghost events", "eta of L1Mu in ghost events", 1000, 0, 10);
	
	phiAllEvents = fs -> make<TH1D>("phi of L1Mu in all events", "phi of L1Mu in all events", 1000, -2*PI, 2*PI);
	phiNonGhostEvents = fs -> make<TH1D>("phi of L1Mu in non ghost events", "phi of L1Mu in non ghost events", 1000, -2*PI, 2*PI);
	phiGhostEvents = fs -> make<TH1D>("phi of L1Mu in ghost events", "phi of L1Mu in ghost events", 1000, -2*PI, 2*PI);
	
	pdgIdAllEvents = fs -> make<TH1D>("pdgId of L1Mu in all events", "pdgId of L1Mu in all events", 1000, -500, 499);
	pdgIdNonGhostEvents  = fs -> make<TH1D>("pdgId of L1Mu in non ghost events", "pdgId of L1Mu in non ghost events", 5, -2, 2);
	pdgIdGhostEvents  = fs -> make<TH1D>("pdgId of L1Mu in ghost events", "pdgId of L1Mu in ghost events", 5, -2, 2);
	
	chargeAllEvents = fs -> make<TH1D>("charge of L1Mu in all events", "charge of L1Mu in all events", 5, -2, 2);
	chargeNonGhostEvents = fs -> make<TH1D>("charge of L1Mu in non ghost events", "charge of L1Mu in non ghost events", 5, -2, 2);
	chargeGhostEvents = fs -> make<TH1D>("charge of L1Mu in ghost events", "charge of L1Mu in ghost events", 5, -2, 2);
	
	deltaPhiAllEvents = fs -> make<TH1D>("Delta phi between L1Mu and GEN mu in all events", "Delta phi between L1Mu and GEN mu in all events", 1000, 0, 2*PI);
	deltaPhiNonGhostEvents = fs -> make<TH1D>("Delta phi between L1Mu and GEN mu in non ghost events", "Delta phi between L1Mu and GEN mu in non ghost events", 1000, 0, 2*PI);
	deltaPhiGhostEvents = fs -> make<TH1D>("Delta phi between L1Mu and GEN mu in ghost events", "Delta phi between L1Mu and GEN mu in ghost events", 1000, 0, 2*PI);
	
	deltaThetaAllEvents = fs -> make<TH1D>("Delta theta between L1Mu and GEN mu in all events", "Delta theta between L1Mu and GEN mu in all events", 1000, 0, 10);
	deltaThetaNonGhostEvents = fs -> make<TH1D>("Delta theta between L1Mu and GEN mu in non ghost events", "Delta theta between L1Mu and GEN mu in non ghost events", 1000, 0, 10);
	deltaThetaGhostEvents = fs -> make<TH1D>("Delta theta between L1Mu and GEN mu in ghost events", "Delta theta between L1Mu and GEN mu in ghost events", 1000, 0, 10);
	
	relPtInAllEvents = fs -> make<TH1D>("Rel Pt between L1Mu and GEN mu in all events", "Rel Pt between L1Mu and GEN mu in all events", 1000, 0, 10);
	relPtInNonGhostEvents = fs -> make<TH1D>("Rel Pt between L1Mu and GEN mu in non ghost events", "Rel Pt between L1Mu and GEN mu in non ghost events", 1000, 0, 10);
	relPtInGhostEvents = fs -> make<TH1D>("Rel Pt between L1Mu and GEN mu in ghost events", "Rel Pt between L1Mu and GEN mu in ghost events", 1000, 0, 10);

	corr1overPtAndDeltaRinAllEvents = fs -> make<TH2D>("2D 1/Pt vs Delta R in all events", "2D 1/Pt vs Delta R in all events", 200, 0, 1, 250, 0, 2*PI);
	corr1overPtAndDeltaRinNonGhostEvents = fs -> make<TH2D>("2D 1/Pt vs Delta R in non ghost events", "2D 1/Pt vs Delta R in non ghost events", 200, 0, 1, 250, 0, 2*PI);
	corr1overPtAndDeltaRinGhostEvents = fs -> make<TH2D>("2D 1/Pt vs Delta R in ghost events", "2D 1/Pt vs Delta R in ghost events", 200, 0, 1, 250, 0, 2*PI);
	
	corrPtAndDeltaRinAllEvents = fs -> make<TH2D>("2D Pt vs Delta R in all events", "2D Pt vs Delta R in all events", 200, 0, 200, 250, 0, 2*PI);
	corrPtAndDeltaRinNonGhostEvents = fs -> make<TH2D>("2D Pt vs Delta R in non ghost events", "2D Pt vs Delta R in non ghost events", 200, 0, 200, 250, 0, 2*PI);
	corrPtAndDeltaRinGhostEvents = fs -> make<TH2D>("2D Pt vs Delta R in ghost events", "2D Pt vs Delta R in ghost events", 200, 0, 200, 250, 0, 2*PI);
	
	corrRelPtAndDeltaRinAllEvents = fs -> make<TH2D>("2D Rel Pt vs Delta R in all events", "2D Rel Pt vs Delta R in all events", 200, 0, 10, 250, 0, 2*PI);
	corrRelPtAndDeltaRinNonGhostEvents = fs -> make<TH2D>("2D Rel Pt vs Delta R in non ghost events", "2D Rel Pt vs Delta R in non ghost events", 200, 0, 10, 250, 0, 2*PI);
	corrRelPtAndDeltaRinGhostEvents = fs -> make<TH2D>("2D Rel Pt vs Delta R in ghost events", "2D Rel Pt vs Delta R in ghost events", 200, 0, 10, 250, 0, 2*PI);
	
	corrDeltaPhiAndDeltaThetaInAllEvents = fs -> make<TH2D>("2D Delta Phi vs Delta Theta in all events", "2D Delta Phi vs Delta Theta in all events", 200, 0, PI, 250, 0, 10);
	corrDeltaPhiAndDeltaThetaInNonGhostEvents = fs -> make<TH2D>("2D Delta Phi vs Delta Theta in non ghost events", "2D Delta Phi vs Delta Theta in non ghost events", 200, 0, PI, 250, 0, 10);
	corrDeltaPhiAndDeltaThetaInGhostEvents = fs -> make<TH2D>("2D Delta Phi vs Delta Theta in ghost events", "2D Delta Phi vs Delta Theta in ghost events", 200, 0, PI, 250, 0, 10);


}

void L1Match::endJob(){
	std::cout << "Not distinguished = " << notDistinguished << std::endl;
	std::cout << "Distinguished = " << distinguished << std::endl;
	std::cout << "Distinguish ratio = " << distinguished * 1. / (notDistinguished + distinguished) << std::endl;
	std::cout << "No of ghosts = " << (distinguished+notDistinguished) << std::endl;
	std::cout << "No of muons = " << noMuons << std::endl;
	std::cout << "Additional double mu trigger = " << over15Below20 << std::endl;
	std::cout << "No of events = " << noEvents << std::endl;
}


//Checks if the GenParticles are equal in condition of charge, p_t
bool L1Match::equalsGenParticles(reco::GenParticle gp1, reco::GenParticle gp2){
	if(gp1.pt() != gp2.pt()){
		return false;
	}
	if(gp1.charge() != gp2.charge()){
		return false;
	}
	if(gp1.px() != gp2.px()){
		return false;
	}
	if(gp1.py() != gp2.py()){
		return false;
	}
	if(gp1.vx() != gp2.vx()){
		return false;
	}
	return true;
}

// ------------ method called when starting to processes a run  ------------
void L1Match::beginRun(edm::Run const&, edm::EventSetup const& iSetup){
	//propagatorToMuon.init(iSetup);
}

// ------------ method called when ending the processing of a run  ------------
void L1Match::endRun(edm::Run const&, edm::EventSetup const&){
}

// ------------ method called when starting to processes a luminosity block  ------------
void L1Match::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
	std::cout << "new Lumi Block" << std::endl;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void L1Match::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1Match::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1Match);


// -*- C++ -*-
//
// Package:    HLTForDiHiggs/HLTForDiHiggsAnalyzer
// Class:      HLTForDiHiggsAnalyzer
//
/**\class HLTForDiHiggsAnalyzer HLTForDiHiggsAnalyzer.cc HLTForDiHiggs/HLTForDiHiggsAnalyzer/plugins/HLTForDiHiggsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marina Kolosova
//         Created:  Wed, 25 Jan 2023 10:23:45 GMT
//
//

// system include files
#include <memory>
#include <random>
#include <queue>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <array>
#include <iostream>
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Trigger include files
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// Physics Objects include files
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Common/interface/ValueMap.h"

// Physics Objects - MiniAOD
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class HLTForDiHiggsAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HLTForDiHiggsAnalyzer(const edm::ParameterSet&);
  ~HLTForDiHiggsAnalyzer() override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void initiateHistVariables();
  
  TTree* initiateOutTree(edm::Service<TFileService> fs);
  
  void beginJob() override;
  double getdR(const double eta1, const double eta2, const double phi1, const double phi2);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::TriggerResults>    trigresultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent>  trigsummaryToken_;
  
  edm::EDGetTokenT<std::vector<reco::Vertex> > vertexToken_;  
  edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken_;
  edm::EDGetTokenT<std::vector<pat::Jet> > fatjetsToken_;
  edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron> > eleToken_;
  
  edm::EDGetToken genJetsToken_;
  edm::EDGetTokenT<reco::GenJetCollection>      genFatjetsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  
  const bool isMC_;
  const bool isSignal_;
  unsigned int run_;
  unsigned int lumi_;
  unsigned int evt_;
  bool isData_;
  float rho_;
  
  int nEvents;
  const int maxResults_ = 500;
  bool firstEvent_;
  
  unsigned int *passtrig_;
  
  // Physics objects
  
  // Muons
  /*
  int nmuons;
  float *muon_pt;
  float *muon_eta;
  float *muon_phi;
  float *muon_mass;
  float *muon_charge;  
  float *muon_IsLoose;
  float *muon_IsMedium;
  float *muon_relIso;
  
  // Electrons
  int nelectrons;
  float *electron_pt;
  float *electron_eta;
  float *electron_phi;
  float *electron_mass;
  float *electron_reliso;
  */
  
  // b-quarks from Higgs
  int nbquarks;
  float *bquark_pt;
  float *bquark_eta;
  float *bquark_phi;
  float *bquark_mass;

  int ngenjets;
  float *genjet_pt;
  float *genjet_eta;
  float *genjet_phi;
  float *genjet_mass;
  
  int ngenfatjets;
  float *genfatjet_pt;
  float *genfatjet_eta;
  float *genfatjet_phi;
  float *genfatjet_mass;
  
  // AK4 Jets
  int njets;
  float *jet_pt;
  float *jet_eta;
  float *jet_phi;
  float *jet_mass;
  float *jet_deepjet_probb;
  float *jet_deepjet_probbb;
  float *jet_deepjet_problepb;
  float *jet_deepjet;
  float *jet_pnet_probb;
  float *jet_pnet_probbb;
  float *jet_pnet_probc;
  float *jet_pnet_probcc;
  float *jet_pnet_probuds;
  float *jet_pnet_probg;
  float *jet_pnet_probundef;
  float *jet_pnet_probpu;
  float *jet_pnet_BvsAll;
  
  // AK8 Jets
  int nfatjets;
  float *fatjet_pt;
  float *fatjet_eta;
  float *fatjet_phi;
  float *fatjet_mass;
  float *fatjet_massSD;
  float *fatjet_PNetXbb;
  
  TTree* outTree_;
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
HLTForDiHiggsAnalyzer::HLTForDiHiggsAnalyzer(const edm::ParameterSet& iConfig)
  :
  // Trigger Results
  trigresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltresults"))),
  trigsummaryToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("trigSummary"))),
  vertexToken_(consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"))),
  jetsToken_(consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsPuppi"))),
  fatjetsToken_(consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsAK8"))),
  muonToken_(consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"))),
  eleToken_(consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons"))),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJets"))),
  genFatjetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJetsAK8"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
  isMC_(iConfig.getParameter<bool>("isMC")),
  isSignal_(iConfig.getParameter<bool>("isSignal")),
  run_(0),
  lumi_(0),
  evt_(0),
  maxResults_(500)
{
  firstEvent_ = true;
  initiateHistVariables();
  
  edm::Service<TFileService> fs;
  outTree_=initiateOutTree(fs);
}

HLTForDiHiggsAnalyzer::~HLTForDiHiggsAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

double HLTForDiHiggsAnalyzer::getdR(const double eta1, const double eta2, const double phi1, const double phi2)
{
  const double pi = 3.14159265358979323846;
  double dPhi = std::abs(std::abs(std::abs(phi1 - phi2) - pi) - pi);
  double dEta = eta1-eta2;
  double dR2 = dPhi*dPhi + dEta*dEta;
  double dR  = sqrt(dR2);
  return dR;
}

// ------------ method called for each event  ------------
void HLTForDiHiggsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Event
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  run_  = iEvent.id().run();
  lumi_ = iEvent.luminosityBlock();
  evt_  = iEvent.id().event();
  isData_ = iEvent.isRealData();
  
  std::cout<<" "<<std::endl;
  std::cout<<"  Event = "<<evt_<<std::endl;
  
  //edm::Handle<reco::VertexCollection> vertices;
  //iEvent.getByToken(vertexToken_, vertices);
  //if (vertices->empty()) return;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // HLT Paths
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  edm::Handle<edm::TriggerResults> trigresults;
  iEvent.getByToken(trigresultsToken_, trigresults);
  
  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*trigresults);
  for(unsigned int itrig = 0; itrig < trigresults->size(); ++itrig) {
    TString trigname = triggerNames.triggerName(itrig);
    
    if (trigname == "HLTriggerFirstPath" || trigname == "HLTriggerFinalPath" || trigname == "p") continue;
    if(firstEvent_) 
      {
	std::string hltname = triggerNames.triggerName(itrig);
	std::vector<std::string> results;
	boost::split(results, hltname, [](char c){return c == '_';});
	std::string newName;
	for (unsigned int i=0; i<results.size()-1; i++){
	  if (i == results.size()-2){
	    newName +=results.at(i)+"_v";
	  }
	  else{
	    newName+=results.at(i)+"_";
	  }
	}
	TString myNewName = newName;
	if (results.size() == 0){
	  outTree_->Branch(trigname, &passtrig_[itrig], trigname+"/i");
	}
	else{
	  outTree_->Branch(myNewName, &passtrig_[itrig], myNewName+"/i");
	}
      }
    bool accept = trigresults->accept(itrig);
    if(accept) passtrig_[itrig] = 1;
    else       passtrig_[itrig] = 0;
    std::cout << "trigName = "<<trigname<<"     accept? "<<accept<<std::endl;
  }
  if(firstEvent_) firstEvent_ = false;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Offline Selections
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);
  
  std::vector<reco::GenParticle> HiggsBosons_LastCopy;
  std::vector<reco::GenParticle> BQuarks_FirstCopy;
  
  nbquarks = 0;
  for(size_t i = 0; i < genParticles->size(); ++i) {
    const reco::GenParticle & gen = genParticles->at(i);
    
    int pdgID       = std::abs(gen.pdgId());
    int status      = gen.status();
    bool isLastCopy = gen.isLastCopy();
    unsigned int nDaus = gen.numberOfDaughters();
    // b-quarks
    if (pdgID == 5)
      {
	const reco::Candidate * mom = gen.mother();
	int mumID = std::abs(mom->pdgId());
	if (mumID == 25)
	  {
	    BQuarks_FirstCopy.push_back(gen);
	    
	    bquark_pt[nbquarks] = gen.pt();
	    bquark_eta[nbquarks] = gen.eta();
	    bquark_phi[nbquarks] = gen.phi();
	    bquark_mass[nbquarks] = gen.mass();
	    nbquarks++;
	  }
      }
    if (pdgID != 25) continue;
    if (isLastCopy) HiggsBosons_LastCopy.push_back(gen);
  } // Closes loop over all gen-particles
  
  edm::Handle<std::vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsToken_, genJets);
  ngenjets = 0;
  for (const auto& genjet : *genJets)
    {
      double pt = genjet.pt();
      double eta = genjet.eta();
      double phi = genjet.phi();
      double mass = genjet.mass();
      
      genjet_pt[ngenjets]   = pt;
      genjet_eta[ngenjets]  = eta;
      genjet_phi[ngenjets]  = phi;
      genjet_mass[ngenjets] = mass;
      ngenjets++;
    }
  
  edm::Handle<std::vector<reco::GenJet> > genFatjets;
  iEvent.getByToken(genFatjetsToken_, genFatjets);
  
  ngenfatjets = 0;
  for (const auto& genfatjet : *genFatjets)
    {
      double pt = genfatjet.pt();
      double eta = genfatjet.eta();
      double phi = genfatjet.phi();
      double mass = genfatjet.mass();
      
      genfatjet_pt[ngenfatjets] = pt;
      genfatjet_eta[ngenfatjets] = eta;
      genfatjet_phi[ngenfatjets] = phi;
      genfatjet_mass[ngenfatjets] = mass;
      ngenfatjets++;
    }
  
  //=========================================================================================================================================================================================
  //                                                                                  Jet Selection
  //=========================================================================================================================================================================================
  njets = 0;
  edm::Handle< std::vector<pat::Jet> > jets;
  iEvent.getByToken(jetsToken_, jets);
  for( std::vector<pat::Jet>::const_iterator jet = (*jets).begin(); jet != (*jets).end(); jet++)
    {
      double pt = jet->pt();
      double eta = jet->eta();
      double phi = jet->phi();
      double mass = jet->mass();
      
      double deepjet_probb  = jet->bDiscriminator("pfDeepFlavourJetTags:probb");
      double deepjet_probbb = jet->bDiscriminator("pfDeepFlavourJetTags:probbb");
      double deepjet_problepb = jet->bDiscriminator("pfDeepFlavourJetTags:problepb");
      double deepjet = deepjet_probb + deepjet_probbb + deepjet_problepb;
      
      double pnet_probb     = jet->bDiscriminator("pfParticleNetAK4JetTags:probb");
      double pnet_probbb    = jet->bDiscriminator("pfParticleNetAK4JetTags:probbb");
      double pnet_probc     = jet->bDiscriminator("pfParticleNetAK4JetTags:probc");
      double pnet_probcc    = jet->bDiscriminator("pfParticleNetAK4JetTags:probcc");
      double pnet_probuds   = jet->bDiscriminator("pfParticleNetAK4JetTags:probuds");
      double pnet_probg     = jet->bDiscriminator("pfParticleNetAK4JetTags:probg");
      double pnet_probundef = jet->bDiscriminator("pfParticleNetAK4JetTags:probundef");
      double pnet_probpu    = jet->bDiscriminator("pfParticleNetAK4JetTags:probpu");
      
      double pnet_BvsAll = jet->bDiscriminator("pfParticleNetAK4DiscriminatorsJetTags:BvsAll");
      
      if (pt < 20.0) continue;
      if (std::abs(eta) > 2.5) continue;
      
      std::cout << "AK4 jet "<<njets<<"   pt="<<pt<<"   eta="<<eta<<"   phi="<<phi<<"   probb="<<pnet_probb<<"   probbb="<<pnet_probbb<<"   BvsAll="<<pnet_BvsAll<<"   DeepJet = "<<deepjet<<std::endl;
      
      //The next following lines just remove e/mu from (semi)leptonic ttbar. 
      //if( jet->muonEnergyFraction() >0.7)continue;
      //if( jet->electronEnergyFraction() >0.7)continue;
      
      jet_pt[njets]   = pt;
      jet_eta[njets]  = eta;
      jet_phi[njets]  = phi;
      jet_mass[njets] = mass;
      
      jet_deepjet_probb[njets]    = deepjet_probb;
      jet_deepjet_probbb[njets]   = deepjet_probbb;
      jet_deepjet_problepb[njets] = deepjet_problepb;
      jet_deepjet[njets]          = deepjet;
      jet_pnet_probb[njets]       = pnet_probb;
      jet_pnet_probbb[njets]      = pnet_probbb;
      jet_pnet_probc[njets]       = pnet_probc;
      jet_pnet_probcc[njets]      = pnet_probcc;
      jet_pnet_probuds[njets]     = pnet_probuds;
      jet_pnet_probg[njets]       = pnet_probg;
      jet_pnet_probundef[njets]   = pnet_probundef;
      jet_pnet_probpu[njets]      = pnet_probpu;
      jet_pnet_BvsAll[njets]      = pnet_BvsAll;
      
      njets++;
    }
  std::cout << "Number of jets = "<<njets<<std::endl;
  
  //=========================================================================================================================================================================================
  //                                                                               Fatjet Selection
  //=========================================================================================================================================================================================
  nfatjets = 0;
  edm::Handle<std::vector<pat::Jet> > fatjets;
  iEvent.getByToken(fatjetsToken_, fatjets);
  for(std::vector<pat::Jet>::const_iterator fatjet = (*fatjets).begin(); fatjet != (*fatjets).end(); fatjet++)
    {
      double pt = fatjet->pt();
      double eta = fatjet->eta();
      double phi = fatjet->phi();
      double mass = fatjet->mass();
      
      double pnet_Xbb = fatjet->bDiscriminator("pfMassDecorrelatedParticleNetDiscriminatorsJetTags:XbbvsQCD");
      double massSD   = fatjet->userFloat("ak8PFJetsPuppiSoftDropMass");
      
      if (pt < 180) continue;
      if (std::abs(eta) > 2.5) continue;

      std::cout << "AK8 jet "<<nfatjets<<"  pt="<<pt<<"   eta="<<eta<<"   phi="<<phi<<"   mass="<<mass<<"  PNet Xbb="<<pnet_Xbb<<"  massSD ="<<massSD<<std::endl;
      
      fatjet_pt[nfatjets] = pt;
      fatjet_eta[nfatjets] = eta;
      fatjet_phi[nfatjets] = phi;
      fatjet_mass[nfatjets] = mass;
      fatjet_massSD[nfatjets] = massSD;
      fatjet_PNetXbb[nfatjets] = pnet_Xbb;
      nfatjets++;
    }
  
  std::cout << "Number of fatjets = "<<nfatjets<<std::endl;
  
  // Fill the output tree 
  outTree_->Fill();
  }
  
void HLTForDiHiggsAnalyzer::initiateHistVariables()
{
  passtrig_ = new unsigned int[maxResults_];
  nEvents = 0;

  // Muons
  /*
  nmuons = 0;
  muon_pt   = new float[maxResults_];
  muon_eta  = new float[maxResults_];
  muon_phi  = new float[maxResults_];
  muon_mass = new float[maxResults_];
  muon_charge = new float[maxResults_];
  muon_IsLoose = new float[maxResults_];
  muon_IsMedium = new float[maxResults_];
  muon_relIso  = new float[maxResults_];

  // Electron
  nelectrons      = 0;
  electron_pt     = new float[maxResults_];
  electron_eta    = new float[maxResults_];
  electron_phi    = new float[maxResults_];
  electron_mass   = new float[maxResults_];
  electron_reliso = new float[maxResults_];
  */
  
  // Gen quarks
  nbquarks    = 0;
  bquark_pt   = new float[maxResults_];
  bquark_eta  = new float[maxResults_];
  bquark_phi  = new float[maxResults_];
  bquark_mass = new float[maxResults_];
  
  // AK4 Gen-jets
  ngenjets = 0;
  genjet_pt = new float[maxResults_];
  genjet_eta = new float[maxResults_];
  genjet_phi = new float[maxResults_];
  genjet_mass = new float[maxResults_];
  
  ngenfatjets = 0;
  genfatjet_pt = new float[maxResults_];
  genfatjet_eta = new float[maxResults_];
  genfatjet_phi = new float[maxResults_];
  genfatjet_mass = new float[maxResults_];
  
  // AK4 Jets
  njets = 0;
  jet_pt          = new float[maxResults_];
  jet_eta         = new float[maxResults_];
  jet_phi         = new float[maxResults_];
  jet_mass        = new float[maxResults_];
  jet_deepjet_probb  = new float[maxResults_];
  jet_deepjet_probbb = new float[maxResults_];
  jet_deepjet_problepb = new float[maxResults_];
  jet_deepjet = new float[maxResults_];
  jet_pnet_probb = new float[maxResults_];
  jet_pnet_probbb = new float[maxResults_];
  jet_pnet_probc = new float[maxResults_];
  jet_pnet_probcc= new float[maxResults_];
  jet_pnet_probuds= new float[maxResults_];
  jet_pnet_probg= new float[maxResults_];
  jet_pnet_probundef= new float[maxResults_];
  jet_pnet_probpu= new float[maxResults_];
  jet_pnet_BvsAll= new float[maxResults_];
  
  // AK8 Jets
  nfatjets = 0;
  fatjet_pt  = new float[maxResults_];
  fatjet_eta = new float[maxResults_];
  fatjet_phi = new float[maxResults_];
  fatjet_mass = new float[maxResults_];
  fatjet_massSD = new float[maxResults_];
  fatjet_PNetXbb = new float[maxResults_];
}
 
TTree* HLTForDiHiggsAnalyzer::initiateOutTree(edm::Service<TFileService> fs)
{
  std::cout << "\n ***** initiateOutTree"<<std::endl;
  
  TTree* outTree = fs->make<TTree>("HLTforHH","");
  
  outTree->Branch("run",  &run_,  "run/i");
  outTree->Branch("lumi", &lumi_, "lumi/i");
  outTree->Branch("evt",  &evt_,  "evt/i");
  outTree->Branch("rho",  &rho_,  "rho/f");

  //outTree->Branch("trueNumInteractions", &trueNumInteractions_, "trueNumInteractions/f");
  //outTree->Branch("puNumInteractions", &puNumInteractions_, "puNumInteractions/i");
  
  // Muons
  /*
  outTree -> Branch("nMuons"         , &nmuons    , "nMuons/I");
  outTree -> Branch("Muon_Pt"        , muon_pt    , "Muon_Pt[nMuons]/F");
  outTree -> Branch("Muon_Eta"       , muon_eta   , "Muon_Eta[nMuons]/F");
  outTree -> Branch("Muon_Phi"       , muon_phi   , "Muon_Phi[nMuons]/F");
  outTree -> Branch("Muon_Mass"      , muon_mass  , "Muon_Mass[nMuons]/F");
  outTree -> Branch("Muon_Charge"    , muon_charge, "Muon_Charge[nMuons]/F");
  outTree -> Branch("Muon_IsLoose"   , muon_IsLoose,  "Muon_IsLoose[nMuons]/O");
  outTree -> Branch("Muon_IsMedium"  , muon_IsMedium, "Muon_IsMediun[nMuons]/O");
  outTree -> Branch("Muon_RelIso"    , muon_relIso,   "Muon_RelIso[nMuons]/F");

  // Electrons
  outTree -> Branch("nelectrons",             &nelectrons,       "nelectrons/I");
  outTree -> Branch("electron_pt",             electron_pt,      "electron_pt[nelectrons]/F");
  outTree -> Branch("electron_eta",            electron_eta,     "electron_eta[nelectrons]/F");
  outTree -> Branch("electron_phi",            electron_phi,     "electron_phi[nelectrons]/F");
  outTree -> Branch("electron_mass",           electron_mass,    "electron_mass[nelectrons]/F");
  outTree -> Branch("electron_reliso",         electron_reliso,  "electron_reliso[nelectrons]/F");
  */
  
  // b-quarks
  outTree -> Branch("nbquarks",    &nbquarks,   "nbquarks/I");
  outTree -> Branch("bquark_pt",   bquark_pt,   "bquark_pt[nbquarks]/F");
  outTree -> Branch("bquark_eta",  bquark_eta,  "bquark_eta[nbquarks]/F");
  outTree -> Branch("bquark_phi",  bquark_phi,  "bquark_phi[nbquarks]/F");
  outTree -> Branch("bquark_mass", bquark_mass, "bquark_mass[nbquarks]/F");

  // AK4 genjets 
  outTree -> Branch("ngenjets", &ngenjets, "ngenjets/I");
  outTree -> Branch("genjet_pt", genjet_pt, "genjet_pt[ngenjets]/F");
  outTree -> Branch("genjet_eta", genjet_eta, "genjet_eta[ngenjets]/F");
  outTree -> Branch("genjet_phi", genjet_phi, "genjet_phi[ngenjets]/F");
  outTree -> Branch("genjet_mass", genjet_mass, "genjet_mass[ngenjets]/F");
  
  // AK8 genfatjets
  outTree -> Branch("ngenfatjets", &ngenfatjets, "ngenfatjets/I");
  outTree -> Branch("genfatjet_pt", genfatjet_pt, "genfatjet_pt[ngenfatjets]/F");
  outTree -> Branch("genfatjet_eta", genfatjet_eta, "genfatjet_eta[ngenfatjets]/F");
  outTree -> Branch("genfatjet_phi", genfatjet_phi, "genfatjet_phi[ngenfatjets]/F");
  outTree -> Branch("genfatjet_mass", genfatjet_mass, "genfatjet_mass[ngenfatjets]/F");
  
  // AK4 Jets
  outTree -> Branch("njets",                  &njets,            "njets/I");
  outTree -> Branch("jet_pt",                  jet_pt,           "jet_pt[njets]/F");
  outTree -> Branch("jet_eta",                 jet_eta,          "jet_eta[njets]/F");
  outTree -> Branch("jet_phi",                 jet_phi,          "jet_phi[njets]/F");
  outTree -> Branch("jet_mass",                jet_mass,         "jet_mass[njets]/F");
  outTree -> Branch("jet_deepjet_probb", jet_deepjet_probb, "jet_deepjet_probb[njets]/F");
  outTree -> Branch("jet_deepjet_probbb", jet_deepjet_probbb, "jet_deepjet_probbb[njets]/F");
  outTree -> Branch("jet_deepjet_problepb", jet_deepjet_problepb, "jet_deepjet_problepb[njets]/F");
  outTree -> Branch("jet_deepjet", jet_deepjet, "jet_deepjet[njets]/F");
  outTree -> Branch("jet_pnet_probb", jet_pnet_probb, "jet_pnet_probb[njets]/F");
  outTree -> Branch("jet_pnet_probbb", jet_pnet_probbb, "jet_pnet_probbb[njets]/F");
  outTree -> Branch("jet_pnet_probc", jet_pnet_probc, "jet_pnet_probc[njets]/F");
  outTree -> Branch("jet_pnet_probcc", jet_pnet_probcc, "jet_pnet_probcc[njets]/F");
  outTree -> Branch("jet_pnet_probuds", jet_pnet_probuds, "jet_pnet_probuds[njets]/F");
  outTree -> Branch("jet_pnet_probg", jet_pnet_probg, "jet_pnet_probg[njets]/F");
  outTree -> Branch("jet_pnet_probundef", jet_pnet_probundef, "jet_pnet_probundef[njets]/F");
  outTree -> Branch("jet_pnet_probpu", jet_pnet_probpu, "jet_pnet_probpu[njets]/F");
  outTree -> Branch("jet_pnet_BvsAll", jet_pnet_BvsAll, "jet_pnet_BvsAll[njets]/F");
  
  // AK8 jets
  outTree -> Branch("nfatjets", &nfatjets, "nfatjets/I");
  outTree -> Branch("fatjet_pt", fatjet_pt, "fatjet_pt[nfatjets]/F");
  outTree -> Branch("fatjet_eta", fatjet_eta, "fatjet_eta[nfatjets]/F");
  outTree -> Branch("fatjet_phi", fatjet_phi, "fatjet_phi[nfatjets]/F");
  outTree -> Branch("fatjet_mass", fatjet_mass, "fatjet_mass[nfatjets]/F");
  outTree -> Branch("fatjet_massSD", fatjet_massSD, "fatjet_massSD[nfatjets]/F");
  outTree -> Branch("fatjet_PNetXbb", fatjet_PNetXbb, "fatjet_PNetXbb[nfatjets]/F");
  return outTree;
}

// ------------ method called once each job just before starting event loop  ------------
void HLTForDiHiggsAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void HLTForDiHiggsAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HLTForDiHiggsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTForDiHiggsAnalyzer);

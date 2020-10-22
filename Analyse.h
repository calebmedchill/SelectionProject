//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 13 15:31:45 2020 by ROOT version 6.23/01
// from TTree EventTree/Event data
// found on file: WGToLNuGPtG130TuneCP513TeVamcatnloFXFXpythia8_000.root
//////////////////////////////////////////////////////////

#ifndef Analyse_h
#define Analyse_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class Analyse {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   UShort_t        lumis;
   UChar_t         nVtx;
   UChar_t         nGoodVtx;
   Bool_t          isPVGood;
   Float_t         vtx;
   Float_t         vty;
   Float_t         vtz;
   Float_t         rho;
   Float_t         rhoCentral;
   ULong64_t       HLTEleMuX;
   ULong64_t       HLTPho;
   ULong64_t       HLTPhoRejectedByPS;
   ULong64_t       HLTJet;
   ULong64_t       HLTEleMuXIsPrescaled;
   ULong64_t       HLTPhoIsPrescaled;
   ULong64_t       HLTJetIsPrescaled;
   Float_t         ecalPrefireW;
   Float_t         ecalPrefireWup;
   Float_t         ecalPrefireWdn;
   UShort_t        beamHaloSummary;
   UChar_t         ntrgObjPho;
   vector<unsigned char> *trgObjPhoBits;
   vector<float>   *trgObjPhoPt;
   vector<float>   *trgObjPhoEta;
   vector<float>   *trgObjPhoPhi;
   vector<float>   *pdf;
   Float_t         genWeight;
   Float_t         genHT;
   Float_t         genPho1;
   Float_t         genPho2;
   UChar_t         nPUInfo;
   UChar_t         nPU;
   UChar_t         puBX;
   UChar_t         puTrue;
   UShort_t        nMC;
   vector<int>     *mcPID;
   vector<float>   *mcVtx;
   vector<float>   *mcVty;
   vector<float>   *mcVtz;
   vector<float>   *mcPt;
   vector<float>   *mcMass;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<int>     *mcGMomPID;
   vector<int>     *mcMomPID;
   vector<float>   *mcMomPt;
   vector<float>   *mcMomMass;
   vector<float>   *mcMomEta;
   vector<float>   *mcMomPhi;
   vector<short>   *mcIndex;
   vector<unsigned short> *mcStatusFlag;
   vector<unsigned char> *mcParentage;
   vector<short>   *mcStatus;
   vector<char>    *mcPromptStatusType;
   Char_t          mcHasDirectPromptPho;
   vector<float>   *mcCalIsoDR03;
   vector<float>   *mcTrkIsoDR03;
   vector<float>   *mcCalIsoDR04;
   vector<float>   *mcTrkIsoDR04;
   UShort_t        nGenAK4Jets;
   vector<float>   *GenAK4JetEn;
   vector<float>   *GenAK4JetPt;
   vector<float>   *GenAK4JetEta;
   vector<float>   *GenAK4JetDetEta;
   vector<float>   *GenAK4JetPhi;
   vector<float>   *GenAK4JetEMenergy;
   vector<float>   *GenAK4JetHADenergy;
   Float_t         genMET;
   Float_t         genMETPhi;
   UShort_t        metFilters;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMET_T1JERUp;
   Float_t         pfMET_T1JERDo;
   Float_t         pfMET_T1JESUp;
   Float_t         pfMET_T1JESDo;
   Float_t         pfMET_T1UESUp;
   Float_t         pfMET_T1UESDo;
   Float_t         pfMETPhi_T1JESUp;
   Float_t         pfMETPhi_T1JESDo;
   Float_t         pfMETPhi_T1UESUp;
   Float_t         pfMETPhi_T1UESDo;
   Float_t         pfMET_metSig;
   Float_t         pfMET_EtSig;
   UShort_t        nPho;
   vector<float>   *phoE;
   vector<float>   *phoSigmaE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoCalibE;
   vector<float>   *phoSigmaCalibE;
   vector<float>   *phoCalibEt;
   vector<short>   *phoSCindex;
   vector<float>   *phoESEnP1;
   vector<float>   *phoESEnP2;
   vector<unsigned char> *phoFiducialRegion;
   vector<unsigned char> *phoQualityBits;
   vector<float>   *phoR9;
   vector<float>   *phoHoverE;
   vector<float>   *phoESEffSigmaRR;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoSigmaIEtaIPhiFull5x5;
   vector<float>   *phoSigmaIPhiIPhiFull5x5;
   vector<float>   *phoE2x2Full5x5;
   vector<float>   *phoE5x5Full5x5;
   vector<float>   *phoMaxEnergyXtal;
   vector<float>   *phoE2ndFull5x5;
   vector<float>   *phoE1x3Full5x5;
   vector<float>   *phoE1x5Full5x5;
   vector<float>   *phoE2x5Full5x5;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoSeedBCE;
   vector<float>   *phoSeedBCEta;
   vector<float>   *phoSeedBCPhi;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoPFChWorstIso;
   vector<float>   *phoPFClusEcalIso;
   vector<float>   *phoPFClusHcalIso;
   vector<unsigned char> *nPhoTrkSolidConeDR03;
   vector<unsigned char> *nPhoTrkHollowConeDR03;
   vector<float>   *phoTrkSumPtSolidConeDR03;
   vector<float>   *phoTrkSumPtHollowConeDR03;
   vector<unsigned char> *nPhoTrkSolidConeDR04;
   vector<unsigned char> *nPhoTrkHollowConeDR04;
   vector<float>   *phoTrkSumPtSolidConeDR04;
   vector<float>   *phoTrkSumPtHollowConeDR04;
   vector<float>   *phoECALIso;
   vector<float>   *phoHCALIso;
   vector<float>   *phoIDMVA;
   vector<ULong64_t> *phoFiredSingleTrgs;
   vector<ULong64_t> *phoFiredDoubleTrgs;
   vector<ULong64_t> *phoFiredTripleTrgs;
   vector<ULong64_t> *phoFiredL1Trgs;
   vector<float>   *phoSeedTime;
   vector<float>   *phoSeedEnergy;
   vector<float>   *phoMIPChi2;
   vector<float>   *phoMIPTotEnergy;
   vector<float>   *phoMIPSlope;
   vector<float>   *phoMIPIntercept;
   vector<short>   *phoMIPNhitCone;
   vector<unsigned char> *phoIDbit;
   vector<float>   *phoScale_stat_up;
   vector<float>   *phoScale_stat_dn;
   vector<float>   *phoScale_syst_up;
   vector<float>   *phoScale_syst_dn;
   vector<float>   *phoScale_gain_up;
   vector<float>   *phoScale_gain_dn;
   vector<float>   *phoResol_rho_up;
   vector<float>   *phoResol_rho_dn;
   vector<float>   *phoResol_phi_up;
   vector<float>   *phoResol_phi_dn;
   vector<int>     *phoNConvLegs;
   vector<float>   *phoZVtxWithConv;
   vector<short>   *pho_gen_index;
   vector<short>   *phoDirectEcalSCindex;
   UShort_t        necalSC;
   vector<float>   *ecalSCeta;
   vector<float>   *ecalSCphi;
   vector<float>   *ecalSCEn;
   vector<float>   *ecalSCRawEn;
   vector<float>   *ecalSCetaWidth;
   vector<float>   *ecalSCphiWidth;
   vector<float>   *ecalSC_LICTD;
   vector<unsigned char> *ecalSC_nL1Spike;
   vector<unsigned char> *ecalSC_nDiweird;
   vector<unsigned char> *ecalSC_nWeird;
   vector<unsigned char> *ecalSC_nSaturated;
   vector<unsigned char> *ecalSC_nOutOfTime;
   vector<unsigned char> *ecalSC_nXtals;
   vector<float>   *ecalSC_maxEnXtalTime;
   vector<float>   *ecalSC_maxEnXtalSwissCross;
   vector<unsigned char> *ecalSC_maxEnXtalBits;
   vector<char>    *ecalSCseedIx;
   vector<char>    *ecalSCseedIy;
   vector<char>    *ecalSCseedIz;
   UShort_t        nootPho;
   vector<float>   *ootPho_E;
   vector<float>   *ootPhoSigmaE;
   vector<float>   *ootPho_Et;
   vector<float>   *ootPhoCalibE;
   vector<float>   *ootPhoSigmaCalibE;
   vector<float>   *ootPhoCalibEt;
   vector<float>   *ootPho_Eta;
   vector<float>   *ootPho_Phi;
   vector<short>   *ootPho_SCindex;
   vector<float>   *ootPhoESEnP1;
   vector<float>   *ootPhoESEnP2;
   vector<unsigned char> *ootPho_FiducialRegion;
   vector<unsigned char> *ootPho_QualityBits;
   vector<float>   *ootPho_R9;
   vector<float>   *ootPho_HoverE;
   vector<float>   *ootPho_ESEffSigmaRR;
   vector<float>   *ootPho_SigmaIEtaIEtaFull5x5;
   vector<float>   *ootPho_SigmaIEtaIPhiFull5x5;
   vector<float>   *ootPho_SigmaIPhiIPhiFull5x5;
   vector<float>   *ootPhoE2x2Full5x5;
   vector<float>   *ootPhoE5x5Full5x5;
   vector<float>   *ootPho_R9Full5x5;
   vector<float>   *ootPhoMaxEnergyXtal;
   vector<float>   *ootPhoE2ndFull5x5;
   vector<float>   *ootPhoE1x3Full5x5;
   vector<float>   *ootPhoE1x5Full5x5;
   vector<float>   *ootPhoE2x5Full5x5;
   vector<float>   *ootPhoPFClusEcalIso;
   vector<float>   *ootPhoPFClusHcalIso;
   vector<unsigned char> *nootPhoTrkHollowConeDR03;
   vector<unsigned char> *nootPhoTrkSolidConeDR03;
   vector<float>   *ootPhoTrkSumPtSolidConeDR03;
   vector<float>   *ootPhoTrkSumPtHollowConeDR03;
   vector<unsigned char> *nootPhoTrkHollowConeDR04;
   vector<unsigned char> *nootPhoTrkSolidConeDR04;
   vector<float>   *ootPhoTrkSumPtSolidConeDR04;
   vector<float>   *ootPhoTrkSumPtHollowConeDR04;
   vector<float>   *ootPhoECALIso;
   vector<float>   *ootPhoHCALIso;
   vector<float>   *ootPhoSeedBCE;
   vector<float>   *ootPhoSeedBCEta;
   vector<float>   *ootPhoSeedBCPhi;
   vector<ULong64_t> *ootPho_FiredSingleTrgs;
   vector<ULong64_t> *ootPho_FiredDoubleTrgs;
   vector<ULong64_t> *ootPho_FiredTripleTrgs;
   vector<ULong64_t> *ootPho_FiredL1Trgs;
   vector<float>   *ootPho_SeedTime;
   vector<float>   *ootPho_SeedEnergy;
   vector<float>   *ootPho_MIPChi2;
   vector<float>   *ootPho_MIPTotEnergy;
   vector<float>   *ootPho_MIPSlope;
   vector<float>   *ootPho_MIPIntercept;
   vector<short>   *ootPho_MIPNhitCone;
   vector<unsigned char> *ootPho_IDbit;
   vector<short>   *ootPhoDirectEcalSCindex;
   UShort_t        necalootSC;
   vector<float>   *ecalootSC_eta;
   vector<float>   *ecalootSC_phi;
   vector<float>   *ecalootSC_En;
   vector<float>   *ecalootSC_RawEn;
   vector<float>   *ecalootSC_etaWidth;
   vector<float>   *ecalootSC_phiWidth;
   vector<float>   *ecalootSC_LICTD;
   vector<unsigned char> *ecalootSC_nL1Spike;
   vector<unsigned char> *ecalootSC_nDiweird;
   vector<unsigned char> *ecalootSC_nWeird;
   vector<unsigned char> *ecalootSC_nSaturated;
   vector<unsigned char> *ecalootSC_nOutOfTime;
   vector<unsigned char> *ecalootSC_nXtals;
   vector<float>   *ecalootSC_maxEnXtalTime;
   vector<float>   *ecalootSC_maxEnXtalSwissCross;
   vector<unsigned char> *ecalootSC_maxEnXtalBits;
   vector<char>    *ecalootSC_seedIx;
   vector<char>    *ecalootSC_seedIy;
   vector<char>    *ecalootSC_seedIz;
   UShort_t        nEle;
   vector<char>    *eleCharge;
   vector<float>   *eleEn;
   vector<float>   *eleEcalEn;
   vector<float>   *elePt;
   vector<float>   *elePtError;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9;
   vector<float>   *eleCalibPt;
   vector<float>   *eleCalibEn;
   vector<short>   *eleSCindex;
   vector<float>   *eleHoverE;
   vector<float>   *eleEoverP;
   vector<float>   *eleEoverPout;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleSigmaIPhiIPhiFull5x5;
   vector<unsigned short> *eleQualityBits;
   vector<unsigned char> *eleMissHits;
   vector<float>   *eleESEffSigmaRR;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *elePFClusEcalIso;
   vector<float>   *elePFClusHcalIso;
   vector<float>   *eleIDMVAIso;
   vector<float>   *eleIDMVANoIso;
   vector<float>   *eleR9Full5x5;
   vector<ULong64_t> *eleFiredSingleTrgs;
   vector<ULong64_t> *eleFiredDoubleTrgs;
   vector<ULong64_t> *eleFiredL1Trgs;
   vector<unsigned char> *eleIDbit;
   vector<short>   *eleGenIndex;
   vector<short>   *eleDirectEcalSCindex;
   UShort_t        nMu;
   vector<float>   *muPt;
   vector<float>   *muEn;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<char>    *muCharge;
   vector<short>   *muType;
   vector<int>     *muIDbit;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muSIP;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<short>   *muTrkLayers;
   vector<short>   *muPixelLayers;
   vector<short>   *muPixelHits;
   vector<short>   *muMuonHits;
   vector<short>   *muStations;
   vector<short>   *muMatches;
   vector<char>    *muTrkQuality;
   vector<float>   *muIsoTrk;
   vector<float>   *muPFChIso;
   vector<float>   *muPFPhoIso;
   vector<float>   *muPFNeuIso;
   vector<float>   *muPFPUIso;
   vector<ULong64_t> *muFiredTrgs;
   vector<ULong64_t> *muFiredL1Trgs;
   vector<float>   *muInnervalidFraction;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muBestTrkPtError;
   vector<float>   *muBestTrkPt;
   vector<char>    *muBestTrkType;
   vector<short>   *muGenIndex;
   UShort_t        nAK4CHSJet;
   vector<char>    *AK4CHSJet_Charge;
   vector<float>   *AK4CHSJet_Pt;
   vector<float>   *AK4CHSJet_En;
   vector<float>   *AK4CHSJet_Eta;
   vector<float>   *AK4CHSJet_Phi;
   vector<float>   *AK4CHSJet_RawPt;
   vector<float>   *AK4CHSJet_RawEn;
   vector<float>   *AK4CHSJet_Mt;
   vector<float>   *AK4CHSJet_Area;
   vector<float>   *AK4CHSJet_LeadTrackPt;
   vector<float>   *AK4CHSJet_LeadTrackEta;
   vector<float>   *AK4CHSJet_LeadTrackPhi;
   vector<int>     *AK4CHSJet_LepTrackPID;
   vector<float>   *AK4CHSJet_LepTrackPt;
   vector<float>   *AK4CHSJet_LepTrackEta;
   vector<float>   *AK4CHSJet_LepTrackPhi;
   vector<float>   *AK4CHSJet_ConsituentEtaPhiSpread;
   vector<float>   *AK4CHSJet_constituentPtDistribution;
   vector<float>   *AK4CHSJet_MaxConstituentDistance;
   vector<float>   *AK4CHSJet_CSV2BJetTags;
   vector<float>   *AK4CHSJet_DeepCSVTags_b;
   vector<float>   *AK4CHSJet_DeepCSVTags_bb;
   vector<float>   *AK4CHSJet_DeepCSVTags_c;
   vector<float>   *AK4CHSJet_DeepCSVTags_udsg;
   vector<float>   *AK4CHSJet_DeepFlavTags_bb;
   vector<float>   *AK4CHSJet_DeepFlavTags_bbb;
   vector<float>   *AK4CHSJet_DeepFlavTags_lepb;
   vector<float>   *AK4CHSJet_DeepFlavTags_c;
   vector<float>   *AK4CHSJet_DeepFlavTags_uds;
   vector<float>   *AK4CHSJet_DeepFlavTags_g;
   vector<float>   *AK4CHSJet_CombMVA2Tags;
   vector<int>     *AK4CHSJet_PartonFlavour;
   vector<int>     *AK4CHSJet_HadronFlavour;
   vector<short>   *AK4CHSJet_GenJetIndex;
   vector<short>   *AK4CHSJet_GenPartonIndex;
   vector<char>    *AK4CHSJet_ID;
   vector<float>   *AK4CHSJet_PUID;
   vector<char>    *AK4CHSJet_PUFullID;
   vector<float>   *AK4CHSJet_qgLikelihood;
   vector<float>   *AK4CHSJet_JECUnc;
   vector<float>   *AK4CHSJet_JER_ptRes;
   vector<float>   *AK4CHSJet_JER_phiRes;
   vector<float>   *AK4CHSJet_JER_SF;
   vector<float>   *AK4CHSJet_JER_SFup;
   vector<float>   *AK4CHSJet_JER_SFdown;
   vector<ULong64_t> *AK4CHSJet_FiredTrgs;
   vector<float>   *AK4CHSJet_CHF;
   vector<float>   *AK4CHSJet_NHF;
   vector<float>   *AK4CHSJet_CEF;
   vector<float>   *AK4CHSJet_NEF;
   vector<unsigned short> *AK4CHSJet_NCH;
   vector<unsigned short> *AK4CHSJet_NNP;
   vector<float>   *AK4CHSJet_MUF;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_isPVGood;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vty;   //!
   TBranch        *b_vtz;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_HLTEleMuX;   //!
   TBranch        *b_HLTPho;   //!
   TBranch        *b_HLTPhoRejectedByPS;   //!
   TBranch        *b_HLTJet;   //!
   TBranch        *b_HLTEleMuXIsPrescaled;   //!
   TBranch        *b_HLTPhoIsPrescaled;   //!
   TBranch        *b_HLTJetIsPrescaled;   //!
   TBranch        *b_ecalPrefireW;   //!
   TBranch        *b_ecalPrefireWup;   //!
   TBranch        *b_ecalPrefireWdn;   //!
   TBranch        *b_beamHaloSummary;   //!
   TBranch        *b_ntrgObjPho;   //!
   TBranch        *b_trgObjPhoBits;   //!
   TBranch        *b_trgObjPhoPt;   //!
   TBranch        *b_trgObjPhoEta;   //!
   TBranch        *b_trgObjPhoPhi;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genHT;   //!
   TBranch        *b_genPho1;   //!
   TBranch        *b_genPho2;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcVty;   //!
   TBranch        *b_mcVtz;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcStatusFlag;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcPromptStatusType;   //!
   TBranch        *b_mcHasDirectPromptPho;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_nGenAK4Jets;   //!
   TBranch        *b_GenAK4JetEn;   //!
   TBranch        *b_GenAK4JetPt;   //!
   TBranch        *b_GenAK4JetEta;   //!
   TBranch        *b_GenAK4JetDetEta;   //!
   TBranch        *b_GenAK4JetPhi;   //!
   TBranch        *b_GenAK4JetEMenergy;   //!
   TBranch        *b_GenAK4JetHADenergy;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMET_T1JERUp;   //!
   TBranch        *b_pfMET_T1JERDo;   //!
   TBranch        *b_pfMET_T1JESUp;   //!
   TBranch        *b_pfMET_T1JESDo;   //!
   TBranch        *b_pfMET_T1UESUp;   //!
   TBranch        *b_pfMET_T1UESDo;   //!
   TBranch        *b_pfMETPhi_T1JESUp;   //!
   TBranch        *b_pfMETPhi_T1JESDo;   //!
   TBranch        *b_pfMETPhi_T1UESUp;   //!
   TBranch        *b_pfMETPhi_T1UESDo;   //!
   TBranch        *b_pfMET_metSig;   //!
   TBranch        *b_pfMET_EtSig;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoSigmaE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoCalibE;   //!
   TBranch        *b_phoSigmaCalibE;   //!
   TBranch        *b_phoCalibEt;   //!
   TBranch        *b_phoSCindex;   //!
   TBranch        *b_phoESEnP1;   //!
   TBranch        *b_phoESEnP2;   //!
   TBranch        *b_phoFiducialRegion;   //!
   TBranch        *b_phoQualityBits;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
   TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_phoE2x2Full5x5;   //!
   TBranch        *b_phoE5x5Full5x5;   //!
   TBranch        *b_phoMaxEnergyXtal;   //!
   TBranch        *b_phoE2ndFull5x5;   //!
   TBranch        *b_phoE1x3Full5x5;   //!
   TBranch        *b_phoE1x5Full5x5;   //!
   TBranch        *b_phoE2x5Full5x5;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoSeedBCE;   //!
   TBranch        *b_phoSeedBCEta;   //!
   TBranch        *b_phoSeedBCPhi;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoPFChWorstIso;   //!
   TBranch        *b_phoPFClusEcalIso;   //!
   TBranch        *b_phoPFClusHcalIso;   //!
   TBranch        *b_nPhoTrkSolidConeDR03;   //!
   TBranch        *b_nPhoTrkHollowConeDR03;   //!
   TBranch        *b_phoTrkSumPtSolidConeDR03;   //!
   TBranch        *b_phoTrkSumPtHollowConeDR03;   //!
   TBranch        *b_nPhoTrkSolidConeDR04;   //!
   TBranch        *b_nPhoTrkHollowConeDR04;   //!
   TBranch        *b_phoTrkSumPtSolidConeDR04;   //!
   TBranch        *b_phoTrkSumPtHollowConeDR04;   //!
   TBranch        *b_phoECALIso;   //!
   TBranch        *b_phoHCALIso;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoFiredSingleTrgs;   //!
   TBranch        *b_phoFiredDoubleTrgs;   //!
   TBranch        *b_phoFiredTripleTrgs;   //!
   TBranch        *b_phoFiredL1Trgs;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoMIPChi2;   //!
   TBranch        *b_phoMIPTotEnergy;   //!
   TBranch        *b_phoMIPSlope;   //!
   TBranch        *b_phoMIPIntercept;   //!
   TBranch        *b_phoMIPNhitCone;   //!
   TBranch        *b_phoIDbit;   //!
   TBranch        *b_phoScale_stat_up;   //!
   TBranch        *b_phoScale_stat_dn;   //!
   TBranch        *b_phoScale_syst_up;   //!
   TBranch        *b_phoScale_syst_dn;   //!
   TBranch        *b_phoScale_gain_up;   //!
   TBranch        *b_phoScale_gain_dn;   //!
   TBranch        *b_phoResol_rho_up;   //!
   TBranch        *b_phoResol_rho_dn;   //!
   TBranch        *b_phoResol_phi_up;   //!
   TBranch        *b_phoResol_phi_dn;   //!
   TBranch        *b_phoNConvLegs;   //!
   TBranch        *b_phoZVtxWithConv;   //!
   TBranch        *b_pho_gen_index;   //!
   TBranch        *b_phoDirectEcalSCindex;   //!
   TBranch        *b_necalSC;   //!
   TBranch        *b_ecalSCeta;   //!
   TBranch        *b_ecalSCphi;   //!
   TBranch        *b_ecalSCEn;   //!
   TBranch        *b_ecalSCRawEn;   //!
   TBranch        *b_ecalSCetaWidth;   //!
   TBranch        *b_ecalSCphiWidth;   //!
   TBranch        *b_ecalSC_LICTD;   //!
   TBranch        *b_ecalSC_nL1Spike;   //!
   TBranch        *b_ecalSC_nDiweird;   //!
   TBranch        *b_ecalSC_nWeird;   //!
   TBranch        *b_ecalSC_nSaturated;   //!
   TBranch        *b_ecalSC_nOutOfTime;   //!
   TBranch        *b_ecalSC_nXtals;   //!
   TBranch        *b_ecalSC_maxEnXtalTime;   //!
   TBranch        *b_ecalSC_maxEnXtalSwissCross;   //!
   TBranch        *b_ecalSC_maxEnXtalBits;   //!
   TBranch        *b_ecalSCseedIx;   //!
   TBranch        *b_ecalSCseedIy;   //!
   TBranch        *b_ecalSCseedIz;   //!
   TBranch        *b_nootPho;   //!
   TBranch        *b_ootPho_E;   //!
   TBranch        *b_ootPhoSigmaE;   //!
   TBranch        *b_ootPho_Et;   //!
   TBranch        *b_ootPhoCalibE;   //!
   TBranch        *b_ootPhoSigmaCalibE;   //!
   TBranch        *b_ootPhoCalibEt;   //!
   TBranch        *b_ootPho_Eta;   //!
   TBranch        *b_ootPho_Phi;   //!
   TBranch        *b_ootPho_SCindex;   //!
   TBranch        *b_ootPhoESEnP1;   //!
   TBranch        *b_ootPhoESEnP2;   //!
   TBranch        *b_ootPho_FiducialRegion;   //!
   TBranch        *b_ootPho_QualityBits;   //!
   TBranch        *b_ootPho_R9;   //!
   TBranch        *b_ootPho_HoverE;   //!
   TBranch        *b_ootPho_ESEffSigmaRR;   //!
   TBranch        *b_ootPho_SigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_ootPho_SigmaIEtaIPhiFull5x5;   //!
   TBranch        *b_ootPho_SigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_ootPhoE2x2Full5x5;   //!
   TBranch        *b_ootPhoE5x5Full5x5;   //!
   TBranch        *b_ootPho_R9Full5x5;   //!
   TBranch        *b_ootPhoMaxEnergyXtal;   //!
   TBranch        *b_ootPhoE2ndFull5x5;   //!
   TBranch        *b_ootPhoE1x3Full5x5;   //!
   TBranch        *b_ootPhoE1x5Full5x5;   //!
   TBranch        *b_ootPhoE2x5Full5x5;   //!
   TBranch        *b_ootPhoPFClusEcalIso;   //!
   TBranch        *b_ootPhoPFClusHcalIso;   //!
   TBranch        *b_nootPhoTrkHollowConeDR03;   //!
   TBranch        *b_nootPhoTrkSolidConeDR03;   //!
   TBranch        *b_ootPhoTrkSumPtSolidConeDR03;   //!
   TBranch        *b_ootPhoTrkSumPtHollowConeDR03;   //!
   TBranch        *b_nootPhoTrkHollowConeDR04;   //!
   TBranch        *b_nootPhoTrkSolidConeDR04;   //!
   TBranch        *b_ootPhoTrkSumPtSolidConeDR04;   //!
   TBranch        *b_ootPhoTrkSumPtHollowConeDR04;   //!
   TBranch        *b_ootPhoECALIso;   //!
   TBranch        *b_ootPhoHCALIso;   //!
   TBranch        *b_ootPhoSeedBCE;   //!
   TBranch        *b_ootPhoSeedBCEta;   //!
   TBranch        *b_ootPhoSeedBCPhi;   //!
   TBranch        *b_ootPho_FiredSingleTrgs;   //!
   TBranch        *b_ootPho_FiredDoubleTrgs;   //!
   TBranch        *b_ootPho_FiredTripleTrgs;   //!
   TBranch        *b_ootPho_FiredL1Trgs;   //!
   TBranch        *b_ootPho_SeedTime;   //!
   TBranch        *b_ootPho_SeedEnergy;   //!
   TBranch        *b_ootPho_MIPChi2;   //!
   TBranch        *b_ootPho_MIPTotEnergy;   //!
   TBranch        *b_ootPho_MIPSlope;   //!
   TBranch        *b_ootPho_MIPIntercept;   //!
   TBranch        *b_ootPho_MIPNhitCone;   //!
   TBranch        *b_ootPho_IDbit;   //!
   TBranch        *b_ootPhoDirectEcalSCindex;   //!
   TBranch        *b_necalootSC;   //!
   TBranch        *b_ecalootSC_eta;   //!
   TBranch        *b_ecalootSC_phi;   //!
   TBranch        *b_ecalootSC_En;   //!
   TBranch        *b_ecalootSC_RawEn;   //!
   TBranch        *b_ecalootSC_etaWidth;   //!
   TBranch        *b_ecalootSC_phiWidth;   //!
   TBranch        *b_ecalootSC_LICTD;   //!
   TBranch        *b_ecalootSC_nL1Spike;   //!
   TBranch        *b_ecalootSC_nDiweird;   //!
   TBranch        *b_ecalootSC_nWeird;   //!
   TBranch        *b_ecalootSC_nSaturated;   //!
   TBranch        *b_ecalootSC_nOutOfTime;   //!
   TBranch        *b_ecalootSC_nXtals;   //!
   TBranch        *b_ecalootSC_maxEnXtalTime;   //!
   TBranch        *b_ecalootSC_maxEnXtalSwissCross;   //!
   TBranch        *b_ecalootSC_maxEnXtalBits;   //!
   TBranch        *b_ecalootSC_seedIx;   //!
   TBranch        *b_ecalootSC_seedIy;   //!
   TBranch        *b_ecalootSC_seedIz;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleEcalEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_elePtError;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleCalibPt;   //!
   TBranch        *b_eleCalibEn;   //!
   TBranch        *b_eleSCindex;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPout;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_eleQualityBits;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_elePFClusEcalIso;   //!
   TBranch        *b_elePFClusHcalIso;   //!
   TBranch        *b_eleIDMVAIso;   //!
   TBranch        *b_eleIDMVANoIso;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleFiredSingleTrgs;   //!
   TBranch        *b_eleFiredDoubleTrgs;   //!
   TBranch        *b_eleFiredL1Trgs;   //!
   TBranch        *b_eleIDbit;   //!
   TBranch        *b_eleGenIndex;   //!
   TBranch        *b_eleDirectEcalSCindex;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEn;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIDbit;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muSIP;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muMatches;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
   TBranch        *b_muFiredTrgs;   //!
   TBranch        *b_muFiredL1Trgs;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muBestTrkPtError;   //!
   TBranch        *b_muBestTrkPt;   //!
   TBranch        *b_muBestTrkType;   //!
   TBranch        *b_muGenIndex;   //!
   TBranch        *b_nAK4CHSJet;   //!
   TBranch        *b_AK4CHSJet_Charge;   //!
   TBranch        *b_AK4CHSJet_Pt;   //!
   TBranch        *b_AK4CHSJet_En;   //!
   TBranch        *b_AK4CHSJet_Eta;   //!
   TBranch        *b_AK4CHSJet_Phi;   //!
   TBranch        *b_AK4CHSJet_RawPt;   //!
   TBranch        *b_AK4CHSJet_RawEn;   //!
   TBranch        *b_AK4CHSJet_Mt;   //!
   TBranch        *b_AK4CHSJet_Area;   //!
   TBranch        *b_AK4CHSJet_LeadTrackPt;   //!
   TBranch        *b_AK4CHSJet_LeadTrackEta;   //!
   TBranch        *b_AK4CHSJet_LeadTrackPhi;   //!
   TBranch        *b_AK4CHSJet_LepTrackPID;   //!
   TBranch        *b_AK4CHSJet_LepTrackPt;   //!
   TBranch        *b_AK4CHSJet_LepTrackEta;   //!
   TBranch        *b_AK4CHSJet_LepTrackPhi;   //!
   TBranch        *b_AK4CHSJet_ConsituentEtaPhiSpread;   //!
   TBranch        *b_AK4CHSJet_constituentPtDistribution;   //!
   TBranch        *b_AK4CHSJet_MaxConstituentDistance;   //!
   TBranch        *b_AK4CHSJet_CSV2BJetTags;   //!
   TBranch        *b_AK4CHSJet_DeepCSVTags_b;   //!
   TBranch        *b_AK4CHSJet_DeepCSVTags_bb;   //!
   TBranch        *b_AK4CHSJet_DeepCSVTags_c;   //!
   TBranch        *b_AK4CHSJet_DeepCSVTags_udsg;   //!
   TBranch        *b_AK4CHSJet_DeepFlavTags_bb;   //!
   TBranch        *b_AK4CHSJet_DeepFlavTags_bbb;   //!
   TBranch        *b_AK4CHSJet_DeepFlavTags_lepb;   //!
   TBranch        *b_AK4CHSJet_DeepFlavTags_c;   //!
   TBranch        *b_AK4CHSJet_DeepFlavTags_uds;   //!
   TBranch        *b_AK4CHSJet_DeepFlavTags_g;   //!
   TBranch        *b_AK4CHSJet_CombMVA2Tags;   //!
   TBranch        *b_AK4CHSJet_PartonFlavour;   //!
   TBranch        *b_AK4CHSJet_HadronFlavour;   //!
   TBranch        *b_AK4CHSJet_GenJetIndex;   //!
   TBranch        *b_AK4CHSJet_GenPartonIndex;   //!
   TBranch        *b_AK4CHSJet_ID;   //!
   TBranch        *b_AK4CHSJet_PUID;   //!
   TBranch        *b_AK4CHSJet_PUFullID;   //!
   TBranch        *b_AK4CHSJet_qgLikelihood;   //!
   TBranch        *b_AK4CHSJet_JECUnc;   //!
   TBranch        *b_AK4CHSJet_JER_ptRes;   //!
   TBranch        *b_AK4CHSJet_JER_phiRes;   //!
   TBranch        *b_AK4CHSJet_JER_SF;   //!
   TBranch        *b_AK4CHSJet_JER_SFup;   //!
   TBranch        *b_AK4CHSJet_JER_SFdown;   //!
   TBranch        *b_AK4CHSJet_FiredTrgs;   //!
   TBranch        *b_AK4CHSJet_CHF;   //!
   TBranch        *b_AK4CHSJet_NHF;   //!
   TBranch        *b_AK4CHSJet_CEF;   //!
   TBranch        *b_AK4CHSJet_NEF;   //!
   TBranch        *b_AK4CHSJet_NCH;   //!
   TBranch        *b_AK4CHSJet_NNP;   //!
   TBranch        *b_AK4CHSJet_MUF;   //!

   Analyse(TTree *tree=0);
   virtual ~Analyse();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analyse_cxx
Analyse::Analyse(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("WGToLNuGPtG130TuneCP513TeVamcatnloFXFXpythia8_000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("WGToLNuGPtG130TuneCP513TeVamcatnloFXFXpythia8_000.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("WGToLNuGPtG130TuneCP513TeVamcatnloFXFXpythia8_000.root:/ggNtuplizer");
      dir->GetObject("EventTree",tree);

   }
   Init(tree);
}

Analyse::~Analyse()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyse::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyse::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analyse::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trgObjPhoBits = 0;
   trgObjPhoPt = 0;
   trgObjPhoEta = 0;
   trgObjPhoPhi = 0;
   pdf = 0;
   mcPID = 0;
   mcVtx = 0;
   mcVty = 0;
   mcVtz = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomMass = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcIndex = 0;
   mcStatusFlag = 0;
   mcParentage = 0;
   mcStatus = 0;
   mcPromptStatusType = 0;
   mcCalIsoDR03 = 0;
   mcTrkIsoDR03 = 0;
   mcCalIsoDR04 = 0;
   mcTrkIsoDR04 = 0;
   GenAK4JetEn = 0;
   GenAK4JetPt = 0;
   GenAK4JetEta = 0;
   GenAK4JetDetEta = 0;
   GenAK4JetPhi = 0;
   GenAK4JetEMenergy = 0;
   GenAK4JetHADenergy = 0;
   phoE = 0;
   phoSigmaE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoCalibE = 0;
   phoSigmaCalibE = 0;
   phoCalibEt = 0;
   phoSCindex = 0;
   phoESEnP1 = 0;
   phoESEnP2 = 0;
   phoFiducialRegion = 0;
   phoQualityBits = 0;
   phoR9 = 0;
   phoHoverE = 0;
   phoESEffSigmaRR = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIPhiFull5x5 = 0;
   phoSigmaIPhiIPhiFull5x5 = 0;
   phoE2x2Full5x5 = 0;
   phoE5x5Full5x5 = 0;
   phoMaxEnergyXtal = 0;
   phoE2ndFull5x5 = 0;
   phoE1x3Full5x5 = 0;
   phoE1x5Full5x5 = 0;
   phoE2x5Full5x5 = 0;
   phoR9Full5x5 = 0;
   phoSeedBCE = 0;
   phoSeedBCEta = 0;
   phoSeedBCPhi = 0;
   phoPFChIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoPFChWorstIso = 0;
   phoPFClusEcalIso = 0;
   phoPFClusHcalIso = 0;
   nPhoTrkSolidConeDR03 = 0;
   nPhoTrkHollowConeDR03 = 0;
   phoTrkSumPtSolidConeDR03 = 0;
   phoTrkSumPtHollowConeDR03 = 0;
   nPhoTrkSolidConeDR04 = 0;
   nPhoTrkHollowConeDR04 = 0;
   nPhoTrkSolidConeDR04 = 0;
   phoTrkSumPtSolidConeDR04 = 0;
   phoTrkSumPtHollowConeDR04 = 0;
   phoECALIso = 0;
   phoHCALIso = 0;
   phoIDMVA = 0;
   phoFiredSingleTrgs = 0;
   phoFiredDoubleTrgs = 0;
   phoFiredTripleTrgs = 0;
   phoFiredL1Trgs = 0;
   phoSeedTime = 0;
   phoSeedEnergy = 0;
   phoMIPChi2 = 0;
   phoMIPTotEnergy = 0;
   phoMIPSlope = 0;
   phoMIPIntercept = 0;
   phoMIPNhitCone = 0;
   phoIDbit = 0;
   phoScale_stat_up = 0;
   phoScale_stat_dn = 0;
   phoScale_syst_up = 0;
   phoScale_syst_dn = 0;
   phoScale_gain_up = 0;
   phoScale_gain_dn = 0;
   phoResol_rho_up = 0;
   phoResol_rho_dn = 0;
   phoResol_phi_up = 0;
   phoResol_phi_dn = 0;
   phoNConvLegs = 0;
   phoZVtxWithConv = 0;
   pho_gen_index = 0;
   phoDirectEcalSCindex = 0;
   ecalSCeta = 0;
   ecalSCphi = 0;
   ecalSCEn = 0;
   ecalSCRawEn = 0;
   ecalSCetaWidth = 0;
   ecalSCphiWidth = 0;
   ecalSC_LICTD = 0;
   ecalSC_nL1Spike = 0;
   ecalSC_nDiweird = 0;
   ecalSC_nWeird = 0;
   ecalSC_nSaturated = 0;
   ecalSC_nOutOfTime = 0;
   ecalSC_nXtals = 0;
   ecalSC_maxEnXtalTime = 0;
   ecalSC_maxEnXtalSwissCross = 0;
   ecalSC_maxEnXtalBits = 0;
   ecalSCseedIx = 0;
   ecalSCseedIy = 0;
   ecalSCseedIz = 0;
   ootPho_E = 0;
   ootPhoSigmaE = 0;
   ootPho_Et = 0;
   ootPhoCalibE = 0;
   ootPhoSigmaCalibE = 0;
   ootPhoCalibEt = 0;
   ootPho_Eta = 0;
   ootPho_Phi = 0;
   ootPho_SCindex = 0;
   ootPhoESEnP1 = 0;
   ootPhoESEnP2 = 0;
   ootPho_FiducialRegion = 0;
   ootPho_QualityBits = 0;
   ootPho_R9 = 0;
   ootPho_HoverE = 0;
   ootPho_ESEffSigmaRR = 0;
   ootPho_SigmaIEtaIEtaFull5x5 = 0;
   ootPho_SigmaIEtaIPhiFull5x5 = 0;
   ootPho_SigmaIPhiIPhiFull5x5 = 0;
   ootPhoE2x2Full5x5 = 0;
   ootPhoE5x5Full5x5 = 0;
   ootPho_R9Full5x5 = 0;
   ootPhoMaxEnergyXtal = 0;
   ootPhoE2ndFull5x5 = 0;
   ootPhoE1x3Full5x5 = 0;
   ootPhoE1x5Full5x5 = 0;
   ootPhoE2x5Full5x5 = 0;
   ootPhoPFClusEcalIso = 0;
   ootPhoPFClusHcalIso = 0;
   nootPhoTrkHollowConeDR03 = 0;
   nootPhoTrkSolidConeDR03 = 0;
   ootPhoTrkSumPtSolidConeDR03 = 0;
   ootPhoTrkSumPtHollowConeDR03 = 0;
   nootPhoTrkHollowConeDR04 = 0;
   nootPhoTrkSolidConeDR04 = 0;
   ootPhoTrkSumPtSolidConeDR04 = 0;
   ootPhoTrkSumPtHollowConeDR04 = 0;
   ootPhoECALIso = 0;
   ootPhoHCALIso = 0;
   ootPhoSeedBCE = 0;
   ootPhoSeedBCEta = 0;
   ootPhoSeedBCPhi = 0;
   ootPho_FiredSingleTrgs = 0;
   ootPho_FiredDoubleTrgs = 0;
   ootPho_FiredTripleTrgs = 0;
   ootPho_FiredL1Trgs = 0;
   ootPho_SeedTime = 0;
   ootPho_SeedEnergy = 0;
   ootPho_MIPChi2 = 0;
   ootPho_MIPTotEnergy = 0;
   ootPho_MIPSlope = 0;
   ootPho_MIPIntercept = 0;
   ootPho_MIPNhitCone = 0;
   ootPho_IDbit = 0;
   ootPhoDirectEcalSCindex = 0;
   ecalootSC_eta = 0;
   ecalootSC_phi = 0;
   ecalootSC_En = 0;
   ecalootSC_RawEn = 0;
   ecalootSC_etaWidth = 0;
   ecalootSC_phiWidth = 0;
   ecalootSC_LICTD = 0;
   ecalootSC_nL1Spike = 0;
   ecalootSC_nDiweird = 0;
   ecalootSC_nWeird = 0;
   ecalootSC_nSaturated = 0;
   ecalootSC_nOutOfTime = 0;
   ecalootSC_nXtals = 0;
   ecalootSC_maxEnXtalTime = 0;
   ecalootSC_maxEnXtalSwissCross = 0;
   ecalootSC_maxEnXtalBits = 0;
   ecalootSC_seedIx = 0;
   ecalootSC_seedIy = 0;
   ecalootSC_seedIz = 0;
   eleCharge = 0;
   eleEn = 0;
   eleEcalEn = 0;
   elePt = 0;
   elePtError = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleCalibPt = 0;
   eleCalibEn = 0;
   eleSCindex = 0;
   eleHoverE = 0;
   eleEoverP = 0;
   eleEoverPout = 0;
   eleEoverPInv = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   eleQualityBits = 0;
   eleMissHits = 0;
   eleESEffSigmaRR = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   elePFClusEcalIso = 0;
   elePFClusHcalIso = 0;
   eleIDMVAIso = 0;
   eleIDMVANoIso = 0;
   eleR9Full5x5 = 0;
   eleFiredSingleTrgs = 0;
   eleFiredDoubleTrgs = 0;
   eleFiredL1Trgs = 0;
   eleIDbit = 0;
   eleGenIndex = 0;
   eleDirectEcalSCindex = 0;
   muPt = 0;
   muEn = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIDbit = 0;
   muD0 = 0;
   muDz = 0;
   muSIP = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muStations = 0;
   muMatches = 0;
   muTrkQuality = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   muFiredTrgs = 0;
   muFiredL1Trgs = 0;
   muInnervalidFraction = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   muBestTrkType = 0;
   muGenIndex = 0;
   AK4CHSJet_Charge = 0;
   AK4CHSJet_Pt = 0;
   AK4CHSJet_En = 0;
   AK4CHSJet_Eta = 0;
   AK4CHSJet_Phi = 0;
   AK4CHSJet_RawPt = 0;
   AK4CHSJet_RawEn = 0;
   AK4CHSJet_Mt = 0;
   AK4CHSJet_Area = 0;
   AK4CHSJet_LeadTrackPt = 0;
   AK4CHSJet_LeadTrackEta = 0;
   AK4CHSJet_LeadTrackPhi = 0;
   AK4CHSJet_LepTrackPID = 0;
   AK4CHSJet_LepTrackPt = 0;
   AK4CHSJet_LepTrackEta = 0;
   AK4CHSJet_LepTrackPhi = 0;
   AK4CHSJet_ConsituentEtaPhiSpread = 0;
   AK4CHSJet_constituentPtDistribution = 0;
   AK4CHSJet_MaxConstituentDistance = 0;
   AK4CHSJet_CSV2BJetTags = 0;
   AK4CHSJet_DeepCSVTags_b = 0;
   AK4CHSJet_DeepCSVTags_bb = 0;
   AK4CHSJet_DeepCSVTags_c = 0;
   AK4CHSJet_DeepCSVTags_udsg = 0;
   AK4CHSJet_DeepFlavTags_bb = 0;
   AK4CHSJet_DeepFlavTags_bbb = 0;
   AK4CHSJet_DeepFlavTags_lepb = 0;
   AK4CHSJet_DeepFlavTags_c = 0;
   AK4CHSJet_DeepFlavTags_uds = 0;
   AK4CHSJet_DeepFlavTags_g = 0;
   AK4CHSJet_CombMVA2Tags = 0;
   AK4CHSJet_PartonFlavour = 0;
   AK4CHSJet_HadronFlavour = 0;
   AK4CHSJet_GenJetIndex = 0;
   AK4CHSJet_GenPartonIndex = 0;
   AK4CHSJet_ID = 0;
   AK4CHSJet_PUID = 0;
   AK4CHSJet_PUFullID = 0;
   AK4CHSJet_qgLikelihood = 0;
   AK4CHSJet_JECUnc = 0;
   AK4CHSJet_JER_ptRes = 0;
   AK4CHSJet_JER_phiRes = 0;
   AK4CHSJet_JER_SF = 0;
   AK4CHSJet_JER_SFup = 0;
   AK4CHSJet_JER_SFdown = 0;
   AK4CHSJet_FiredTrgs = 0;
   AK4CHSJet_CHF = 0;
   AK4CHSJet_NHF = 0;
   AK4CHSJet_CEF = 0;
   AK4CHSJet_NEF = 0;
   AK4CHSJet_NCH = 0;
   AK4CHSJet_NNP = 0;
   AK4CHSJet_MUF = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("isPVGood", &isPVGood, &b_isPVGood);
   fChain->SetBranchAddress("vtx", &vtx, &b_vtx);
   fChain->SetBranchAddress("vty", &vty, &b_vty);
   fChain->SetBranchAddress("vtz", &vtz, &b_vtz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
   fChain->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   fChain->SetBranchAddress("HLTPhoRejectedByPS", &HLTPhoRejectedByPS, &b_HLTPhoRejectedByPS);
   fChain->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
   fChain->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
   fChain->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   fChain->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
   fChain->SetBranchAddress("ecalPrefireW", &ecalPrefireW, &b_ecalPrefireW);
   fChain->SetBranchAddress("ecalPrefireWup", &ecalPrefireWup, &b_ecalPrefireWup);
   fChain->SetBranchAddress("ecalPrefireWdn", &ecalPrefireWdn, &b_ecalPrefireWdn);
   fChain->SetBranchAddress("beamHaloSummary", &beamHaloSummary, &b_beamHaloSummary);
   fChain->SetBranchAddress("ntrgObjPho", &ntrgObjPho, &b_ntrgObjPho);
   fChain->SetBranchAddress("trgObjPhoBits", &trgObjPhoBits, &b_trgObjPhoBits);
   fChain->SetBranchAddress("trgObjPhoPt", &trgObjPhoPt, &b_trgObjPhoPt);
   fChain->SetBranchAddress("trgObjPhoEta", &trgObjPhoEta, &b_trgObjPhoEta);
   fChain->SetBranchAddress("trgObjPhoPhi", &trgObjPhoPhi, &b_trgObjPhoPhi);
   fChain->SetBranchAddress("pdf", &pdf, &b_pdf);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("genPho1", &genPho1, &b_genPho1);
   fChain->SetBranchAddress("genPho2", &genPho2, &b_genPho2);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", &mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcVty", &mcVty, &b_mcVty);
   fChain->SetBranchAddress("mcVtz", &mcVtz, &b_mcVtz);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcPromptStatusType", &mcPromptStatusType, &b_mcPromptStatusType);
   fChain->SetBranchAddress("mcHasDirectPromptPho", &mcHasDirectPromptPho, &b_mcHasDirectPromptPho);
   fChain->SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03, &b_mcCalIsoDR03);
   fChain->SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03, &b_mcTrkIsoDR03);
   fChain->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04, &b_mcCalIsoDR04);
   fChain->SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04, &b_mcTrkIsoDR04);
   fChain->SetBranchAddress("nGenAK4Jets", &nGenAK4Jets, &b_nGenAK4Jets);
   fChain->SetBranchAddress("GenAK4JetEn", &GenAK4JetEn, &b_GenAK4JetEn);
   fChain->SetBranchAddress("GenAK4JetPt", &GenAK4JetPt, &b_GenAK4JetPt);
   fChain->SetBranchAddress("GenAK4JetEta", &GenAK4JetEta, &b_GenAK4JetEta);
   fChain->SetBranchAddress("GenAK4JetDetEta", &GenAK4JetDetEta, &b_GenAK4JetDetEta);
   fChain->SetBranchAddress("GenAK4JetPhi", &GenAK4JetPhi, &b_GenAK4JetPhi);
   fChain->SetBranchAddress("GenAK4JetEMenergy", &GenAK4JetEMenergy, &b_GenAK4JetEMenergy);
   fChain->SetBranchAddress("GenAK4JetHADenergy", &GenAK4JetHADenergy, &b_GenAK4JetHADenergy);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
   fChain->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
   fChain->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
   fChain->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
   fChain->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
   fChain->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
   fChain->SetBranchAddress("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp, &b_pfMETPhi_T1JESUp);
   fChain->SetBranchAddress("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo, &b_pfMETPhi_T1JESDo);
   fChain->SetBranchAddress("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp, &b_pfMETPhi_T1UESUp);
   fChain->SetBranchAddress("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo, &b_pfMETPhi_T1UESDo);
   fChain->SetBranchAddress("pfMET_metSig", &pfMET_metSig, &b_pfMET_metSig);
   fChain->SetBranchAddress("pfMET_EtSig", &pfMET_EtSig, &b_pfMET_EtSig);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoSigmaE", &phoSigmaE, &b_phoSigmaE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoCalibE", &phoCalibE, &b_phoCalibE);
   fChain->SetBranchAddress("phoSigmaCalibE", &phoSigmaCalibE, &b_phoSigmaCalibE);
   fChain->SetBranchAddress("phoCalibEt", &phoCalibEt, &b_phoCalibEt);
   fChain->SetBranchAddress("phoSCindex", &phoSCindex, &b_phoSCindex);
   fChain->SetBranchAddress("phoESEnP1", &phoESEnP1, &b_phoESEnP1);
   fChain->SetBranchAddress("phoESEnP2", &phoESEnP2, &b_phoESEnP2);
   fChain->SetBranchAddress("phoFiducialRegion", &phoFiducialRegion, &b_phoFiducialRegion);
   fChain->SetBranchAddress("phoQualityBits", &phoQualityBits, &b_phoQualityBits);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoE2x2Full5x5", &phoE2x2Full5x5, &b_phoE2x2Full5x5);
   fChain->SetBranchAddress("phoE5x5Full5x5", &phoE5x5Full5x5, &b_phoE5x5Full5x5);
   fChain->SetBranchAddress("phoMaxEnergyXtal", &phoMaxEnergyXtal, &b_phoMaxEnergyXtal);
   fChain->SetBranchAddress("phoE2ndFull5x5", &phoE2ndFull5x5, &b_phoE2ndFull5x5);
   fChain->SetBranchAddress("phoE1x3Full5x5", &phoE1x3Full5x5, &b_phoE1x3Full5x5);
   fChain->SetBranchAddress("phoE1x5Full5x5", &phoE1x5Full5x5, &b_phoE1x5Full5x5);
   fChain->SetBranchAddress("phoE2x5Full5x5", &phoE2x5Full5x5, &b_phoE2x5Full5x5);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoSeedBCE", &phoSeedBCE, &b_phoSeedBCE);
   fChain->SetBranchAddress("phoSeedBCEta", &phoSeedBCEta, &b_phoSeedBCEta);
   fChain->SetBranchAddress("phoSeedBCPhi", &phoSeedBCPhi, &b_phoSeedBCPhi);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   fChain->SetBranchAddress("phoPFClusEcalIso", &phoPFClusEcalIso, &b_phoPFClusEcalIso);
   fChain->SetBranchAddress("phoPFClusHcalIso", &phoPFClusHcalIso, &b_phoPFClusHcalIso);
   fChain->SetBranchAddress("nPhoTrkSolidConeDR03", &nPhoTrkSolidConeDR03, &b_nPhoTrkSolidConeDR03);
   fChain->SetBranchAddress("nPhoTrkHollowConeDR03", &nPhoTrkHollowConeDR03, &b_nPhoTrkHollowConeDR03);
   fChain->SetBranchAddress("phoTrkSumPtSolidConeDR03", &phoTrkSumPtSolidConeDR03, &b_phoTrkSumPtSolidConeDR03);
   fChain->SetBranchAddress("phoTrkSumPtHollowConeDR03", &phoTrkSumPtHollowConeDR03, &b_phoTrkSumPtHollowConeDR03);
   fChain->SetBranchAddress("nPhoTrkSolidConeDR04", &nPhoTrkSolidConeDR04, &b_nPhoTrkSolidConeDR04);
   fChain->SetBranchAddress("nPhoTrkHollowConeDR04", &nPhoTrkHollowConeDR04, &b_nPhoTrkHollowConeDR04);
//    fChain->SetBranchAddress("nPhoTrkSolidConeDR04", &nPhoTrkSolidConeDR04, &b_nPhoTrkSolidConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04, &b_phoTrkSumPtSolidConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04, &b_phoTrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("phoECALIso", &phoECALIso, &b_phoECALIso);
   fChain->SetBranchAddress("phoHCALIso", &phoHCALIso, &b_phoHCALIso);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   fChain->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   fChain->SetBranchAddress("phoFiredTripleTrgs", &phoFiredTripleTrgs, &b_phoFiredTripleTrgs);
   fChain->SetBranchAddress("phoFiredL1Trgs", &phoFiredL1Trgs, &b_phoFiredL1Trgs);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedEnergy", &phoSeedEnergy, &b_phoSeedEnergy);
   fChain->SetBranchAddress("phoMIPChi2", &phoMIPChi2, &b_phoMIPChi2);
   fChain->SetBranchAddress("phoMIPTotEnergy", &phoMIPTotEnergy, &b_phoMIPTotEnergy);
   fChain->SetBranchAddress("phoMIPSlope", &phoMIPSlope, &b_phoMIPSlope);
   fChain->SetBranchAddress("phoMIPIntercept", &phoMIPIntercept, &b_phoMIPIntercept);
   fChain->SetBranchAddress("phoMIPNhitCone", &phoMIPNhitCone, &b_phoMIPNhitCone);
   fChain->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   fChain->SetBranchAddress("phoScale_stat_up", &phoScale_stat_up, &b_phoScale_stat_up);
   fChain->SetBranchAddress("phoScale_stat_dn", &phoScale_stat_dn, &b_phoScale_stat_dn);
   fChain->SetBranchAddress("phoScale_syst_up", &phoScale_syst_up, &b_phoScale_syst_up);
   fChain->SetBranchAddress("phoScale_syst_dn", &phoScale_syst_dn, &b_phoScale_syst_dn);
   fChain->SetBranchAddress("phoScale_gain_up", &phoScale_gain_up, &b_phoScale_gain_up);
   fChain->SetBranchAddress("phoScale_gain_dn", &phoScale_gain_dn, &b_phoScale_gain_dn);
   fChain->SetBranchAddress("phoResol_rho_up", &phoResol_rho_up, &b_phoResol_rho_up);
   fChain->SetBranchAddress("phoResol_rho_dn", &phoResol_rho_dn, &b_phoResol_rho_dn);
   fChain->SetBranchAddress("phoResol_phi_up", &phoResol_phi_up, &b_phoResol_phi_up);
   fChain->SetBranchAddress("phoResol_phi_dn", &phoResol_phi_dn, &b_phoResol_phi_dn);
   fChain->SetBranchAddress("phoNConvLegs", &phoNConvLegs, &b_phoNConvLegs);
   fChain->SetBranchAddress("phoZVtxWithConv", &phoZVtxWithConv, &b_phoZVtxWithConv);
   fChain->SetBranchAddress("pho_gen_index", &pho_gen_index, &b_pho_gen_index);
   fChain->SetBranchAddress("phoDirectEcalSCindex", &phoDirectEcalSCindex, &b_phoDirectEcalSCindex);
   fChain->SetBranchAddress("necalSC", &necalSC, &b_necalSC);
   fChain->SetBranchAddress("ecalSCeta", &ecalSCeta, &b_ecalSCeta);
   fChain->SetBranchAddress("ecalSCphi", &ecalSCphi, &b_ecalSCphi);
   fChain->SetBranchAddress("ecalSCEn", &ecalSCEn, &b_ecalSCEn);
   fChain->SetBranchAddress("ecalSCRawEn", &ecalSCRawEn, &b_ecalSCRawEn);
   fChain->SetBranchAddress("ecalSCetaWidth", &ecalSCetaWidth, &b_ecalSCetaWidth);
   fChain->SetBranchAddress("ecalSCphiWidth", &ecalSCphiWidth, &b_ecalSCphiWidth);
   fChain->SetBranchAddress("ecalSC_LICTD", &ecalSC_LICTD, &b_ecalSC_LICTD);
   fChain->SetBranchAddress("ecalSC_nL1Spike", &ecalSC_nL1Spike, &b_ecalSC_nL1Spike);
   fChain->SetBranchAddress("ecalSC_nDiweird", &ecalSC_nDiweird, &b_ecalSC_nDiweird);
   fChain->SetBranchAddress("ecalSC_nWeird", &ecalSC_nWeird, &b_ecalSC_nWeird);
   fChain->SetBranchAddress("ecalSC_nSaturated", &ecalSC_nSaturated, &b_ecalSC_nSaturated);
   fChain->SetBranchAddress("ecalSC_nOutOfTime", &ecalSC_nOutOfTime, &b_ecalSC_nOutOfTime);
   fChain->SetBranchAddress("ecalSC_nXtals", &ecalSC_nXtals, &b_ecalSC_nXtals);
   fChain->SetBranchAddress("ecalSC_maxEnXtalTime", &ecalSC_maxEnXtalTime, &b_ecalSC_maxEnXtalTime);
   fChain->SetBranchAddress("ecalSC_maxEnXtalSwissCross", &ecalSC_maxEnXtalSwissCross, &b_ecalSC_maxEnXtalSwissCross);
   fChain->SetBranchAddress("ecalSC_maxEnXtalBits", &ecalSC_maxEnXtalBits, &b_ecalSC_maxEnXtalBits);
   fChain->SetBranchAddress("ecalSCseedIx", &ecalSCseedIx, &b_ecalSCseedIx);
   fChain->SetBranchAddress("ecalSCseedIy", &ecalSCseedIy, &b_ecalSCseedIy);
   fChain->SetBranchAddress("ecalSCseedIz", &ecalSCseedIz, &b_ecalSCseedIz);
   fChain->SetBranchAddress("nootPho", &nootPho, &b_nootPho);
   fChain->SetBranchAddress("ootPho_E", &ootPho_E, &b_ootPho_E);
   fChain->SetBranchAddress("ootPhoSigmaE", &ootPhoSigmaE, &b_ootPhoSigmaE);
   fChain->SetBranchAddress("ootPho_Et", &ootPho_Et, &b_ootPho_Et);
   fChain->SetBranchAddress("ootPhoCalibE", &ootPhoCalibE, &b_ootPhoCalibE);
   fChain->SetBranchAddress("ootPhoSigmaCalibE", &ootPhoSigmaCalibE, &b_ootPhoSigmaCalibE);
   fChain->SetBranchAddress("ootPhoCalibEt", &ootPhoCalibEt, &b_ootPhoCalibEt);
   fChain->SetBranchAddress("ootPho_Eta", &ootPho_Eta, &b_ootPho_Eta);
   fChain->SetBranchAddress("ootPho_Phi", &ootPho_Phi, &b_ootPho_Phi);
   fChain->SetBranchAddress("ootPho_SCindex", &ootPho_SCindex, &b_ootPho_SCindex);
   fChain->SetBranchAddress("ootPhoESEnP1", &ootPhoESEnP1, &b_ootPhoESEnP1);
   fChain->SetBranchAddress("ootPhoESEnP2", &ootPhoESEnP2, &b_ootPhoESEnP2);
   fChain->SetBranchAddress("ootPho_FiducialRegion", &ootPho_FiducialRegion, &b_ootPho_FiducialRegion);
   fChain->SetBranchAddress("ootPho_QualityBits", &ootPho_QualityBits, &b_ootPho_QualityBits);
   fChain->SetBranchAddress("ootPho_R9", &ootPho_R9, &b_ootPho_R9);
   fChain->SetBranchAddress("ootPho_HoverE", &ootPho_HoverE, &b_ootPho_HoverE);
   fChain->SetBranchAddress("ootPho_ESEffSigmaRR", &ootPho_ESEffSigmaRR, &b_ootPho_ESEffSigmaRR);
   fChain->SetBranchAddress("ootPho_SigmaIEtaIEtaFull5x5", &ootPho_SigmaIEtaIEtaFull5x5, &b_ootPho_SigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("ootPho_SigmaIEtaIPhiFull5x5", &ootPho_SigmaIEtaIPhiFull5x5, &b_ootPho_SigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("ootPho_SigmaIPhiIPhiFull5x5", &ootPho_SigmaIPhiIPhiFull5x5, &b_ootPho_SigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("ootPhoE2x2Full5x5", &ootPhoE2x2Full5x5, &b_ootPhoE2x2Full5x5);
   fChain->SetBranchAddress("ootPhoE5x5Full5x5", &ootPhoE5x5Full5x5, &b_ootPhoE5x5Full5x5);
   fChain->SetBranchAddress("ootPho_R9Full5x5", &ootPho_R9Full5x5, &b_ootPho_R9Full5x5);
   fChain->SetBranchAddress("ootPhoMaxEnergyXtal", &ootPhoMaxEnergyXtal, &b_ootPhoMaxEnergyXtal);
   fChain->SetBranchAddress("ootPhoE2ndFull5x5", &ootPhoE2ndFull5x5, &b_ootPhoE2ndFull5x5);
   fChain->SetBranchAddress("ootPhoE1x3Full5x5", &ootPhoE1x3Full5x5, &b_ootPhoE1x3Full5x5);
   fChain->SetBranchAddress("ootPhoE1x5Full5x5", &ootPhoE1x5Full5x5, &b_ootPhoE1x5Full5x5);
   fChain->SetBranchAddress("ootPhoE2x5Full5x5", &ootPhoE2x5Full5x5, &b_ootPhoE2x5Full5x5);
   fChain->SetBranchAddress("ootPhoPFClusEcalIso", &ootPhoPFClusEcalIso, &b_ootPhoPFClusEcalIso);
   fChain->SetBranchAddress("ootPhoPFClusHcalIso", &ootPhoPFClusHcalIso, &b_ootPhoPFClusHcalIso);
   fChain->SetBranchAddress("nootPhoTrkHollowConeDR03", &nootPhoTrkHollowConeDR03, &b_nootPhoTrkHollowConeDR03);
   fChain->SetBranchAddress("nootPhoTrkSolidConeDR03", &nootPhoTrkSolidConeDR03, &b_nootPhoTrkSolidConeDR03);
   fChain->SetBranchAddress("ootPhoTrkSumPtSolidConeDR03", &ootPhoTrkSumPtSolidConeDR03, &b_ootPhoTrkSumPtSolidConeDR03);
   fChain->SetBranchAddress("ootPhoTrkSumPtHollowConeDR03", &ootPhoTrkSumPtHollowConeDR03, &b_ootPhoTrkSumPtHollowConeDR03);
   fChain->SetBranchAddress("nootPhoTrkHollowConeDR04", &nootPhoTrkHollowConeDR04, &b_nootPhoTrkHollowConeDR04);
   fChain->SetBranchAddress("nootPhoTrkSolidConeDR04", &nootPhoTrkSolidConeDR04, &b_nootPhoTrkSolidConeDR04);
   fChain->SetBranchAddress("ootPhoTrkSumPtSolidConeDR04", &ootPhoTrkSumPtSolidConeDR04, &b_ootPhoTrkSumPtSolidConeDR04);
   fChain->SetBranchAddress("ootPhoTrkSumPtHollowConeDR04", &ootPhoTrkSumPtHollowConeDR04, &b_ootPhoTrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("ootPhoECALIso", &ootPhoECALIso, &b_ootPhoECALIso);
   fChain->SetBranchAddress("ootPhoHCALIso", &ootPhoHCALIso, &b_ootPhoHCALIso);
   fChain->SetBranchAddress("ootPhoSeedBCE", &ootPhoSeedBCE, &b_ootPhoSeedBCE);
   fChain->SetBranchAddress("ootPhoSeedBCEta", &ootPhoSeedBCEta, &b_ootPhoSeedBCEta);
   fChain->SetBranchAddress("ootPhoSeedBCPhi", &ootPhoSeedBCPhi, &b_ootPhoSeedBCPhi);
   fChain->SetBranchAddress("ootPho_FiredSingleTrgs", &ootPho_FiredSingleTrgs, &b_ootPho_FiredSingleTrgs);
   fChain->SetBranchAddress("ootPho_FiredDoubleTrgs", &ootPho_FiredDoubleTrgs, &b_ootPho_FiredDoubleTrgs);
   fChain->SetBranchAddress("ootPho_FiredTripleTrgs", &ootPho_FiredTripleTrgs, &b_ootPho_FiredTripleTrgs);
   fChain->SetBranchAddress("ootPho_FiredL1Trgs", &ootPho_FiredL1Trgs, &b_ootPho_FiredL1Trgs);
   fChain->SetBranchAddress("ootPho_SeedTime", &ootPho_SeedTime, &b_ootPho_SeedTime);
   fChain->SetBranchAddress("ootPho_SeedEnergy", &ootPho_SeedEnergy, &b_ootPho_SeedEnergy);
   fChain->SetBranchAddress("ootPho_MIPChi2", &ootPho_MIPChi2, &b_ootPho_MIPChi2);
   fChain->SetBranchAddress("ootPho_MIPTotEnergy", &ootPho_MIPTotEnergy, &b_ootPho_MIPTotEnergy);
   fChain->SetBranchAddress("ootPho_MIPSlope", &ootPho_MIPSlope, &b_ootPho_MIPSlope);
   fChain->SetBranchAddress("ootPho_MIPIntercept", &ootPho_MIPIntercept, &b_ootPho_MIPIntercept);
   fChain->SetBranchAddress("ootPho_MIPNhitCone", &ootPho_MIPNhitCone, &b_ootPho_MIPNhitCone);
   fChain->SetBranchAddress("ootPho_IDbit", &ootPho_IDbit, &b_ootPho_IDbit);
   fChain->SetBranchAddress("ootPhoDirectEcalSCindex", &ootPhoDirectEcalSCindex, &b_ootPhoDirectEcalSCindex);
   fChain->SetBranchAddress("necalootSC", &necalootSC, &b_necalootSC);
   fChain->SetBranchAddress("ecalootSC_eta", &ecalootSC_eta, &b_ecalootSC_eta);
   fChain->SetBranchAddress("ecalootSC_phi", &ecalootSC_phi, &b_ecalootSC_phi);
   fChain->SetBranchAddress("ecalootSC_En", &ecalootSC_En, &b_ecalootSC_En);
   fChain->SetBranchAddress("ecalootSC_RawEn", &ecalootSC_RawEn, &b_ecalootSC_RawEn);
   fChain->SetBranchAddress("ecalootSC_etaWidth", &ecalootSC_etaWidth, &b_ecalootSC_etaWidth);
   fChain->SetBranchAddress("ecalootSC_phiWidth", &ecalootSC_phiWidth, &b_ecalootSC_phiWidth);
   fChain->SetBranchAddress("ecalootSC_LICTD", &ecalootSC_LICTD, &b_ecalootSC_LICTD);
   fChain->SetBranchAddress("ecalootSC_nL1Spike", &ecalootSC_nL1Spike, &b_ecalootSC_nL1Spike);
   fChain->SetBranchAddress("ecalootSC_nDiweird", &ecalootSC_nDiweird, &b_ecalootSC_nDiweird);
   fChain->SetBranchAddress("ecalootSC_nWeird", &ecalootSC_nWeird, &b_ecalootSC_nWeird);
   fChain->SetBranchAddress("ecalootSC_nSaturated", &ecalootSC_nSaturated, &b_ecalootSC_nSaturated);
   fChain->SetBranchAddress("ecalootSC_nOutOfTime", &ecalootSC_nOutOfTime, &b_ecalootSC_nOutOfTime);
   fChain->SetBranchAddress("ecalootSC_nXtals", &ecalootSC_nXtals, &b_ecalootSC_nXtals);
   fChain->SetBranchAddress("ecalootSC_maxEnXtalTime", &ecalootSC_maxEnXtalTime, &b_ecalootSC_maxEnXtalTime);
   fChain->SetBranchAddress("ecalootSC_maxEnXtalSwissCross", &ecalootSC_maxEnXtalSwissCross, &b_ecalootSC_maxEnXtalSwissCross);
   fChain->SetBranchAddress("ecalootSC_maxEnXtalBits", &ecalootSC_maxEnXtalBits, &b_ecalootSC_maxEnXtalBits);
   fChain->SetBranchAddress("ecalootSC_seedIx", &ecalootSC_seedIx, &b_ecalootSC_seedIx);
   fChain->SetBranchAddress("ecalootSC_seedIy", &ecalootSC_seedIy, &b_ecalootSC_seedIy);
   fChain->SetBranchAddress("ecalootSC_seedIz", &ecalootSC_seedIz, &b_ecalootSC_seedIz);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleEcalEn", &eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("elePtError", &elePtError, &b_elePtError);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleCalibPt", &eleCalibPt, &b_eleCalibPt);
   fChain->SetBranchAddress("eleCalibEn", &eleCalibEn, &b_eleCalibEn);
   fChain->SetBranchAddress("eleSCindex", &eleSCindex, &b_eleSCindex);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("eleQualityBits", &eleQualityBits, &b_eleQualityBits);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleESEffSigmaRR", &eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso, &b_elePFClusEcalIso);
   fChain->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso, &b_elePFClusHcalIso);
   fChain->SetBranchAddress("eleIDMVAIso", &eleIDMVAIso, &b_eleIDMVAIso);
   fChain->SetBranchAddress("eleIDMVANoIso", &eleIDMVANoIso, &b_eleIDMVANoIso);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleFiredSingleTrgs", &eleFiredSingleTrgs, &b_eleFiredSingleTrgs);
   fChain->SetBranchAddress("eleFiredDoubleTrgs", &eleFiredDoubleTrgs, &b_eleFiredDoubleTrgs);
   fChain->SetBranchAddress("eleFiredL1Trgs", &eleFiredL1Trgs, &b_eleFiredL1Trgs);
   fChain->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
   fChain->SetBranchAddress("eleGenIndex", &eleGenIndex, &b_eleGenIndex);
   fChain->SetBranchAddress("eleDirectEcalSCindex", &eleDirectEcalSCindex, &b_eleDirectEcalSCindex);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEn", &muEn, &b_muEn);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIDbit", &muIDbit, &b_muIDbit);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muSIP", &muSIP, &b_muSIP);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   fChain->SetBranchAddress("muFiredTrgs", &muFiredTrgs, &b_muFiredTrgs);
   fChain->SetBranchAddress("muFiredL1Trgs", &muFiredL1Trgs, &b_muFiredL1Trgs);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("muBestTrkType", &muBestTrkType, &b_muBestTrkType);
   fChain->SetBranchAddress("muGenIndex", &muGenIndex, &b_muGenIndex);
   fChain->SetBranchAddress("nAK4CHSJet", &nAK4CHSJet, &b_nAK4CHSJet);
   fChain->SetBranchAddress("AK4CHSJet_Charge", &AK4CHSJet_Charge, &b_AK4CHSJet_Charge);
   fChain->SetBranchAddress("AK4CHSJet_Pt", &AK4CHSJet_Pt, &b_AK4CHSJet_Pt);
   fChain->SetBranchAddress("AK4CHSJet_En", &AK4CHSJet_En, &b_AK4CHSJet_En);
   fChain->SetBranchAddress("AK4CHSJet_Eta", &AK4CHSJet_Eta, &b_AK4CHSJet_Eta);
   fChain->SetBranchAddress("AK4CHSJet_Phi", &AK4CHSJet_Phi, &b_AK4CHSJet_Phi);
   fChain->SetBranchAddress("AK4CHSJet_RawPt", &AK4CHSJet_RawPt, &b_AK4CHSJet_RawPt);
   fChain->SetBranchAddress("AK4CHSJet_RawEn", &AK4CHSJet_RawEn, &b_AK4CHSJet_RawEn);
   fChain->SetBranchAddress("AK4CHSJet_Mt", &AK4CHSJet_Mt, &b_AK4CHSJet_Mt);
   fChain->SetBranchAddress("AK4CHSJet_Area", &AK4CHSJet_Area, &b_AK4CHSJet_Area);
   fChain->SetBranchAddress("AK4CHSJet_LeadTrackPt", &AK4CHSJet_LeadTrackPt, &b_AK4CHSJet_LeadTrackPt);
   fChain->SetBranchAddress("AK4CHSJet_LeadTrackEta", &AK4CHSJet_LeadTrackEta, &b_AK4CHSJet_LeadTrackEta);
   fChain->SetBranchAddress("AK4CHSJet_LeadTrackPhi", &AK4CHSJet_LeadTrackPhi, &b_AK4CHSJet_LeadTrackPhi);
   fChain->SetBranchAddress("AK4CHSJet_LepTrackPID", &AK4CHSJet_LepTrackPID, &b_AK4CHSJet_LepTrackPID);
   fChain->SetBranchAddress("AK4CHSJet_LepTrackPt", &AK4CHSJet_LepTrackPt, &b_AK4CHSJet_LepTrackPt);
   fChain->SetBranchAddress("AK4CHSJet_LepTrackEta", &AK4CHSJet_LepTrackEta, &b_AK4CHSJet_LepTrackEta);
   fChain->SetBranchAddress("AK4CHSJet_LepTrackPhi", &AK4CHSJet_LepTrackPhi, &b_AK4CHSJet_LepTrackPhi);
   fChain->SetBranchAddress("AK4CHSJet_ConsituentEtaPhiSpread", &AK4CHSJet_ConsituentEtaPhiSpread, &b_AK4CHSJet_ConsituentEtaPhiSpread);
   fChain->SetBranchAddress("AK4CHSJet_constituentPtDistribution", &AK4CHSJet_constituentPtDistribution, &b_AK4CHSJet_constituentPtDistribution);
   fChain->SetBranchAddress("AK4CHSJet_MaxConstituentDistance", &AK4CHSJet_MaxConstituentDistance, &b_AK4CHSJet_MaxConstituentDistance);
   fChain->SetBranchAddress("AK4CHSJet_CSV2BJetTags", &AK4CHSJet_CSV2BJetTags, &b_AK4CHSJet_CSV2BJetTags);
   fChain->SetBranchAddress("AK4CHSJet_DeepCSVTags_b", &AK4CHSJet_DeepCSVTags_b, &b_AK4CHSJet_DeepCSVTags_b);
   fChain->SetBranchAddress("AK4CHSJet_DeepCSVTags_bb", &AK4CHSJet_DeepCSVTags_bb, &b_AK4CHSJet_DeepCSVTags_bb);
   fChain->SetBranchAddress("AK4CHSJet_DeepCSVTags_c", &AK4CHSJet_DeepCSVTags_c, &b_AK4CHSJet_DeepCSVTags_c);
   fChain->SetBranchAddress("AK4CHSJet_DeepCSVTags_udsg", &AK4CHSJet_DeepCSVTags_udsg, &b_AK4CHSJet_DeepCSVTags_udsg);
   fChain->SetBranchAddress("AK4CHSJet_DeepFlavTags_bb", &AK4CHSJet_DeepFlavTags_bb, &b_AK4CHSJet_DeepFlavTags_bb);
   fChain->SetBranchAddress("AK4CHSJet_DeepFlavTags_bbb", &AK4CHSJet_DeepFlavTags_bbb, &b_AK4CHSJet_DeepFlavTags_bbb);
   fChain->SetBranchAddress("AK4CHSJet_DeepFlavTags_lepb", &AK4CHSJet_DeepFlavTags_lepb, &b_AK4CHSJet_DeepFlavTags_lepb);
   fChain->SetBranchAddress("AK4CHSJet_DeepFlavTags_c", &AK4CHSJet_DeepFlavTags_c, &b_AK4CHSJet_DeepFlavTags_c);
   fChain->SetBranchAddress("AK4CHSJet_DeepFlavTags_uds", &AK4CHSJet_DeepFlavTags_uds, &b_AK4CHSJet_DeepFlavTags_uds);
   fChain->SetBranchAddress("AK4CHSJet_DeepFlavTags_g", &AK4CHSJet_DeepFlavTags_g, &b_AK4CHSJet_DeepFlavTags_g);
   fChain->SetBranchAddress("AK4CHSJet_CombMVA2Tags", &AK4CHSJet_CombMVA2Tags, &b_AK4CHSJet_CombMVA2Tags);
   fChain->SetBranchAddress("AK4CHSJet_PartonFlavour", &AK4CHSJet_PartonFlavour, &b_AK4CHSJet_PartonFlavour);
   fChain->SetBranchAddress("AK4CHSJet_HadronFlavour", &AK4CHSJet_HadronFlavour, &b_AK4CHSJet_HadronFlavour);
   fChain->SetBranchAddress("AK4CHSJet_GenJetIndex", &AK4CHSJet_GenJetIndex, &b_AK4CHSJet_GenJetIndex);
   fChain->SetBranchAddress("AK4CHSJet_GenPartonIndex", &AK4CHSJet_GenPartonIndex, &b_AK4CHSJet_GenPartonIndex);
   fChain->SetBranchAddress("AK4CHSJet_ID", &AK4CHSJet_ID, &b_AK4CHSJet_ID);
   fChain->SetBranchAddress("AK4CHSJet_PUID", &AK4CHSJet_PUID, &b_AK4CHSJet_PUID);
   fChain->SetBranchAddress("AK4CHSJet_PUFullID", &AK4CHSJet_PUFullID, &b_AK4CHSJet_PUFullID);
   fChain->SetBranchAddress("AK4CHSJet_qgLikelihood", &AK4CHSJet_qgLikelihood, &b_AK4CHSJet_qgLikelihood);
   fChain->SetBranchAddress("AK4CHSJet_JECUnc", &AK4CHSJet_JECUnc, &b_AK4CHSJet_JECUnc);
   fChain->SetBranchAddress("AK4CHSJet_JER_ptRes", &AK4CHSJet_JER_ptRes, &b_AK4CHSJet_JER_ptRes);
   fChain->SetBranchAddress("AK4CHSJet_JER_phiRes", &AK4CHSJet_JER_phiRes, &b_AK4CHSJet_JER_phiRes);
   fChain->SetBranchAddress("AK4CHSJet_JER_SF", &AK4CHSJet_JER_SF, &b_AK4CHSJet_JER_SF);
   fChain->SetBranchAddress("AK4CHSJet_JER_SFup", &AK4CHSJet_JER_SFup, &b_AK4CHSJet_JER_SFup);
   fChain->SetBranchAddress("AK4CHSJet_JER_SFdown", &AK4CHSJet_JER_SFdown, &b_AK4CHSJet_JER_SFdown);
   fChain->SetBranchAddress("AK4CHSJet_FiredTrgs", &AK4CHSJet_FiredTrgs, &b_AK4CHSJet_FiredTrgs);
   fChain->SetBranchAddress("AK4CHSJet_CHF", &AK4CHSJet_CHF, &b_AK4CHSJet_CHF);
   fChain->SetBranchAddress("AK4CHSJet_NHF", &AK4CHSJet_NHF, &b_AK4CHSJet_NHF);
   fChain->SetBranchAddress("AK4CHSJet_CEF", &AK4CHSJet_CEF, &b_AK4CHSJet_CEF);
   fChain->SetBranchAddress("AK4CHSJet_NEF", &AK4CHSJet_NEF, &b_AK4CHSJet_NEF);
   fChain->SetBranchAddress("AK4CHSJet_NCH", &AK4CHSJet_NCH, &b_AK4CHSJet_NCH);
   fChain->SetBranchAddress("AK4CHSJet_NNP", &AK4CHSJet_NNP, &b_AK4CHSJet_NNP);
   fChain->SetBranchAddress("AK4CHSJet_MUF", &AK4CHSJet_MUF, &b_AK4CHSJet_MUF);
   Notify();
}

Bool_t Analyse::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyse::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyse::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analyse_cxx

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 24 09:58:11 2020 by ROOT version 6.23/01
// from TTree EventTree/Intermediate TTree
// found on file: ztonunu.root
//////////////////////////////////////////////////////////

#ifndef Analyse_h
#define Analyse_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Analyse {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Long64_t        event;
   Int_t           isPVGood;
   Int_t           nVtx;
   Float_t         genWeight;
   Float_t         puweisj;
   Float_t         puweiup;
   Float_t         puweidown;
   Float_t         mcwei;
   Int_t           isphoHLT175;
   Int_t           isphoHLT200;
   Int_t           nPho;
   Float_t         phoPt[9];   //[nPho]
   Float_t         phoCorrPt[9];   //[nPho]
   Float_t         phoEta[9];   //[nPho]
   Float_t         phoPhi[9];   //[nPho]
   Float_t         phoSCEta[9];   //[nPho]
   Float_t         phoSCPhi[9];   //[nPho]
   Float_t         phoR9[9];   //[nPho]
   Float_t         phoMVA[9];   //[nPho]
   Float_t         phoIDLoose[9];   //[nPho]
   Float_t         phogenPt[9];   //[nPho]
   Float_t         phoReso[9];   //[nPho]
   Float_t         phoIDMedium[9];   //[nPho]
   Float_t         phoIDTight[9];   //[nPho]
   Float_t         phoIDMVA[9];   //[nPho]
   Float_t         phoIDhighpt[9];   //[nPho]
   Int_t           phoPresel[9];   //[nPho]
   Int_t           phoisMCmatch[9];   //[nPho]
   Int_t           phoMCID[9];   //[nPho]
   Int_t           phoMCmomID[9];   //[nPho]
   Float_t         phogenEta[9];   //[nPho]
   Float_t         phogenPhi[9];   //[nPho]
   Int_t           keepEv;
   Float_t         hasDirectPromptPho;
   Float_t         phoHoverE[9];   //[nPho]
   Float_t         phoSieie[9];   //[nPho]
   Float_t         phoChwiso[9];   //[nPho]
   Float_t         phoNhiso[9];   //[nPho]
   Float_t         phoPiso[9];   //[nPho]
   Float_t         phoSeedTime[9];   //[nPho]
   Float_t         phoSipip[9];   //[nPho]
   Float_t         phoSieip[9];   //[nPho]
   Float_t         phoIDmonopho[9];   //[nPho]
   Bool_t          phoeleVeto[9];   //[nPho]
   Bool_t          phopixVeto[9];   //[nPho]
   Float_t         phoMIPChi2[9];   //[nPho]
   Float_t         phoMIPTotEnergy[9];   //[nPho]
   Float_t         phoMIPSlope[9];   //[nPho]
   Float_t         phoMIPIntercept[9];   //[nPho]
   Float_t         phoMIPNhitCone[9];   //[nPho]
   Float_t         phoESEffSRR[9];   //[nPho]
   Float_t         phoESEnP1[9];   //[nPho]
   Float_t         phoESEnP2[9];   //[nPho]
   Float_t         phoSCEtaW[9];   //[nPho]
   Float_t         phoSCPhiW[9];   //[nPho]
   Float_t         phoSCBrem[9];   //[nPho]
   Float_t         phoSCRawE[9];   //[nPho]
   Float_t         phoSCE[9];   //[nPho]
   Int_t           metFilters;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMET_metSig;
   Float_t         pfMET_EtSig;
   Int_t           nEle;
   Float_t         elePt[8];   //[nEle]
   Float_t         eleCorrPt[8];   //[nEle]
   Float_t         eleEta[8];   //[nEle]
   Float_t         elePhi[8];   //[nEle]
   Float_t         eleR9Full5x5[8];   //[nEle]
   Float_t         eleSCEta[8];   //[nEle]
   Float_t         eleSCPhi[8];   //[nEle]
   Float_t         eleEnRes[8];   //[nEle]
   Float_t         eleEnResRegr[8];   //[nEle]
   Int_t           eleCharge[8];   //[nEle]
   Float_t         eleIDVeto[8];   //[nEle]
   Float_t         eleIDLoose[8];   //[nEle]
   Float_t         eleIDMedium[8];   //[nEle]
   Float_t         eleIDTight[8];   //[nEle]
   Float_t         eleIDMVATrig[8];   //[nEle]
   Float_t         eleIDHZZMVA[8];   //[nEle]
   Float_t         eleID3Lepton[8];   //[nEle]
   Float_t         elePassPresel[8];   //[nEle]
   Int_t           eleisMCmatch[8];   //[nEle]
   Int_t           eleMCID[8];   //[nEle]
   Int_t           eleMCmomID[8];   //[nEle]
   Float_t         elegenPt[8];   //[nEle]
   Float_t         elegenEta[8];   //[nEle]
   Float_t         elegenPhi[8];   //[nEle]
   Int_t           nMu;
   Int_t           muNtrk[5];   //[nMu]
   Float_t         muPt[5];   //[nMu]
   Float_t         muCorrPt[5];   //[nMu]
   Float_t         muEta[5];   //[nMu]
   Float_t         muPhi[5];   //[nMu]
   Float_t         muEnRes[5];   //[nMu]
   Float_t         muCharge[5];   //[nMu]
   Float_t         muIDLoose[5];   //[nMu]
   Float_t         muIDTight[5];   //[nMu]
   Float_t         muIDVeryTight[5];   //[nMu]
   Float_t         muIDHZZ[5];   //[nMu]
   Int_t           muisMCmatch[5];   //[nMu]
   Int_t           muMCID[5];   //[nMu]
   Int_t           muMCmomID[5];   //[nMu]
   Float_t         mugenPt[5];   //[nMu]
   Float_t         mugenEta[5];   //[nMu]
   Float_t         mugenPhi[5];   //[nMu]
   Int_t           nJet;
   Float_t         jetPt[1];   //[nJet]
   Float_t         jetEta[1];   //[nJet]
   Float_t         jetPhi[1];   //[nJet]
   Float_t         jetEn[1];   //[nJet]
   Int_t           nGen;
   Float_t         genEta[9];   //[nGen]
   Float_t         genPhi[9];   //[nGen]
   Float_t         genPt[9];   //[nGen]
   Float_t         genE[9];   //[nGen]
   Int_t           genID[9];   //[nGen]
   Int_t           genStatus[9];   //[nGen]
   Int_t           nhfEle;
   Float_t         hfelePt[6];   //[nhfEle]
   Float_t         hfeleEn[6];   //[nhfEle]
   Float_t         hfeleEta[6];   //[nhfEle]
   Float_t         hfelePhi[6];   //[nhfEle]
   Float_t         hfeleIso[6];   //[nhfEle]
   Float_t         hfeleECALEn[6];   //[nhfEle]
   Float_t         hfeleHCALEn[6];   //[nhfEle]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_isPVGood;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_puweisj;   //!
   TBranch        *b_puweiup;   //!
   TBranch        *b_puweidown;   //!
   TBranch        *b_mcwei;   //!
   TBranch        *b_isphoHLT175;   //!
   TBranch        *b_isphoHLT200;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoPt;   //!
   TBranch        *b_phoCorrPt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoMVA;   //!
   TBranch        *b_phoIDLoose;   //!
   TBranch        *b_phogenPt;   //!
   TBranch        *b_phoReso;   //!
   TBranch        *b_phoIDMedium;   //!
   TBranch        *b_phoIDTight;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoIDhighpt;   //!
   TBranch        *b_phoPresel;   //!
   TBranch        *b_phoisMCmatch;   //!
   TBranch        *b_phoMCID;   //!
   TBranch        *b_phoMCmomID;   //!
   TBranch        *b_phogenEta;   //!
   TBranch        *b_phogenPhi;   //!
   TBranch        *b_keepEv;   //!
   TBranch        *b_hasDirectPromptPho;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSieie;   //!
   TBranch        *b_phoChwiso;   //!
   TBranch        *b_phoNhiso;   //!
   TBranch        *b_phoPiso;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSipip;   //!
   TBranch        *b_phoSieip;   //!
   TBranch        *b_phoIDmonopho;   //!
   TBranch        *b_phoeleVeto;   //!
   TBranch        *b_phopixVeto;   //!
   TBranch        *b_phoMIPChi2;   //!
   TBranch        *b_phoMIPTotEnergy;   //!
   TBranch        *b_phoMIPSlope;   //!
   TBranch        *b_phoMIPIntercept;   //!
   TBranch        *b_phoMIPNhitCone;   //!
   TBranch        *b_phoESEffSRR;   //!
   TBranch        *b_phoESEnP1;   //!
   TBranch        *b_phoESEnP2;   //!
   TBranch        *b_phoSCEtaW;   //!
   TBranch        *b_phoSCPhiW;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMET_metSig;   //!
   TBranch        *b_pfMET_EtSig;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleCorrPt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleEnRes;   //!
   TBranch        *b_eleEnResRegr;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleIDVeto;   //!
   TBranch        *b_eleIDLoose;   //!
   TBranch        *b_eleIDMedium;   //!
   TBranch        *b_eleIDTight;   //!
   TBranch        *b_eleIDMVATrig;   //!
   TBranch        *b_eleIDHZZMVA;   //!
   TBranch        *b_eleID3Lepton;   //!
   TBranch        *b_elePassPresel;   //!
   TBranch        *b_eleisMCmatch;   //!
   TBranch        *b_eleMCID;   //!
   TBranch        *b_eleMCmomID;   //!
   TBranch        *b_elegenPt;   //!
   TBranch        *b_elegenEta;   //!
   TBranch        *b_elegenPhi;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muNtrk;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muCorrPt;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muEnRes;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muIDLoose;   //!
   TBranch        *b_muIDTight;   //!
   TBranch        *b_muIDVeryTight;   //!
   TBranch        *b_muIDHZZ;   //!
   TBranch        *b_muisMCmatch;   //!
   TBranch        *b_muMCID;   //!
   TBranch        *b_muMCmomID;   //!
   TBranch        *b_mugenPt;   //!
   TBranch        *b_mugenEta;   //!
   TBranch        *b_mugenPhi;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_nGen;   //!
   TBranch        *b_genEta;   //!
   TBranch        *b_genPhi;   //!
   TBranch        *b_genPt;   //!
   TBranch        *b_genE;   //!
   TBranch        *b_genID;   //!
   TBranch        *b_genStatus;   //!
   TBranch        *b_nhfEle;   //!
   TBranch        *b_hfelePt;   //!
   TBranch        *b_hfeleEn;   //!
   TBranch        *b_hfeleEta;   //!
   TBranch        *b_hfelePhi;   //!
   TBranch        *b_hfeleIso;   //!
   TBranch        *b_hfeleECALEn;   //!
   TBranch        *b_hfeleHCALEn;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("wgjets.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("wgjets.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("wgjets.root:/ggNtuplizer");
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

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("isPVGood", &isPVGood, &b_isPVGood);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("puweisj", &puweisj, &b_puweisj);
   fChain->SetBranchAddress("puweiup", &puweiup, &b_puweiup);
   fChain->SetBranchAddress("puweidown", &puweidown, &b_puweidown);
   fChain->SetBranchAddress("mcwei", &mcwei, &b_mcwei);
   fChain->SetBranchAddress("isphoHLT175", &isphoHLT175, &b_isphoHLT175);
   fChain->SetBranchAddress("isphoHLT200", &isphoHLT200, &b_isphoHLT200);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoPt", phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoCorrPt", phoCorrPt, &b_phoCorrPt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoSCEta", phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoMVA", phoMVA, &b_phoMVA);
   fChain->SetBranchAddress("phoIDLoose", phoIDLoose, &b_phoIDLoose);
   fChain->SetBranchAddress("phogenPt", phogenPt, &b_phogenPt);
   fChain->SetBranchAddress("phoReso", phoReso, &b_phoReso);
   fChain->SetBranchAddress("phoIDMedium", phoIDMedium, &b_phoIDMedium);
   fChain->SetBranchAddress("phoIDTight", phoIDTight, &b_phoIDTight);
   fChain->SetBranchAddress("phoIDMVA", phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoIDhighpt", phoIDhighpt, &b_phoIDhighpt);
   fChain->SetBranchAddress("phoPresel", phoPresel, &b_phoPresel);
   fChain->SetBranchAddress("phoisMCmatch", phoisMCmatch, &b_phoisMCmatch);
   fChain->SetBranchAddress("phoMCID", phoMCID, &b_phoMCID);
   fChain->SetBranchAddress("phoMCmomID", phoMCmomID, &b_phoMCmomID);
   fChain->SetBranchAddress("phogenEta", phogenEta, &b_phogenEta);
   fChain->SetBranchAddress("phogenPhi", phogenPhi, &b_phogenPhi);
   fChain->SetBranchAddress("keepEv", &keepEv, &b_keepEv);
   fChain->SetBranchAddress("hasDirectPromptPho", &hasDirectPromptPho, &b_hasDirectPromptPho);
   fChain->SetBranchAddress("phoHoverE", phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSieie", phoSieie, &b_phoSieie);
   fChain->SetBranchAddress("phoChwiso", phoChwiso, &b_phoChwiso);
   fChain->SetBranchAddress("phoNhiso", phoNhiso, &b_phoNhiso);
   fChain->SetBranchAddress("phoPiso", phoPiso, &b_phoPiso);
   fChain->SetBranchAddress("phoSeedTime", phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSipip", phoSipip, &b_phoSipip);
   fChain->SetBranchAddress("phoSieip", phoSieip, &b_phoSieip);
   fChain->SetBranchAddress("phoIDmonopho", phoIDmonopho, &b_phoIDmonopho);
   fChain->SetBranchAddress("phoeleVeto", phoeleVeto, &b_phoeleVeto);
   fChain->SetBranchAddress("phopixVeto", phopixVeto, &b_phopixVeto);
   fChain->SetBranchAddress("phoMIPChi2", phoMIPChi2, &b_phoMIPChi2);
   fChain->SetBranchAddress("phoMIPTotEnergy", phoMIPTotEnergy, &b_phoMIPTotEnergy);
   fChain->SetBranchAddress("phoMIPSlope", phoMIPSlope, &b_phoMIPSlope);
   fChain->SetBranchAddress("phoMIPIntercept", phoMIPIntercept, &b_phoMIPIntercept);
   fChain->SetBranchAddress("phoMIPNhitCone", phoMIPNhitCone, &b_phoMIPNhitCone);
   fChain->SetBranchAddress("phoESEffSRR", phoESEffSRR, &b_phoESEffSRR);
   fChain->SetBranchAddress("phoESEnP1", phoESEnP1, &b_phoESEnP1);
   fChain->SetBranchAddress("phoESEnP2", phoESEnP2, &b_phoESEnP2);
   fChain->SetBranchAddress("phoSCEtaW", phoSCEtaW, &b_phoSCEtaW);
   fChain->SetBranchAddress("phoSCPhiW", phoSCPhiW, &b_phoSCPhiW);
   fChain->SetBranchAddress("phoSCBrem", phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phoSCRawE", phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoSCE", phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMET_metSig", &pfMET_metSig, &b_pfMET_metSig);
   fChain->SetBranchAddress("pfMET_EtSig", &pfMET_EtSig, &b_pfMET_EtSig);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("eleCorrPt", eleCorrPt, &b_eleCorrPt);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9Full5x5", eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleEnRes", eleEnRes, &b_eleEnRes);
   fChain->SetBranchAddress("eleEnResRegr", eleEnResRegr, &b_eleEnResRegr);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleIDVeto", eleIDVeto, &b_eleIDVeto);
   fChain->SetBranchAddress("eleIDLoose", eleIDLoose, &b_eleIDLoose);
   fChain->SetBranchAddress("eleIDMedium", eleIDMedium, &b_eleIDMedium);
   fChain->SetBranchAddress("eleIDTight", eleIDTight, &b_eleIDTight);
   fChain->SetBranchAddress("eleIDMVATrig", eleIDMVATrig, &b_eleIDMVATrig);
   fChain->SetBranchAddress("eleIDHZZMVA", eleIDHZZMVA, &b_eleIDHZZMVA);
   fChain->SetBranchAddress("eleID3Lepton", eleID3Lepton, &b_eleID3Lepton);
   fChain->SetBranchAddress("elePassPresel", elePassPresel, &b_elePassPresel);
   fChain->SetBranchAddress("eleisMCmatch", eleisMCmatch, &b_eleisMCmatch);
   fChain->SetBranchAddress("eleMCID", eleMCID, &b_eleMCID);
   fChain->SetBranchAddress("eleMCmomID", eleMCmomID, &b_eleMCmomID);
   fChain->SetBranchAddress("elegenPt", elegenPt, &b_elegenPt);
   fChain->SetBranchAddress("elegenEta", elegenEta, &b_elegenEta);
   fChain->SetBranchAddress("elegenPhi", elegenPhi, &b_elegenPhi);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muNtrk", muNtrk, &b_muNtrk);
   fChain->SetBranchAddress("muPt", muPt, &b_muPt);
   fChain->SetBranchAddress("muCorrPt", muCorrPt, &b_muCorrPt);
   fChain->SetBranchAddress("muEta", muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", muPhi, &b_muPhi);
   fChain->SetBranchAddress("muEnRes", muEnRes, &b_muEnRes);
   fChain->SetBranchAddress("muCharge", muCharge, &b_muCharge);
   fChain->SetBranchAddress("muIDLoose", muIDLoose, &b_muIDLoose);
   fChain->SetBranchAddress("muIDTight", muIDTight, &b_muIDTight);
   fChain->SetBranchAddress("muIDVeryTight", muIDVeryTight, &b_muIDVeryTight);
   fChain->SetBranchAddress("muIDHZZ", muIDHZZ, &b_muIDHZZ);
   fChain->SetBranchAddress("muisMCmatch", muisMCmatch, &b_muisMCmatch);
   fChain->SetBranchAddress("muMCID", muMCID, &b_muMCID);
   fChain->SetBranchAddress("muMCmomID", muMCmomID, &b_muMCmomID);
   fChain->SetBranchAddress("mugenPt", mugenPt, &b_mugenPt);
   fChain->SetBranchAddress("mugenEta", mugenEta, &b_mugenEta);
   fChain->SetBranchAddress("mugenPhi", mugenPhi, &b_mugenPhi);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("nGen", &nGen, &b_nGen);
   fChain->SetBranchAddress("genEta", genEta, &b_genEta);
   fChain->SetBranchAddress("genPhi", genPhi, &b_genPhi);
   fChain->SetBranchAddress("genPt", genPt, &b_genPt);
   fChain->SetBranchAddress("genE", genE, &b_genE);
   fChain->SetBranchAddress("genID", genID, &b_genID);
   fChain->SetBranchAddress("genStatus", genStatus, &b_genStatus);
   fChain->SetBranchAddress("nhfEle", &nhfEle, &b_nhfEle);
   fChain->SetBranchAddress("hfelePt", hfelePt, &b_hfelePt);
   fChain->SetBranchAddress("hfeleEn", hfeleEn, &b_hfeleEn);
   fChain->SetBranchAddress("hfeleEta", hfeleEta, &b_hfeleEta);
   fChain->SetBranchAddress("hfelePhi", hfelePhi, &b_hfelePhi);
   fChain->SetBranchAddress("hfeleIso", hfeleIso, &b_hfeleIso);
   fChain->SetBranchAddress("hfeleECALEn", hfeleECALEn, &b_hfeleECALEn);
   fChain->SetBranchAddress("hfeleHCALEn", hfeleHCALEn, &b_hfeleHCALEn);
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

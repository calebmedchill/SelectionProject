#define Analyse_cxx
#include "Analyse.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Analyse::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Analyse.C
//      root> Analyse t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
	double countVeto = 0; 
	double countEle = 0;
	double occPerEvent = 0;
	double occZero = 0;
	double occOne = 0;
	double occTwo = 0;
	double occMore = 0;
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		for (int i =0; i < nGen; i++){
			//cout<<"genID: " << abs(genID[i]) << endl;
			if (abs(genID[i]) == 11){
				countEle++;
				occPerEvent = 0;
				for (int j = 0; j <nGen; j++){
					if (abs(genID[j]) == 22 && (genPt[j] > 15 && genPt[j] < 200)){
						countVeto++;
						cout<< 	"Current Count:" << countVeto << endl;
						occPerEvent++;
					}
				}
				if (occPerEvent == 0){
					occZero++;
				}
				if(occPerEvent == 1){
					occOne++;
				}
				if (occPerEvent == 2){
					occTwo++;
				}
				if (occPerEvent > 2){
					occMore++;
				}
				
			}
		}
	}
	cout << "countEle: " << countEle << " countVeto: " << countVeto << " Ratio: " << countVeto/countEle << endl;
	cout << "Events with 0: " << occZero << endl;
	cout << "Events with 1: " << occOne << endl;
	cout << "Events with 2: " << occTwo << endl;
	cout << "Events with >2: " << occMore << endl;
}

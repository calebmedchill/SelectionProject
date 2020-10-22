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
	
	//Type casting the variable and assinging it to a new variable and trying to use that
	vector<float> *elePtAc = (vector<float> *)elePt;
	double veto = 0;
	double total = 0;
	//cout << *elePt.Pop() << " " << endl;


	//Dereferencing and then calling size on the array
	//cout << *elePt.size() << " " << endl;	
	
	//Dereferencing and then calling size on the array
	//cout << *elePt->size() << " " << endl;		

	//This was to see if it was an array.	
	//cout << elePt[0] << " " << endl;

	//This works but gives a size of 0
	//cout << eleEn->size() << " " << endl;
		
	//Using iterators
	//cout << *elePt->begin() << endl;
	//cout << *elePt->end() << endl;

	//vector<float> elePtVec = *elePT; 
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		//cout << elePt->size() << " " << endl;		
		cout << "Vector Size: " << phoEt->size() << " " << endl;		
		//This loop does not end		
				
		//cout << phoE->at(0) << " " << endl;
		total = total + 1;	
		for (auto i = 0; i < phoEt->size(); i++){
			//cout << "Count: " << i << " "<< endl;
			cout << phoEt->at(i) << " "<< endl;
			if(phoEt->at(i) > 15 && nEle>0){
				veto=veto + 1;
				i=1000;
			}
		}
		// if (Cut(ientry) < 0) continue;
	}
	cout<< "Total Events: " << total << " Vetoed Events: " << veto << " Ratio: " << veto/total << endl; 
}

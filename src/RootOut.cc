#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "RootOut.hh"
using namespace std;
RootOut::RootOut(string name0) {
	TTree::SetMaxTreeSize(100000000000LL);      // Set maximum size of the tree file
                                              // to 100 GB (instead of 1.9 GB).
	pointerToRoot=this;
	strcpy(rootOutputDirectoryName,"data");
	file_in = name0;
}

RootOut::~RootOut() {}

RootOut* RootOut::pointerToRoot=0;
RootOut* RootOut::GetRootInstance() {
	return pointerToRoot;
}

void RootOut::BeginOfRunAction(){
	string name1;
	G4cout << "musrRootOutput::BeginOfRunAction()  Defining the Root tree and branches:"<<G4endl;
	G4int tmpRunNr=(G4RunManager::GetRunManager())->GetCurrentRun()->GetRunID();
	char RootOutputFileName[200];
	//  sprintf(RootOutputFileName, "data/musr_%i.root", tmpRunNr);
	sprintf(RootOutputFileName, "%s/pal_%i.root",rootOutputDirectoryName,tmpRunNr);
	name1 = RootOutputFileName;
	cout<<name1<<endl;
	fstream filein(file_in,ios::in);
	string line;
	E_thr = 0.3;
	while(getline(filein,line)){
		istringstream sin(line);
		string command;
		sin>>command;
		if(command == "#/file_out"){
			sin>>name1;
		}
		if(command == "#/energy_thr"){
			sin>>E_thr;
		}
	}
	//cout<<E_thr<<endl;
	rootFile=new TFile(name1.data(),"recreate");
	
	rootTree=new TTree("t1","a simple Tree with simple variables");
	//BRANCH
	rootTree->Branch("Run_id",&Run_id,"Run_id/I");
	rootTree->Branch("Event_id",&Event_id,"Event_id/I");
	rootTree->Branch("Scint_n",&Scint_n,"Scint_n/I");
	rootTree->Branch("Scint_CopyN",Scint_CopyN,"Scint_CopyN[Scint_n]/I");
	rootTree->Branch("Edep",Edep,"Edep[Scint_n]/D");
	rootTree->Branch("Op_detect_num",&Op_detect_num,"Op_detect_num/I");
	rootTree->Branch("Op_num_0",&Op_num_0,"Op_num_0/I");
	rootTree->Branch("Op_num_1",&Op_num_1,"Op_num_1/I");
	rootTree->Branch("Anni_Pos_P",&Anni_Pos_P,"Anni_Pos_P/I");
	rootTree->Branch("kin_n",&kin_n,"kin_n/I");
	rootTree->Branch("Pos_P",Pos_P,"Pos_P[kin_n]/I");
	
	
	kin_n = 2;
	/*rootTree->Branch("kin_n",&kin_n,"kin_n/I");
	rootTree->Branch("Kin_E",Kin_E,"Kin_E[kin_n]/D");
	
	rootTree->Branch("Anni_Pos_X",&Anni_Pos_X,"Anni_Pos_X/D");
	rootTree->Branch("Anni_Pos_Y",&Anni_Pos_Y,"Anni_Pos_Y/D");
	rootTree->Branch("Anni_Pos_Z",&Anni_Pos_Z,"Anni_Pos_Z/D");
	rootTree->Branch("Anni_Pos_P",&Anni_Pos_P,"Anni_Pos_P/I");
	rootTree->Branch("Pos_P_X",Pos_P_X,"Pos_P_X[kin_n]/D");
	rootTree->Branch("Pos_P_Y",Pos_P_Y,"Pos_P_Y[kin_n]/D");
	rootTree->Branch("Pos_P_Z",Pos_P_Z,"Pos_P_Z[kin_n]/D");*/
}

void RootOut::EndOfRunAction(){
	rootFile->Write();
	rootFile->Close();
}

void RootOut::FillEvent(){
	//G4cout<<Edep[0]<<G4endl;
	rootTree->Fill();
	//cout<<Anni_Pos_P<<endl;
}

int RootOut::E_judge(){
	if(Edep[2]>=0.4&&Edep[2]<=0.6&&Edep[3]>=1.0&&Edep[3]<=1.42){
		return 1;
	}
	else if(Edep[3]>=0.4&&Edep[3]<=0.6&&Edep[2]>=1.0&&Edep[2]<=1.42){
		return 1;
	}
	else return 0;
}

void RootOut::Reset(){
	//Edep[0] = 0;
	//Edep[1] = 0;
	Pos_P[0] = 0;
	Pos_P[1] = 0;
	//kin_n = 0;
	Op_detect_num = 0;
	Op_num_0 = 0;
	Op_num_1 = 0;
}
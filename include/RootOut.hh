#ifndef RootOut_h
#define RootOut_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
//  ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
//
#include <map>

#include "G4ios.hh"
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

class RootOut{
	public:
	RootOut(string);
	~RootOut();
	static RootOut* GetRootInstance();
	
	void BeginOfRunAction();
    void EndOfRunAction();
	void FillEvent();
    void ClearAllRootVariables();
	
	
	private:
    TFile* rootFile;
    TTree* rootTree;
    static RootOut* pointerToRoot;
	string file_in;
	
	char   rootOutputDirectoryName[1000];
	
	//BRANCH
	static const Int_t num_Max = 200;  //闪烁探测器数量最大值 
	
	Int_t Run_id;
	Int_t Event_id;
	
	Int_t Scint_n; //闪烁体数量
	Int_t Scint_CopyN[num_Max]; //闪烁体copynum
	Double_t Edep[num_Max]; //闪烁体沉积能量
	Double_t Anni_Pos_X; //正电子湮没位置
	Double_t Anni_Pos_Y;
	Double_t Anni_Pos_Z;
	Int_t kin_n;//正电子过程数量
	Int_t Pos_P[1000];
	Double_t Kin_E[1000];
	Double_t Pos_P_X[1000];
	Double_t Pos_P_Y[1000];
	Double_t Pos_P_Z[1000];
	Int_t Anni_Pos_P;
	double E_thr;
	Int_t Op_detect_num; //探测光子数
	Int_t Op_num_0;
	Int_t Op_num_1;
	
	
	public:
	void SetRun_id(G4int id) {Run_id = id;}
	void SetEvent_id(G4int id) {Event_id = id;}
	void SetAnni_Pos(G4ThreeVector pol) {
		Anni_Pos_X = pol.x();
		Anni_Pos_Y = pol.y();
		Anni_Pos_Z = pol.z();
	}
	void SetPos(G4ThreeVector pol) {
		Pos_P_X[kin_n] = pol.x();
		Pos_P_Y[kin_n] = pol.y();
		Pos_P_Z[kin_n] = pol.z();
	}
	void Set_PosP(G4int a, G4int id) {
		Pos_P[a] = id;
	}
	void Set_Kin_E(G4double id) {
		Kin_E[kin_n] = id;
	}
	void kn_pp() {
		kin_n++;
	}
	void Op_pp_0(){
		Op_num_0++;
	}
	void Op_pp_1(){
		Op_num_1++;
	}
	void Op_pp(){
		Op_detect_num++;
	}
	void SetEdep(G4double edep, Int_t n, int copyno){
		Edep[n] = edep;
		Scint_CopyN[n] = copyno;
	}
	void SetP(int a){
		Anni_Pos_P = a;
	}
	void SetScint_N(Int_t n){Scint_n=n;}
	int E_judge();
	void Reset();
};

#endif
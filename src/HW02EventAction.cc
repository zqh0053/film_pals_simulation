//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: HW02EventAction.cc 69899 2013-05-17 10:05:33Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 #include "G4SystemOfUnits.hh"
#include "HW02EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "RunAction.hh"
#include "TGraph.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "RootOut.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 using namespace std;
HW02EventAction::HW02EventAction(RunAction* runact)
:frunaction(runact)
{
	posenergy = 0;
	starttrackerScintID = -1;	
	stoptrackerScintID = -2;
	antitrackerScintID = -2;
	minitrackerID = -1;
	Photon_num=0;
	Photon_numtotal=0;
	sture=0;
	sfalse=0;
	eff=0;
	fPriPhotonCerenkov=0;
	fPriPhotonScintillation=0;
	process=0;
	pos[3]={0};
	stopenergy=0;
	antienergy=0;
	startenergy=0;
	minienergy[256]={0};
	starttrigger=0;
	stoptrigger=0;  
	event_num = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
HW02EventAction::~HW02EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void HW02EventAction::BeginOfEventAction(const G4Event*)
{
	RootOut* aRoot = RootOut::GetRootInstance();
	aRoot->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void HW02EventAction::EndOfEventAction(const G4Event* evt)
{
	RootOut* aRoot = RootOut::GetRootInstance();
	int a = aRoot->E_judge();
	event_num++;
	if(event_num%10000 == 0) cout<<event_num<<endl;
	//cout<<a<<endl;
	if(a==1) aRoot->FillEvent();
}

void HW02EventAction::HW02const()
{
	

	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

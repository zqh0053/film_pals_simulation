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
// $Id: HW02SteppingAction.cc 73915 2013-09-17 07:32:26Z gcosmo $
//
/// \file optical/HW02/src/HW02SteppingAction.cc
/// \brief Implementation of the HW02SteppingAction class
//
//
#include "HW02SteppingAction.hh"
#include "HW02EventAction.hh"

 
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "RootOut.hh"
#include <fstream>
#include <cmath>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HW02SteppingAction::HW02SteppingAction(HW02EventAction* evtAct)
:fEventAction(evtAct)
{
   
  fExpectedNextStatus = Undefined;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HW02SteppingAction::~HW02SteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HW02SteppingAction::UserSteppingAction(const G4Step * theStep){
	RootOut* aRoot = RootOut::GetRootInstance();
	G4Track* theTrack = theStep->GetTrack();
	if ( theTrack->GetCurrentStepNumber() == 1 ) fExpectedNextStatus = Undefined;


	G4int parentID = theTrack->GetParentID();
	G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
	G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
	G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
	if(!thePostPV){//out of world
		fExpectedNextStatus=Undefined;
		return;
	}
 //获取湮没位置
	if(theTrack->GetDefinition() == G4Positron::PositronDefinition()&&parentID==1){
		const G4VProcess* process = thePostPoint->GetProcessDefinedStep(); 
		G4String processName = process->GetProcessName();
		G4ThreeVector postion = thePostPoint->GetPosition();
		/*if(processName=="Scintillation" && thePostPV->GetName() == "physample" ){
			fEventAction->sture = 1;  
		}

		if(processName=="Scintillation" && thePostPV->GetName() == "phyanti" ){
			fEventAction->sfalse = 1;		  
		}*/
		if(thePostPV->GetName() == "psource"){
			//aRoot->Set_PosP(0);
			aRoot->SetP(0);
		}
			
		if(thePostPV->GetName() == "world"){
			//aRoot->Set_PosP(1);
			aRoot->SetP(1);
		}
			
		if(thePostPV->GetName() == "pscint0"){
			aRoot->Set_PosP(0,1);
			aRoot->SetP(2);
		}
		if(thePostPV->GetName() == "pscint1"){
			aRoot->Set_PosP(1,1);
			aRoot->SetP(3);
		}
		G4double ke = thePostPoint->GetKineticEnergy()/MeV;
		//aRoot->Set_Kin_E(ke);
		//aRoot->SetPos(postion);
		//aRoot->kn_pp();
	}
  
	G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
	//find the boundary process only once
	G4OpBoundaryProcessStatus boundaryStatus;
    static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;
	if(!boundary)
	{
		G4ProcessManager* pm= theStep->GetTrack()->GetDefinition()->GetProcessManager();
		G4int nprocesses = pm->GetProcessListLength();
		G4ProcessVector* pv = pm->GetProcessList();
		G4int i;
		for( i=0;i<nprocesses;i++)
		{
			if((*pv)[i]->GetProcessName()=="OpBoundary")
			{
				boundary = (G4OpBoundaryProcess*)(*pv)[i];
				break;
			}
		}
	}
	
	G4ParticleDefinition* particleType = theTrack->GetDefinition();
	if(particleType==G4OpticalPhoton::OpticalPhotonDefinition())
	{

		boundaryStatus=boundary->GetStatus();
		if(thePostPoint->GetStepStatus()==fGeomBoundary)
		{

		  // fExpectedNextStatus=Undefined;
		  
			if(boundaryStatus == Detection)
			{  

				if(thePostPV->GetName()=="p_pc" )
				{ 
					aRoot->Op_pp();
	    
 

				}          
				  
			}       
		}
	}

}

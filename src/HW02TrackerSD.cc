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
// $Id: ExN02TrackerSD.cc 69899 2013-05-17 10:05:33Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HW02TrackerSD.hh"
#include "HW02TrackerHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <iostream>
#include <fstream>
#include "G4SystemOfUnits.hh"
#include "RootOut.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
     using namespace std;
HW02TrackerSD::HW02TrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
  evtID=0;
   event_num=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HW02TrackerSD::~HW02TrackerSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HW02TrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new HW02TrackerHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, trackerCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HW02TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

 // if(edep==0.) return false;

  HW02TrackerHit* newHit = new HW02TrackerHit();
  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetEdep     (edep);
  newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetParticleName (aStep->GetTrack()->GetParticleDefinition()->GetParticleName ());
  
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
  G4int copyNo = thePhysical->GetCopyNo();
  newHit->SetCopyNum (copyNo);
  
  
  trackerCollection->insert( newHit );
  //newHit->Print();
  //newHit->Draw();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 
 
void HW02TrackerSD::EndOfEvent(G4HCofThisEvent*)
{   
	RootOut* aRoot = RootOut::GetRootInstance();
    G4double eff1 = 0;
    evtID++;
	G4double energy[100] = {0.};
	//if(evtID%10000==0)
	//G4cout<<"The No. of event = "<<evtID<<G4endl;
	if (trackerCollection) { 
		G4int NbHits = trackerCollection->entries();
		for (G4int i=0;i<NbHits;i++) {
			HW02TrackerHit* hit=(*trackerCollection)[i];
			if(hit->GetParticleName() != "opticalphoton") energy[hit->GetCopyNum()] += hit->GetEdep();
		}
	}
	for(int i = 0;i<4;i++){
		aRoot->SetEdep(energy[i],i,i);
		//cout<<energy[i]<<" ";
	}
	//cout<<endl;
	aRoot->SetScint_N(4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "HW02TrackingAction.hh"
#include "HW02EventAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4EventManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "RootOut.hh"
#include <fstream>
#include <cmath>
using namespace std;

HW02TrackingAction::HW02TrackingAction(HW02EventAction* evtAct)
:fEventAction(evtAct)
{}

HW02TrackingAction::~HW02TrackingAction() {}


void HW02TrackingAction::PreUserTrackingAction(const G4Track* aTrack){
	G4int trackID= aTrack->GetTrackID();
	RootOut* aRoot = RootOut::GetRootInstance();

	G4int parentID = aTrack->GetParentID();
	G4ThreeVector postion = aTrack->GetPosition();
	G4VPhysicalVolume* thePostPV = aTrack->GetVolume();
	if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
		if(thePostPV->GetName() == "pscint0") aRoot->Op_pp_0();
		if(thePostPV->GetName() == "pscint1") aRoot->Op_pp_1();
	} 
}

void HW02TrackingAction::PostUserTrackingAction(const G4Track* aTrack){
	G4int trackID= aTrack->GetTrackID();
	RootOut* aRoot = RootOut::GetRootInstance();

	G4int parentID = aTrack->GetParentID();
	G4ThreeVector postion = aTrack->GetPosition();
	if(parentID == 1&&aTrack->GetDefinition() == G4Positron::PositronDefinition()) aRoot->SetAnni_Pos(postion);
}
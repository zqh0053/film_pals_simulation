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
// $Id: HW02EventAction.hh 69899 2013-05-17 10:05:33Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef HW02EventAction_h
#define HW02EventAction_h 1
#include "G4ios.hh"
#include "G4UserEventAction.hh"
#include "RunAction.hh"
class G4Event;
class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HW02EventAction : public G4UserEventAction
{
  public:
    HW02EventAction(RunAction*);
   ~HW02EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    G4int Photon_num;
	G4int Photon_numtotal;
	G4int sture;
	G4int sfalse;
	G4int process;
	G4int fPriPhotonCerenkov;
	G4int fPriPhotonScintillation;
	G4double pos[3];
	G4double posenergy;
	
  private:
   G4int starttrackerScintID;
   G4int stoptrackerScintID;
   G4int antitrackerScintID;
   G4int minitrackerID;
   G4int starttrigger;
   G4int stoptrigger;
   int event_num;
   	G4double stopenergy;
	G4double antienergy;
	G4double startenergy;	
	G4double minienergy[256];
  G4double eff;
  void HW02const();
  RunAction* frunaction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    

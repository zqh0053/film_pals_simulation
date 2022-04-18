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
// $Id: ExN02TrackerHit.hh 69899 2013-05-17 10:05:33Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HW02TrackerHit_h
#define HW02TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HW02TrackerHit : public G4VHit
{
  public:

      HW02TrackerHit();
     ~HW02TrackerHit();
      HW02TrackerHit(const HW02TrackerHit&);
      const HW02TrackerHit& operator=(const HW02TrackerHit&);
      G4int operator==(const HW02TrackerHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  
      void SetTrackID  (G4int track)      { trackID = track; };
      void SetEdep     (G4double de)      { edep = de; };
      void SetPos      (G4ThreeVector xyz){ pos = xyz; };
	  void SetCopyNum (G4int copyn) { CopyNum = copyn;};
	  void SetParticleName (G4String pname) { p_name = pname;};
      
      G4int GetTrackID()    { return trackID; };
	  G4int GetCopyNum()    { return CopyNum; };
      G4double GetEdep()    { return edep; };      
      G4ThreeVector GetPos(){ return pos; };
	  G4String GetParticleName() { return p_name;};
      
  private:
  
      G4int         trackID;
	  G4int 		CopyNum;
      G4double      edep;
      G4ThreeVector pos;
	  G4String p_name;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<HW02TrackerHit> HW02TrackerHitsCollection;

extern G4Allocator<HW02TrackerHit> HW02TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* HW02TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) HW02TrackerHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void HW02TrackerHit::operator delete(void *aHit)
{
  HW02TrackerHitAllocator.FreeSingle((HW02TrackerHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

#ifndef HW02NovicePrimaryGeneratorAction_h
#define HW02NovicePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
//class OpNovicePrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HW02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    HW02PrimaryGeneratorAction();
    virtual ~HW02PrimaryGeneratorAction();

  public:
     void GeneratePrimaries(G4Event*);



  private:
    G4ParticleGun* particleGun;
	G4ParticleGun* fParticleGun;
    //OpNovicePrimaryGeneratorMessenger* fGunMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*HW02NovicePrimaryGeneratorAction_h*/

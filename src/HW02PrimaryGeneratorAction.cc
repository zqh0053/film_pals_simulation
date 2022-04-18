#include "HW02PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include <cmath>
#include <fstream>

/*HW02PrimaryGeneratorAction::HW02PrimaryGeneratorAction()
{
	G4int n_particle = 1;
	particleGun = new G4ParticleGun(n_particle);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	//粒子类型
	particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="e+"));
	//能量
	particleGun->SetParticleEnergy(10*keV);
	//动量方向
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}


HW02PrimaryGeneratorAction::~HW02PrimaryGeneratorAction()
{
  delete particleGun;
}

void HW02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  
	particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,-10*cm));

	particleGun->GeneratePrimaryVertex(anEvent);
}*/
HW02PrimaryGeneratorAction::HW02PrimaryGeneratorAction()
{
	G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
  
    fParticleGun->SetParticleEnergy(0*eV);
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}

HW02PrimaryGeneratorAction::~HW02PrimaryGeneratorAction()
{
  delete fParticleGun;
}
void HW02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4double rndm, rndm2;
	G4double px,py;
	int ra = 0;
	while(ra == 0){
		rndm = G4UniformRand();
		rndm2 = G4UniformRand();
		px = rndm - 0.5;
		py = rndm2 - 0.5;
		if(sqrt(px*px+py*py)<0.5){
			ra = 1;
		}
	}
	
	fParticleGun->SetParticlePosition(G4ThreeVector(px,py,0.));
	
	if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {  
		G4int Z = 11, A = 22;
		G4double ionCharge   = 0.*eplus;
		G4double excitEnergy = 0.*keV;
		
		G4ParticleDefinition* ion
		   = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
		fParticleGun->SetParticleDefinition(ion);
		fParticleGun->SetParticleCharge(ionCharge);
	}  
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

/*HW02PrimaryGeneratorAction::HW02PrimaryGeneratorAction()
{

  
	
	G4int n_particle = 1;
	particleGun = new G4ParticleGun(n_particle);

	
}

HW02PrimaryGeneratorAction::~HW02PrimaryGeneratorAction()
{
	delete particleGun;
}



void HW02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  //G4bool iso=false;
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="e+"));
  particleGun->SetParticleEnergy(50*keV);
 // if(iso)
  //{
	    //isotropic source
  G4double rndm, rndm2;
  G4double px, py, pz;
  G4double Phi;
  G4double MinPhi =0, MaxPhi=twopi;
  G4double MinTheta=0, MaxTheta=pi;
  G4double sintheta, sinphi, costheta, cosphi;
  rndm = G4UniformRand();
  costheta = std::cos(MinTheta) - rndm * (std::cos(MinTheta) - std::cos(MaxTheta));
  sintheta = std::sqrt(1. - costheta*costheta);
  
  rndm2 = G4UniformRand();
  Phi = MinPhi + (MaxPhi - MinPhi) * rndm2; 
  sinphi = std::sin(Phi);
  cosphi = std::cos(Phi);

  px = -sintheta * cosphi;
  py = -sintheta * sinphi;
  pz = -costheta;

  G4double ResMag = std::sqrt((px*px) + (py*py) + (pz*pz));
  px = px/ResMag;
  py = py/ResMag;
  pz = pz/ResMag;
 

   particleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
   particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  particleGun->GeneratePrimaryVertex(anEvent);   
  //}
}*/



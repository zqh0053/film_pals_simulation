//Ì½²âÆ÷Í·ÎÄ¼þ



#ifndef HW02DetectorConstruction_H
#define HW02DetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class HW02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    HW02DetectorConstruction(G4String);
    ~HW02DetectorConstruction();

    G4VPhysicalVolume* Construct();

  private:
    
    // Logical volumes
    //
	G4LogicalVolume* logicWorld;
	G4LogicalVolume* logicScint;
	G4LogicalVolume* logicTeflon;

    // Physical volumes
    //
	G4VPhysicalVolume* physiWorld;
    G4VPhysicalVolume* physiTeflon;
	G4VPhysicalVolume* physiScint;
	G4VPhysicalVolume* sioil1_phys;
	G4VPhysicalVolume* photocathod_phys;
	G4String file_in;
	

};

#endif

#ifndef HW02TrackingAction_H
#define HW02TrackingAction_H 1
#include "HW02TrackingAction.hh"
#include "globals.hh"
#include "G4UserTrackingAction.hh"
#include "HW02EventAction.hh"

class HW02TrackingAction : public G4UserTrackingAction
{
public:
  HW02TrackingAction(HW02EventAction* );
  virtual ~HW02TrackingAction();
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);
private:
    //G4int gammaCounter;
	HW02EventAction* fEventAction;

};

#endif
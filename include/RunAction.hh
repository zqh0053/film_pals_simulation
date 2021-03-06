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
/// \file electromagnetic/TestEm8/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id: RunAction.hh 85243 2014-10-27 08:22:42Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   RunAction
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//
#ifndef RunAction_h
#define RunAction_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4UserRunAction.hh"
#include "globals.hh"

#include "g4root.hh"
#include "RootOut.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
public: 

	RunAction();
	virtual ~RunAction();

	// virtual G4Run* GenerateRun(); 

	virtual void BeginOfRunAction(const G4Run*);
	// In this method histogramms are booked

	virtual void EndOfRunAction(const G4Run*);
	// In this method bookHisto method is called in which histogramms are filled

public:
	G4int eventnumber;
	G4int truenum;
	G4int falsenum;
	G4int totalnum;
	G4int sturenum;
	G4int pnum;
	G4int spnum;
	G4int sfnum;

private:


};

#endif


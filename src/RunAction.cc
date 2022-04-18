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
/// \file electromagnetic/TestEm8/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 85243 2014-10-27 08:22:42Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   RunAction
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"


#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 using namespace std;
RunAction::RunAction()
  : G4UserRunAction()
{
	eventnumber = 0;
	truenum=0;
falsenum=0;
totalnum=0;
sturenum=0;
pnum=0;
spnum=0;
sfnum=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void RunAction::Book()
{
  // Always creating analysis manager
  fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->SetActivation(true);
  fAnalysisManager->SetVerboseLevel(1);

  // Open file histogram file
  fAnalysisManager->OpenFile(fHistName);
  fAnalysisManager->SetFirstHistoId(1);

  TestParameters* param = TestParameters::GetPointer();
  G4int nBinsE = param->GetNumberBins();
  G4int nBinsCluster = param->GetNumberBinsCluster();
  G4double maxEnergy = param->GetMaxEnergy();
  G4double factorALICE = param->GetFactorALICE();

  // Creating an 1-dimensional histograms in the root directory of the tree
  fAnalysisManager->CreateH1("h1","Energy deposition in detector (keV)",
                             nBinsE,0.0,maxEnergy/keV);
  fAnalysisManager->CreateH1("h2","Number of primary clusters",
                             nBinsE,-0.5,G4double(nBinsCluster)-0.5);
  fAnalysisManager->CreateH1("h3","Energy deposition in detector (ADC)",
                             nBinsE,0.0,maxEnergy* factorALICE/keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(); 
  return fRun;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
	RootOut::GetRootInstance()->BeginOfRunAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
	RootOut::GetRootInstance()->EndOfRunAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

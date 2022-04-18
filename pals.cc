#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "HW02DetectorConstruction.hh"
#include "HW02PhysicsList.hh"
#include "HW02PrimaryGeneratorAction.hh"
#include "HW02EventAction.hh"
#include "HW02StackingAction.hh"
#include "HW02SteppingAction.hh"
#include "HW02TrackingAction.hh"
#include "RunAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


#include "G4MPImanager.hh"
#include "G4MPIsession.hh"
#include "G4VMPIseedGenerator.hh"
#include "RootOut.hh"

#include <iostream>

using namespace std;

int main(int argc,char** argv)
{
	//
	G4MPImanager* g4MPI = new G4MPImanager(argc, argv);
	G4MPIsession* session = g4MPI-> GetMPIsession();


 // Set the Random engine
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
	G4int rank_ = MPI::COMM_WORLD.Get_rank(); //获取并行进程号
	time_t systime = time(NULL);
	G4long seed =systime + rank_*10.0;
	G4Random::setTheSeed(seed);
	G4cout<<seed<<"_____________________"<<G4endl;
  
	G4RunManager * runManager = new G4RunManager;
	
	string filename;
	if(access(argv[1],0) == -1){
		cout<<"input is null"<<endl;
		filename = "vis.mac";
	}
	else{
		filename = argv[1];
	}
	
	RootOut* myRootOutput = new RootOut(filename);
	
	//初始化探测器
	G4VUserDetectorConstruction* detector = new HW02DetectorConstruction(filename);
	runManager->SetUserInitialization(detector);
	//初始化物理过程
	G4VUserPhysicsList* physics = new HW02PhysicsList;
	runManager->SetUserInitialization(physics);
	//初始化粒子发生器
	G4VUserPrimaryGeneratorAction* gen_action = new HW02PrimaryGeneratorAction;
	runManager->SetUserAction(gen_action);
	RunAction* runaction = new RunAction();
	runManager->SetUserAction(runaction);	   
	//
	HW02EventAction* evtAction = new HW02EventAction(runaction);
	runManager->SetUserAction(evtAction);
	runManager->SetUserAction(new HW02StackingAction(evtAction));
	runManager->SetUserAction(new HW02SteppingAction(evtAction));
	runManager->SetUserAction(new HW02TrackingAction(evtAction));
    

	runManager->Initialize();
	// 获取指向 UI manager 的指针
    //
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/run/verbose 0");
    UI->ApplyCommand("/event/verbose 0");
    UI->ApplyCommand("/tracking/verbose 0");

    // Start a run
    //
    //G4int numberOfEvent = 100;
    //runManager->BeamOn(numberOfEvent);
  
    //
#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif    
     
    // Get the pointer to the User Interface manager
    //
    G4UImanager * UImanager = G4UImanager::GetUIpointer();  

    if (argc!=1)   // batch mode  
      {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
      }
    else           // interactive mode : define UI session
      { 
#ifdef G4UI_USE
        G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
        ui->SessionStart();
        delete ui;
#endif
     
#ifdef G4VIS_USE
        delete visManager;
#endif     
      }

    // Job termination
    //
    // Free the store: user actions, physics_list and detector_description are
    //                 owned and deleted by the run manager, so they should not
    //                 be deleted in the main() program !
    //
	delete g4MPI;
	delete myRootOutput;
    delete runManager;

    return 0;
}

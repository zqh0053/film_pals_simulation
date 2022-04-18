#include "HW02DetectorConstruction.hh"
#include "HW02TrackerSD.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Colour.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4EllipticalCone.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "math.h"
#include <cstring>
#include <fstream>

using namespace std;

HW02DetectorConstruction::HW02DetectorConstruction(G4String filename)
 :  logicWorld(0),  logicScint(0),  logicTeflon(0),
 physiWorld(0), 
 physiTeflon(0),physiScint(0),sioil1_phys(0),photocathod_phys(0)
{
	file_in = filename;
}

HW02DetectorConstruction::~HW02DetectorConstruction()
{;}

G4VPhysicalVolume* HW02DetectorConstruction::Construct()
{
	G4NistManager* man = G4NistManager::Instance();
	G4double a, z;
	G4double density, temperature, pressure,fractionmass;
	G4int nel,ncomponents;

	//Air
	G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
	G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);

	G4Material* Air = man->FindOrBuildMaterial("G4_AIR");
	//Air->AddElement(N, 70*perCent);
	//Air->AddElement(O, 30*perCent);
	
	//vacuum
	G4Material* Vacuum = new G4Material("Vacuum", density= 0.001*mg/cm3, nel=2);
	Vacuum->AddElement(N, 70*perCent);
	Vacuum->AddElement(O, 30*perCent);

	//Lead
	G4Material* Pb = 
	new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

	//Fe
	G4Material* Fe = new G4Material("Iron",z=26.,a=55.845*g/mole,density=7.86*g/cm3);
	 
	//Teflon
	G4Element* elC=new G4Element("Carbon","Ba",z=6., 12.011*g/mole);
	G4Element* elH=new G4Element("Hydrogen" , "H", z=1. , 1.008*g/mole );
	G4Element* elF = new G4Element("Flourine", "F", z = 9.0, 18.998*g/mole);
	G4Material* Teflon = new G4Material("Teflon", density = 2.20*g/cm3, 3);
	Teflon->AddElement(elC, 2);
	Teflon->AddElement(elH, 2);
	Teflon->AddElement(elF, 2);

	//BaF2
	//G4Material* Barium    = man->FindOrBuildMaterial("G4_Ba");
	G4Element* elBa = new G4Element("Barium","Ba",z=56.,137.372*g/mole);
	G4Material* BaF2 = new G4Material("BaF2",density = 4.893*g/cm3, 2);
	BaF2->AddElement(elBa, 1);
	BaF2->AddElement(elF, 2);
	
	//kapton
	G4Material* kapton = new G4Material("kapton",density = 1.401*g/cm3, 4);
	kapton->AddElement(N, 2);
	kapton->AddElement(elH, 10);
	kapton->AddElement(elC, 22);
	kapton->AddElement(O, 5);


	//ej-200 发光
	G4Material* ej_faguang = new G4Material("EJ-faguang", density=1.023*g/cm3,2);
	ej_faguang -> AddElement(elH, 10);
	ej_faguang -> AddElement(elC, 11);

	 //ej-200 不发光
	G4Material* ej_wuguang = new G4Material("EJ-wuguang", density=1.023*g/cm3,2);
	ej_wuguang -> AddElement(elH, 10);
	ej_wuguang -> AddElement(elC, 11);
  
  	//silicon grease BC630 CH3-SiO-C6H5   https://github.com/rat-pac/rat-pac/blob/master/data/MATERIALS.ratdb
	G4Element* elO = man->FindOrBuildElement("O");
	G4Element* elSi = man->FindOrBuildElement("Si");
	G4Material* Grease = new G4Material("Grease", density = 1.06*g/cm3,4);
	Grease->AddElement(elSi,fractionmass=0.2060);
	Grease->AddElement(elC,fractionmass=0.6173);
	Grease->AddElement(elH,fractionmass=0.0592);
	Grease->AddElement(elO,fractionmass=0.1175);
	
    //Bialkali photocathode
	G4Element* K = man->FindOrBuildElement("K");  
    G4Element* Cs = man->FindOrBuildElement("Cs");
    G4Element* Sb = man->FindOrBuildElement("Sb");
	G4Material* Bialkali = new G4Material("Bialkali", density= 4.28*g/cm3, ncomponents=3);
	Bialkali->AddElement(K, fractionmass=0.133);
	Bialkali->AddElement(Cs, fractionmass=0.452);
	Bialkali->AddElement(Sb, fractionmass=0.415); 
	
	//LaBr3
	G4Element* elLa = man->FindOrBuildElement("La");
	G4Element* elBr = man->FindOrBuildElement("Br");
	G4Element* elCe = man->FindOrBuildElement("Ce");
	G4Material* LaBr_Ce = new G4Material("labr", density = 5.08*g/cm3, 3);
	LaBr_Ce->AddElement(elLa, fractionmass=34.85*perCent);
	LaBr_Ce->AddElement(elBr, fractionmass=5.*perCent);
	LaBr_Ce->AddElement(elCe, fractionmass=60.15*perCent);	
	
	//LYSO
	G4Element* elLu = man->FindOrBuildElement("Lu");
	G4Element* elY = man->FindOrBuildElement("Y");
	G4Material* LYSO = new G4Material("lyso", density = 7.15*g/cm3, 5);
	LYSO->AddElement(elLu, fractionmass=24.9375*perCent);
	LYSO->AddElement(elY, fractionmass=13.0625*perCent);
	LYSO->AddElement(elSi, fractionmass=9.5*perCent);
	LYSO->AddElement(elO, fractionmass=47.5*perCent);
	LYSO->AddElement(elCe, fractionmass=5*perCent);
	
	
	//pmt glass
	G4Material* Oxigen      = man->FindOrBuildMaterial("G4_O");	
	G4Material* Silicon     = man->FindOrBuildMaterial("G4_Si");	
	G4Material* Potassium   = man->FindOrBuildMaterial("G4_K");	
	G4Material* Sodium      = man->FindOrBuildMaterial("G4_Na");	
	G4Material* lead        = man->FindOrBuildMaterial("G4_Pb");	
	G4Material* Boron       = man->FindOrBuildMaterial("G4_B");	
	G4Material* Aluminum    = man->FindOrBuildMaterial("G4_Al");
	G4Material* Glass = new G4Material("Borosicilate",density=2.23*g/cm3,6);
	Glass->AddMaterial(Oxigen,     0.539562); 
	Glass->AddMaterial(Boron,      0.040064);
	Glass->AddMaterial(Sodium,     0.028191);
	Glass->AddMaterial(Aluminum,   0.011644);
	Glass->AddMaterial(Silicon,    0.377220);
	Glass->AddMaterial(Potassium,  0.003321);	
	
	G4Material* NaI = man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
	
	//314
	G4Element* elFe = man->FindOrBuildElement("Fe");
	G4Element* elCr = man->FindOrBuildElement("Cr");
	G4Element* elNi = man->FindOrBuildElement("Ni");
	G4Element* elMn = man->FindOrBuildElement("Mn");
	G4Element* elS = man->FindOrBuildElement("S");
	G4Element* elP = man->FindOrBuildElement("P");
	G4Material* steel314 = new G4Material("steel-314",density=7.85*g/cm3,8);
	steel314->AddElement(elC, fractionmass=0.25*perCent);
	steel314->AddElement(elMn, fractionmass=2.0*perCent);
	steel314->AddElement(elP, fractionmass=0.045*perCent);
	steel314->AddElement(elS, fractionmass=0.03*perCent);
	steel314->AddElement(elSi, fractionmass=2.0*perCent);
	steel314->AddElement(elCr, fractionmass=25.0*perCent);
	steel314->AddElement(elNi, fractionmass=20.0*perCent);
	steel314->AddElement(elFe, fractionmass=50.675*perCent);
	
	const G4int num = 2;
    G4double Ephoton[num]={2.5*eV, 3.1*eV};
	
	//ej-faguang
	const G4int ej_num = 20;
	G4double ej_energy[ej_num] = {2.5*eV,2.53158*eV,2.56316*eV,2.59474*eV,2.62632*eV,2.65789*eV,2.68947*eV,2.72105*eV,2.75263*eV,
									2.78421*eV,2.81579*eV,2.84737*eV,2.87895*eV,2.91053*eV,2.94211*eV,2.97368*eV,3.00526*eV,3.03684*eV,3.06842*eV,3.1*eV};
	G4double ej_rindex[ej_num] = {1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58};
	G4double ej_absorption[ej_num] = {70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,
                                     	70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm,70.*cm};
	G4double ej_fastcompnent[ej_num] = {0.06861,0.10003,0.13918,0.19388,0.26824,0.35968,0.42221,0.47399,0.54315,0.63267,0.74095,0.85994,0.94624,0.99999,
                                        0.93294,0.66705,0.38466,0.19654,0.05713,0.00961};
	
	G4MaterialPropertiesTable* ejPMT1 = new G4MaterialPropertiesTable();
	ejPMT1->AddProperty("RINDEX", ej_energy, ej_rindex,ej_num);
	ejPMT1->AddProperty("ABSLENGTH",ej_energy, ej_absorption, ej_num);
	ejPMT1->AddProperty("FASTCOMPONENT",ej_energy, ej_fastcompnent, ej_num);	
    ejPMT1->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);  
	ejPMT1->AddConstProperty("RESOLUTIONSCALE",1.0);
	ejPMT1->AddConstProperty("FASTTIMECONSTANT",2.1*ns);  
	ejPMT1->AddConstProperty("FASTSCINTILLATIONRISETIME",0.9*ns);  
	ejPMT1->AddConstProperty("YIELDRATIO",1.0); 
	ej_faguang->SetMaterialPropertiesTable(ejPMT1);	
	ej_faguang->GetIonisation()->SetBirksConstant(0.133*mm/MeV); //from PPT: Photon Transport Simulations
	
	//Silicon Grease	1.5
	G4double greaseRefractiveIndex[num] = { 1.5,1.5};
	G4double greaseAbsorption [ num ] =	{ 300*mm, 300*mm};

	G4MaterialPropertiesTable* GreaseMPT = new G4MaterialPropertiesTable();
	GreaseMPT->AddProperty("RINDEX", Ephoton, greaseRefractiveIndex, num);
	GreaseMPT->AddProperty("ABSLENGTH",    Ephoton, greaseAbsorption, num);  
	Grease->SetMaterialPropertiesTable(GreaseMPT);
	
	G4double AirRefractiveIndex[num] = { 1.00029,1.00029};
	G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
	AirMPT->AddProperty("RINDEX", Ephoton, AirRefractiveIndex, num);
	AirMPT->AddProperty("ABSLENGTH",    Ephoton, greaseAbsorption, num);  
	Air->SetMaterialPropertiesTable(AirMPT);
	
	G4bool overloops = true;
	//世界体定义
	G4double worldxy=50*cm;
    G4double worldz=50*cm;
  

    G4Box* solidWorld= new G4Box("world",worldxy/2,worldxy/2,worldz/2);
    logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);
  
    physiWorld = new G4PVPlacement(0,               // no rotation
                                   G4ThreeVector(), // at (0,0,0)
                                   logicWorld,      // its logical volume
                                   "World",         // its name
                                   0,               // its mother  volume
                                   false,           // no boolean operations
                                   0,               // copy number
								   overloops);              
                                 
                                
  //------------------------------ 
  // Scint
  //------------------------------  
	//scint
	fstream filein(file_in,ios::in);
	string line;
	G4double r_0,deg_0;  //探测器角度与位置
	G4double x_0,y_0,z_0;  //探测器中心偏移
	x_0 = y_0 = z_0 = 0;
	G4double sizex,sizey,sizez;
	G4double reflect; //反射率
	reflect = 0.99;
	sizex = sizey = sizez = 30;
	G4double sample_thick, source_thick, air_thick;
	sample_thick = 1.5;
	source_thick = 0.0075;
	air_thick = 0.0001;
	G4Material* sample_m;
	while(getline(filein,line)){
		istringstream sin(line);
		string command;
		G4String el_name;
		int el_n;
		sin>>command;
		if(command == "#/scint/pos_offset"){
			sin>>x_0>>y_0>>z_0;
		}
		if(command == "#/scint/pos"){
			sin>>r_0>>deg_0;
		}
		if(command == "#/scint/size"){
			sin>>sizex>>sizez;
		}
		if(command == "#/material"){
			string m_name;
			double density_m;
			int el_num;
			sin>>m_name>>density_m>>el_num;
			sample_m = new G4Material(m_name, density = density_m*1*g/cm3, el_num);
		}
		if(command == "#/element"){
			sin>>el_name>>el_n;
			G4Element* el_X = man->FindOrBuildElement(el_name);
			sample_m->AddElement(el_X, el_n);
		}
		if(command == "#/material_name"){
			string m_name;
			sin>>m_name;
			sample_m = man->FindOrBuildMaterial(m_name);	
		}	
		if(command == "#/thick"){
			sin>>sample_thick>>source_thick>>air_thick;
		}
		if(command == "#/reflect"){
			sin>>reflect;
		}
	}
	
	
	//塑闪
	G4Tubs* scint_ps = new G4Tubs("scint_tub",0.,26.*mm/2,1.*mm/2,0.*deg,360.*deg);
	G4LogicalVolume* logicscint = new G4LogicalVolume(scint_ps, ej_faguang ,"logicscint",0,0,0);
	
	//teflon
	G4Tubs* scint_tfl = new G4Tubs("scint_tub_tfl",0.,(26.*mm+0.02*mm)/2,(1.*mm+0.01*mm)/2,0.*deg,360.*deg);
	G4LogicalVolume* logicscint_tfl = new G4LogicalVolume(scint_tfl, Teflon ,"logicscint",0,0,0);
	G4PVPlacement* phy_tfl = new G4PVPlacement(0,G4ThreeVector(0,0,source_thick/2 + air_thick + 0.5*mm + 0.005*mm),logicscint_tfl,"ptfl",logicWorld,false,0,overloops);
	
	//scint_0
	G4PVPlacement* phy_scint0 = new G4PVPlacement(0,G4ThreeVector(0, 0, -0.005*mm),logicscint,"pscint0",logicscint_tfl,false,0,overloops);
	//scint_1
	G4PVPlacement* phy_scint1 = new G4PVPlacement(0,G4ThreeVector(0, 0, -source_thick/2 - air_thick - 0.5*mm),logicscint,"pscint1",logicWorld,false,1,overloops);
	
	//硅油
	G4Tubs* grease_t = new G4Tubs("scint_tub",0.,26.*mm/2,0.1*mm/2,0.*deg,360.*deg);
	G4LogicalVolume* logic_grease = new G4LogicalVolume(grease_t, Grease ,"logic_grease",0,0,0);
	G4PVPlacement* phy_grease = new G4PVPlacement(0,G4ThreeVector(0, 0, -source_thick/2 - air_thick - 1.*mm - 0.1*mm/2),logic_grease,"pgrease",logicWorld,false,0,overloops);
	
	//光阴极
	G4Tubs* photocathode = new G4Tubs("photocathode _tub",0.,28.*mm/2,1.*mm/2,0.*deg,360.*deg);
	G4LogicalVolume* logic_pc = new G4LogicalVolume(photocathode, Bialkali ,"logic_pc",0,0,0);
	G4PVPlacement* phy_pc = new G4PVPlacement(0,G4ThreeVector(0, 0, -source_thick/2 - air_thick - 1.*mm - 0.1*mm - 0.5*mm),logic_pc,"p_pc",logicWorld,false,0,overloops);
  
	//样品和源
	G4Box* solidsource = new G4Box("source",0.5*cm,0.5*cm, 0.5*source_thick);
	G4LogicalVolume* logicsource = new G4LogicalVolume(solidsource, kapton ,"logicsource",0,0,0);
	G4PVPlacement* phy_source = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicsource,"psource",logicWorld,false,0,overloops);
  
	G4double distance = 0.1*cm;
	//
	G4double Rgamma = 3*cm;
	//通用
	G4double tef_thick=0.25*mm;
	
	
	
	//sample
	/*G4Box* solidsample = new G4Box("sample",0.5*cm,0.5*cm,sample_thick);
	G4LogicalVolume* logicsample = new G4LogicalVolume(solidsample, sample_m ,"logicsample",0,0,0);
	G4PVPlacement* phy_sample = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicsample,"psample",logicWorld,false,0,overloops);
	//air
	G4Box* solidair = new G4Box("sample",0.5*cm,0.5*cm,source_thick + air_thick);
	G4LogicalVolume* logicair = new G4LogicalVolume(solidair, Air ,"logicair",0,0,0);
	G4PVPlacement* phy_air = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicair,"pair",logicsample,false,0,overloops);
	//source
	G4Box* solidsource = new G4Box("source",0.5*cm,0.5*cm,source_thick);
	G4LogicalVolume* logicsource = new G4LogicalVolume(solidsource, kapton ,"logicsource",0,0,0);
	G4PVPlacement* phy_source = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicsource,"psource",logicair,false,0,overloops);*/
	
	
	deg_0 = deg_0*1.*deg;
	G4RotationMatrix* rm1= new G4RotationMatrix();
	rm1->rotateX(deg_0);
	G4RotationMatrix* rm2= new G4RotationMatrix();
	rm2->rotateX(-deg_0);
	
	
	G4Tubs* scint_1 = new G4Tubs("scint_tub",0.,sizex/2,sizez/2,0.*deg,360.*deg);
	G4LogicalVolume* logicscint2 = new G4LogicalVolume(scint_1, LaBr_Ce ,"logicscint2",0,0,0);
	G4PVPlacement* phy_pscint0 = new G4PVPlacement(rm1,G4ThreeVector(x_0,-y_0+r_0*cos(90.*deg - deg_0),-z_0+r_0*sin(90.*deg - deg_0)),logicscint2,"ppscint0",logicWorld,false,2,overloops);
	G4PVPlacement* phy_pscint1 = new G4PVPlacement(rm2,G4ThreeVector(x_0,y_0-r_0*cos(90.*deg - deg_0),z_0+r_0*sin(90.*deg - deg_0)),logicscint2,"ppscint1",logicWorld,false,3,overloops);
	G4cout<<G4ThreeVector(x_0,-r_0*cos(90.*deg - deg_0),-r_0*sin(90.*deg - deg_0))<<endl;
	
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String trackerSDname = "TrackerSD";
    HW02TrackerSD* aTrackerSD = new HW02TrackerSD( trackerSDname );
    SDman->AddNewDetector(aTrackerSD);
	logicscint->SetSensitiveDetector(aTrackerSD);
	logicscint2->SetSensitiveDetector(aTrackerSD);
   
    //color
	G4VisAttributes* coverVis3= new G4VisAttributes(G4Colour(0.,0.9,0.1,0.8));
	G4VisAttributes* coverVis1= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.7));
	G4VisAttributes* coverVis2= new G4VisAttributes(G4Colour(0.0,0.1,0.9,0.5));
	G4VisAttributes* coverVis4= new G4VisAttributes(G4Colour(1,1,1,1));
	
	coverVis3->SetForceSolid(true);
	//coverVis3->SetForceWireframe(true);
	coverVis1->SetForceSolid(true);
	coverVis2->SetForceSolid(true);
	coverVis4->SetForceSolid(true);
	
	logicscint->SetVisAttributes(coverVis3);
	logic_pc->SetVisAttributes(coverVis2);
	logicsource->SetVisAttributes(coverVis1);
	//logicWorld->SetVisAttributes(coverVis4);
	
	//光学界面
	G4OpticalSurface* OpTeflonSurface = new G4OpticalSurface("TeflonSurface");
	OpTeflonSurface->SetModel(unified);	//if ground, unified model should be chosen.
	OpTeflonSurface->SetType(dielectric_dielectric);
	OpTeflonSurface->SetFinish(groundfrontpainted);

	G4MaterialPropertiesTable *Teflonsurface = new G4MaterialPropertiesTable();
	G4double TeflonsurfaceReflectivity[2] = {reflect,reflect};	 
	G4double TeflonsurfaceEfficiency[2] = {0.,0.};	 
	Teflonsurface->AddProperty("REFLECTIVITY", Ephoton, TeflonsurfaceReflectivity, 2);
	Teflonsurface->AddProperty("EFFICIENCY", Ephoton, TeflonsurfaceEfficiency, 2);
	OpTeflonSurface->SetMaterialPropertiesTable(Teflonsurface);
	
	//MPPC
	/*G4OpticalSurface* OpCathodSurface = new G4OpticalSurface("CathodSurface");
	OpCathodSurface->SetModel(unified);	//if ground, unified model should be chosen.
	OpCathodSurface->SetType(dielectric_metal);
	OpCathodSurface->SetFinish(polished);
	G4MaterialPropertiesTable *CathodSurface = new G4MaterialPropertiesTable();
	const G4int Cathod_num = 20;
	
	G4double CathodEnergy[Cathod_num]={1.5834*eV, 1.6800*eV, 1.7928*eV, 1.9438*eV, 2.0812*eV, 
									2.2395*eV, 2.3700*eV, 2.4802*eV, 2.5775*eV, 2.7083*eV, 
									2.8917*eV, 3.0036*eV, 3.2066*eV, 3.3188*eV, 3.5388*eV, 
									3.7734*eV, 4.0222*eV, 4.2632*eV, 4.4872*eV, 4.8435*eV};										
	G4double CathodSurfaceReflectivity[Cathod_num] = {  0.,0.,0,0,0,
														0,0,0,0,0,
														0,0,0,0,0,
														0,0,0,0,0
														};	
	G4double CathodSurfaceEfficiency[Cathod_num]={0.0552, 0.0822, 0.1166, 0.1669, 0.2148, 
												0.2774, 0.3363, 0.3928, 0.4358, 0.4775, 
												0.5107, 0.5021, 0.4542, 0.4088, 0.3756, 
												0.3400, 0.2786, 0.1853, 0.1104, 0.0564 };											

  												
	CathodSurface->AddProperty("REFLECTIVITY", CathodEnergy, CathodSurfaceReflectivity, Cathod_num);
	CathodSurface->AddProperty("EFFICIENCY",CathodEnergy,CathodSurfaceEfficiency,Cathod_num);
	OpCathodSurface->SetMaterialPropertiesTable(CathodSurface); */
	
	G4OpticalSurface* OpCathodSurface = new G4OpticalSurface("CathodSurface");
		OpCathodSurface->SetModel(unified);	//if ground, unified model should be chosen.
		OpCathodSurface->SetType(dielectric_metal);
		OpCathodSurface->SetFinish(polished);
	G4MaterialPropertiesTable *CathodSurface = new G4MaterialPropertiesTable();
		const G4int Cathod_num = 47; 
		G4double CathodEnergy[Cathod_num]={ 1.7534*eV,1.7608*eV,1.7721*eV,1.7797*eV,1.7913*eV,1.7990*eV,1.8108*eV,1.8188*eV,1.8308*eV,1.8430*eV,1.8554*eV,1.8638*eV,1.8764*eV,1.8893*eV,1.9023*eV,1.9199*eV,1.9333*eV,1.9515*eV,1.9701*eV,1.9938*eV,2.0180*eV,2.0429*eV,2.0788*eV,2.1159*eV,2.1600*eV,2.2060*eV,2.2602*eV,2.2978*eV,2.3433*eV,2.3906*eV,2.4837*eV,2.5844*eV,2.6935*eV,2.8316*eV,2.9845*eV,3.1791*eV,3.4009*eV,3.6398*eV,3.8965*eV,4.2134*eV,4.5363*eV,4.9716*eV,5.4632*eV,6.0186*eV,6.7542*eV,7.4863*eV,7.6948*eV};										
		G4double CathodSurfaceReflectivity[Cathod_num] = {  0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000};	
															
		G4double CathodSurfaceEfficiency[Cathod_num]={  0.0001,0.0002,0.0002,0.0003,0.0003,0.0004,0.0005,0.0007,0.0009,0.0011,0.0015,0.0019,0.0024,0.0030,0.0039,0.0049,0.0064,0.0081,0.0102,0.0130,0.0166,0.0209,0.0274,0.0350,0.0426,0.0535,0.0661,0.0843,0.1043,0.1309,0.1547,0.1773,0.2095,0.2295,0.2514,0.2514,0.2514,0.2591,0.2514,0.2295,0.2033,0.1884,0.1747,0.2002,0.2095,0.1773,0.1391 };														
															
		CathodSurface->AddProperty("REFLECTIVITY", CathodEnergy, CathodSurfaceReflectivity, Cathod_num);
		CathodSurface->AddProperty("EFFICIENCY",CathodEnergy,CathodSurfaceEfficiency,Cathod_num);
		OpCathodSurface->SetMaterialPropertiesTable(CathodSurface); 
		
		
	new G4LogicalBorderSurface("ScintTeflonSurface_log", phy_scint0, phy_tfl,OpTeflonSurface);
	new G4LogicalBorderSurface("ScintTeflonSurface_log", phy_grease, phy_pc,OpCathodSurface);
	new G4LogicalBorderSurface("ScintTeflonSurface_log", physiWorld, phy_pc,OpCathodSurface);
	
    return physiWorld;
}

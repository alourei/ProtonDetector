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
// $Id: ProtonConstruction.hh,v 0.1 2013/11/27 17:47:13 D. Perez-Loureiro Exp $
// GEANT4 tag $Name: geant4-09-06p02 $
//
//
// --------------------------------------------------------------
// Based on     GEANT 4 - AstroBox2 and ActarSim2013
// --------------------------------------------------------------
//
//

#ifndef PROTONDETECTORCONSTRUCTION_HH_
#define PROTONDETECTORCONSTRUCTION_HH_

#include <G4VUserDetectorConstruction.hh>

#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class G4Box;
class G4Tubs;
class G4Sphere;
class G4Material;

class G4UnionSolid;
class G4SubtractionSolid;
class G4IntersectionSolid;
class ProtonDetectorConstructionMessenger;
class ProtonDetectorSD;
class IonChamberSD;
class DegraderSD;

class ProtonDetectorConstruction: public G4VUserDetectorConstruction {

private:

	//Volumes
	G4Box *worldVolume;

	//Logical Volumes
	G4LogicalVolume *world_log;
	G4LogicalVolume *gasChamber_log;
	G4LogicalVolume *gasVolume_log;
	G4LogicalVolume *chamberWindow_log;
	G4LogicalVolume *AlDegrader_log;

	//Physical Volumes
	G4VPhysicalVolume *world_phys;
	G4VPhysicalVolume *gasChamber_phys;
	G4VPhysicalVolume *gasVolume_phys;
	G4VPhysicalVolume *chamberWindow_phys;
	G4VPhysicalVolume *AlDegrader_phys;

    //Select the type of volume
	G4String detectorGeometry;
	G4String degraderIncludedFlag;
	void DefineMaterials();
	G4VPhysicalVolume *ConstructDetector();

	//Box Parameters
	G4double XboxLength;
	G4double YboxLength;
	G4double ZboxLength;

	//Tube Parameters
	G4double radiusGasTube;
	G4double lengthGasTube;

	//Degrader Parameters
	G4ThreeVector degraderPosition;
	G4double degraderThickness;
	G4double degraderAngle;

	//The messenger
	ProtonDetectorConstructionMessenger *theMessenger;

	//Sensitive Detectors
	ProtonDetectorSD *theProtonDetectorSD;
	IonChamberSD *theIonChamberSD ;
	DegraderSD *theDegraderSD;

public:
	ProtonDetectorConstruction();
	~ProtonDetectorConstruction();
	G4VPhysicalVolume  *Construct();
	void SetXGasBox(G4double val){XboxLength = val;}
	void SetYGasBox(G4double val){YboxLength = val;}
	void SetZGasBox(G4double val){ZboxLength = val;}
	void SetRadiusGasTub(G4double val){radiusGasTube = val;}
	void SetLengthGasTub(G4double val){lengthGasTube = val;}
	void SetDegraderPosition(G4ThreeVector pos){degraderPosition = pos;}
	void SetDegraderThickness(G4double val){degraderThickness = val;}
	void SetDegraderAngle(G4double val){degraderAngle = val;}

	void SetDetectorGeometry(G4String geometry){detectorGeometry = geometry;}
	void SetDegraderIncludedFlag(G4String value){degraderIncludedFlag = value;}

	G4String GetDetectorGeometry(){return detectorGeometry;}
	G4String GetDegraderIncludedFlag(){return degraderIncludedFlag;}
	G4double GetXGasBox(void){return XboxLength;}
	G4double GetYGasBox(void){return YboxLength;}
	G4double GetZGasBox(void){return ZboxLength;}
	G4double GetRadiusGasTub(void){return radiusGasTube;}
	G4double GetLengthGasTub(void){return lengthGasTube;}
	G4ThreeVector GetDegraderPosition() {return degraderPosition;}
	G4double GetDegradetThickness() {return degraderThickness;}
	G4double GetDegraderAngle(){return degraderAngle;}

	void UpdateGeometry();
	void PrintDetectorParameters();

};

#endif /* PROTONDETECTORCONSTRUCTION_HH_ */

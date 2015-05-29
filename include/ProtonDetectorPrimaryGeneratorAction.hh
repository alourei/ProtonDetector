/*
 * ProtonDetectorPrimaryGeneratorAction.hh
 *
 *  Created on: Dec 4, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORPRIMARYGENERATORACTION_HH_
#define PROTONDETECTORPRIMARYGENERATORACTION_HH_

#include <globals.hh>
#include <G4LorentzVector.hh>
#include "MargotDataRecordTree.hh"
#include <G4VUserPrimaryGeneratorAction.hh>

#include <TFile.h>
#include <TH1F.h>


class G4ParticleGun;
class G4Event;
class G4GeneralParticleSource;
class G4SingleParticleSource;
class G4SPSEneDistribution;
class G4SPSPosDistribution;
class G4SPSAngDistribution;

class TFile;
class TH1F;



class ProtonDetectorPrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction {
public:
	ProtonDetectorPrimaryGeneratorAction();
	ProtonDetectorPrimaryGeneratorAction(G4int numSources,G4String pName, G4String sourceType, G4double beamEnergy );
	~ProtonDetectorPrimaryGeneratorAction();
	void GeneratePrimaries(G4Event *anEvent);

private:

	G4ParticleGun *particleGun;
	G4GeneralParticleSource *theArchitect;
	G4SingleParticleSource *theSource[2];
	G4int dataEvent;
	G4double particleEnergy;

	MargotDataRecordTree* MargotDataOutPG;

	G4int numberOfSources;
	G4String sourceType;
	G4String beamName;

	G4SPSEneDistribution *GPS_ParticleEnergy;
	G4SPSPosDistribution *GPS_ParticlePosition;
	G4SPSAngDistribution *GPS_ParticleMomentum;

	TFile *rootFile;
	TH1F *theHistogram;


};

#endif /* PROTONDETECTORPRIMARYGENERATORACTION_HH_ */

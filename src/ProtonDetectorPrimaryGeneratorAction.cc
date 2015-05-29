/*
 * ProtonDetectorPrimaryGeneratorAction.cc
 *
 *  Created on: Dec 4, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorPrimaryGeneratorAction.hh>
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include "G4GeneralParticleSource.hh"
#include "G4SingleParticleSource.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"

#include "Randomize.hh"
#include "G4UnitsTable.hh"

#include "ProtonDetectorROOTAnalysis.hh"

ProtonDetectorPrimaryGeneratorAction::ProtonDetectorPrimaryGeneratorAction(G4int sources, G4String particleName, G4String type, G4double energy = 5.0*MeV):
numberOfSources(sources), sourceType(type)
{
	G4cout<< "##################################################################" <<G4endl
			<<"######  ProtonDetectorPrimaryGenerationAction::ProtonDetectorPrimaryGeneratorAction()  ########" <<G4endl;
	G4cout<< "##################################################################"<<G4endl;


	//MargotDataOutPG = MargotDataRecordTree::MargotPointer;
	G4int nParticle = 1;
	particleGun = new G4ParticleGun(nParticle);
	theArchitect = new G4GeneralParticleSource();
	theArchitect->SetNumberOfParticles(nParticle);
	theArchitect->SetMultipleVertex(true);
	for(G4int i = 0; i<numberOfSources -1; i++){
		theArchitect->AddaSource(1.0);
	}
	theArchitect->ListSource();
	beamName = particleName;
	particleEnergy = energy;
	dataEvent=0;
	GPS_ParticleEnergy = NULL;
	GPS_ParticlePosition = NULL;
	GPS_ParticleMomentum = NULL;

	rootFile=new TFile("Zpos_al23.root");
	theHistogram=(TH1F*)rootFile->Get("h_pos");
}

ProtonDetectorPrimaryGeneratorAction::~ProtonDetectorPrimaryGeneratorAction() {
	delete theArchitect;
	delete particleGun;

	//delete rootFile;
	//delete theHistogram;
}

void ProtonDetectorPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent){

	//G4cout<< "##################################################################" <<G4endl
	//		<<"######  ProtonDetectorPrimaryGenerationAction::GeneratePrimaries()  ########" <<G4endl;
	//G4cout<< "##################################################################"<<G4endl;


	dataEvent++;
	G4int Z,A;

	//Define the Beam Particle depending on particleName
	//TODO Change this for a messenger in future
	if(beamName == "Al23"){
		Z=13;
		A=23;
	}
	else if(beamName == "Mg20"){
		Z=12;
		A=20;
	}
	else if(beamName == "Cl31"){
		Z=17;
		A=31;
	}
	else if(beamName == "proton"){
		Z=1;
		A=1;
	}

	//Creates a Beam Ion with no Excitation Energy

	if(beamName != "geantino"){
		theArchitect->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));
		particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));
	}
	else{
		G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();
		theArchitect->SetParticleDefinition(theParticleTable->FindParticle("geantino"));
		particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));
	}



	theSource[0] = theArchitect->GetCurrentSource();


	GPS_ParticleEnergy = theSource[0]->GetEneDist();
	GPS_ParticleEnergy->SetEnergyDisType("Mono");
	GPS_ParticleEnergy->SetMonoEnergy(particleEnergy);

	GPS_ParticleMomentum = theSource[0]->GetAngDist();

	GPS_ParticlePosition = theSource[0]->GetPosDist();
	GPS_ParticlePosition ->SetPosDisType("Point");

	//****************************************************
	// Select type of beam
	//****************************************************
	// SourceType = "test" -> Beam on z-axis at beam energy given in the main()
	// SourceType = "pencil" -> Pencil beam with +- 5% energy variance
	// SourceType = "RealBeam" -> 5mm Gaussian in x-y beam spot with 1 deg. divergence

	  G4double momentum_X = 0;   // Default momentum is z-direction for
	  G4double momentum_Y = 0;   // Pencil beams
	  G4double momentum_Z = 1;

	  G4double theta = 0.;
	  G4double phi= 0.;

	  const G4double Pi = CLHEP::pi;
	  const G4double twopi = 2*Pi;
	  G4ThreeVector theStartPosition(0.,0.,-25.*cm);

	  GPS_ParticlePosition->SetCentreCoords(theStartPosition);


	  if(sourceType == "test")
	  {
	     GPS_ParticlePosition->SetCentreCoords(theStartPosition);
	     momentum_X = 0.;
	     momentum_Y = 0.;
	     momentum_Z = 1.;
	     //MargotDataOutPG->senddataPG(particleEnergy,theStartPosition);
	     GPS_ParticleMomentum = theSource[0]->GetAngDist();
	     G4ThreeVector v(momentum_X,momentum_Y,momentum_Z);
	     GPS_ParticleMomentum->SetParticleMomentumDirection(v);
	     theArchitect->GeneratePrimaryVertex(anEvent);
	    }

	  else if(sourceType == "pencil"){

		  //Straight Pencil Beam
		  momentum_X = 0.;
		  momentum_Y = 0.;
		  momentum_Z = 1.;
		  GPS_ParticlePosition->SetCentreCoords(theStartPosition);

		  // Now create Beam Energy Spread of +- 0.5 % !
		  G4double energySpread = 0.01*particleEnergy;
		  G4double energyLowLimit = particleEnergy-0.5*energySpread;
		  G4double deltaEnergy = energySpread*G4UniformRand();
		  G4double energyOut= energyLowLimit+deltaEnergy;
		  GPS_ParticleEnergy->SetMonoEnergy(energyOut);
		  //MargotDataOutPG->senddataPG(energyOut,theStartPosition);
	  }
	  else if(sourceType == "RealBeam"){

	      // +- .5% energy Spread
	      // +- 2.9% energy Spread LISE++
		  G4double energySpread = 0.01*particleEnergy;
		  G4double energyLowLimit = particleEnergy-2.9*energySpread;
		  G4double deltaEnergy = energySpread*G4UniformRand();
		  G4double energyOut= energyLowLimit+deltaEnergy;
		  //GPS_ParticleEnergy->SetMonoEnergy(energyOut);
		  particleGun->SetParticleEnergy(energyOut);


	      // 5 mm gaussian random starting position in x-y
		  // Positions from LISE++ calculation
		  G4double FWHM_X = 22.6*mm;
		  G4double FWHM_Y = 22.6*mm;
		  G4double mean_X = -0.1005*mm;
		  //G4double mean_Y = 6.355*mm;
		  G4double mean_Y = 0*mm;
		  G4double conversionFactor = 2.35482;   // 2*sqrt(2*log(2)) - converts FWHM to Sigma
		  G4double sigmaX =  FWHM_X/conversionFactor;
		  G4double sigmaY =  FWHM_Y/conversionFactor;
		  theStartPosition[0]=CLHEP::RandGauss::shoot(mean_X,sigmaX)*mm;
		  theStartPosition[1]=CLHEP::RandGauss::shoot(mean_Y,sigmaY)*mm;
		  theStartPosition[2]=0.5*cm;
		  //theStartPosition[2]=-30.5*cm;
		  //GPS_ParticlePosition->SetCentreCoords(theStartPosition);
		  particleGun->SetParticlePosition(theStartPosition);

	      // 1 deg. divergence (from starting point)
		  G4double beamDivergence = 1.*deg;
		  G4double ran = G4UniformRand();
		  // isotropic in cosine for OpenAngle (in degrees)
		  theta = acos(1+(cos(beamDivergence)-1)*ran);
		  ran = G4UniformRand();
		  phi = twopi*ran;

		  momentum_X = sin(theta)*cos(phi);
		  momentum_Y = sin(theta)*sin(phi);
		  momentum_Z = cos(theta);

		  G4ThreeVector v(momentum_X ,momentum_Y ,momentum_Z);
		  //GPS_ParticleMomentum->SetParticleMomentumDirection(v);
		  particleGun->SetParticleMomentumDirection(v);
		  //MargotDataOutPG->senddataPG(energyOut,theStartPosition);
		  particleGun->GeneratePrimaryVertex(anEvent);
	  }
	  else if(sourceType == "General"){

		  for(G4int i=0; i<numberOfSources;i++){
			  theArchitect->SetCurrentSourceto(i);
			  theSource[i]=theArchitect->GetCurrentSource();
			  theSource[i]->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));
			  //particleEnergy = 1800.*keV;
			  particleEnergy = 0.05*keV;
			  GPS_ParticleEnergy->SetEnergyDisType("Mono");
			  GPS_ParticleEnergy->SetMonoEnergy(particleEnergy);
			  // 9 mm gaussian random starting position in x-y
			  G4double FWHM_X = 26.33*mm;
			  G4double FWHM_Y = 5.178*mm;
			  //G4double FWHM_X = 0*mm;
			  //G4double FWHM_Y = 0*mm;
			  G4double mean_X = -2.513*mm;
			  G4double mean_Y =  0.02015*mm;
			  //G4double mean_Y = 0*mm;
			  G4double conversionFactor = 2.35482;   // 2*sqrt(2*log(2)) - converts FWHM to Sigma
			  //G4double sigma =  FWHM/conversionFactor;
			  G4double sigmaX =  FWHM_X/conversionFactor;
			  G4double sigmaY =  FWHM_Y/conversionFactor;
			  //G4double centroidz = 0.*mm;
			  theStartPosition[0]=CLHEP::RandGauss::shoot(mean_X,sigmaX)*mm;
			  theStartPosition[1]=CLHEP::RandGauss::shoot(mean_Y,sigmaY)*mm;
			  // 24 mm gaussian random starting position in z 800 torr
			  //G4double FWHMz = 53.54*mm;
			  G4double FWHMz = 0*mm;
			  G4double DeltaZ =  FWHMz*(1-2*G4UniformRand());
			  G4double sigmaz = FWHMz/conversionFactor;
			  //G4double centroidz = 49.35*mm;
			  G4double centroidz = -1.01*mm;
			  //theStartPosition[2]=CLHEP::RandGauss::shoot(centroidz,sigmaz)*mm;
			  //theStartPosition[2]=centroidz;
			  theStartPosition[2]=theHistogram->GetRandom()*mm;

			  GPS_ParticlePosition = theSource[i]->GetPosDist();
			  GPS_ParticlePosition->SetPosDisType("Point");
			  GPS_ParticlePosition->SetCentreCoords(theStartPosition);


			  if(beamName == "proton"){
				  // Use isotropic momentum dist. for proton eff. checks
				  G4double OpenAngle = 180.*deg;
				  G4double ran = G4UniformRand();
				  // isotropic in cosine for OpenAngle (in degrees)
				  theta = acos(1+(cos(OpenAngle)-1)*ran);
				  ran = G4UniformRand();
				  phi = twopi*ran;
			  }
			  momentum_X = sin(theta)*cos(phi);
			  momentum_Y = sin(theta)*sin(phi);
			  momentum_Z = cos(theta);


			  G4ThreeVector v(momentum_X,momentum_Y,momentum_Z);

			  //MargotDataOutPG->senddataPG(particleEnergy,theStartPosition);

			  GPS_ParticleMomentum = theSource[i]->GetAngDist();
			  GPS_ParticleMomentum->SetParticleMomentumDirection(v);
		  }
		  theArchitect->GeneratePrimaryVertex(anEvent);

		  }

	  if(gProtonDetectorROOTAnalysis){
		//G4cout<<"HERE!!"<<G4endl;
		//G4PrimaryVertex *theVertex=anEvent->GetPrimaryVertex();
		//G4cout<<"theVertex "<<theVertex<<G4endl;
		gProtonDetectorROOTAnalysis->GenerateBeam(anEvent);
	  }
}



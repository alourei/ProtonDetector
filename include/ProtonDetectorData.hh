/*
 * ProtonDetectorData.hh
 *
 *  Created on: Dec 19, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORDATA_HH_
#define PROTONDETECTORDATA_HH_

#include "TROOT.h"

class ProtonDetectorData {

private:

	Double_t energyOnGas_beam;
	Double_t energyOnGas_protons;
	Double_t energyOnGas_betas;
	Double_t energyOnGas_total;
	Double_t energyOnDegrader;
	Double_t stepSumLengthOnGas_beam;
	Double_t stepSumLengthOnGas_protons;
	Double_t stepSumLengthOnGas_betas;
	Double_t beam_last_positionX;
	Double_t beam_last_positionY;
	Double_t beam_last_positionZ;
	Int_t eventID;
	Int_t runID;

public:
	ProtonDetectorData();
	~ProtonDetectorData();

	Double_t GetEnergyOnGas_beam(){return energyOnGas_beam;}
	Double_t GetEnergyOnGas_protons(){return energyOnGas_protons;}
	Double_t GetEnergyOnGas_betas(){return energyOnGas_betas;}
	Double_t GetEnergyOnGas_total(){return energyOnGas_total;}
	Double_t GetEnergyOnDegrader(){return energyOnDegrader;}
	Double_t GetStepSumLengthOnGas_beam(){return stepSumLengthOnGas_beam;}
	Double_t GetStepSumLengthOnGas_protons(){return stepSumLengthOnGas_protons;}
	Double_t GetStepSumLengthOnGas_betas(){return stepSumLengthOnGas_betas;}
	Double_t GetBeamLastPositionX(){return beam_last_positionX;}
	Double_t GetBeamLastPositionY(){return beam_last_positionY;}
	Double_t GetBeamLastPositionZ(){return beam_last_positionZ;}
	Int_t GetEventID(){return eventID;}
	Int_t GetRunID(){return runID;}

	void SetEnergyOnGas_beam(Double_t val){energyOnGas_beam=val;}
	void SetEnergyOnGas_protons(Double_t val){energyOnGas_protons=val;}
	void SetEnergyOnGas_betas(Double_t val){energyOnGas_betas=val;}
	void SetEnergyOnGas_total(Double_t val){energyOnGas_total=val;}
	void SetEnergyOnDegrader(Double_t val){energyOnDegrader=val;}
	void SetStepSumLengthOnGas_beam(Double_t val){  stepSumLengthOnGas_beam=val;}
	void SetStepSumLengthOnGas_protons(Double_t val){  stepSumLengthOnGas_protons=val;}
	void SetStepSumLengthOnGas_betas(Double_t val){  stepSumLengthOnGas_betas=val;}
	void SetBeamLastPositionX(Double_t val){  beam_last_positionX=val;}
	void SetBeamLastPositionY(Double_t val){  beam_last_positionY=val;}
	void SetBeamLastPositionZ(Double_t val){  beam_last_positionZ=val;}
	void SetEventID(Int_t val){  eventID =val;}
	void SetRunID(Int_t val){  runID =val;}



ClassDef(ProtonDetectorData,1);
};

#endif /* PROTONDETECTORDATA_HH_ */

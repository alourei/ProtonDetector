/*
 * ProtonDetectorData.cc
 *
 *  Created on: Dec 19, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorData.hh>

ClassImp(ProtonDetectorData)

ProtonDetectorData::ProtonDetectorData() {
	// TODO Auto-generated constructor stub

	 energyOnGas_beam =0;
	 energyOnGas_protons =0;
	 energyOnGas_betas =0;
	 energyOnGas_total =0;
	 energyOnDegrader =0;
	 stepSumLengthOnGas_beam =0;
	 stepSumLengthOnGas_protons =0;
	 stepSumLengthOnGas_betas =0;
	 beam_last_positionX =-999;
	 beam_last_positionY =-999;
	 beam_last_positionZ =-999;
	 eventID =0;
	 runID =0;

}

ProtonDetectorData::~ProtonDetectorData() {
	// TODO Auto-generated destructor stub
}


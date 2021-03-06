/*
 * ProtonDetectorBeamInfo.hh
 *
 *  Created on: Dec 16, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORBEAMINFO_HH_
#define PROTONDETECTORBEAMINFO_HH_

#include <TObject.h>
#include <TROOT.h>

class ProtonDetectorBeamInfo: public TObject {
private:

	Double_t energyEntrance;

	Double_t thetaEntrance;      // theta emission angle

	Double_t phiEntrance;        // phi emission angle

	Double_t xInitial;          // beam entrance position
	Double_t yInitial;
	Double_t zInitial;

	//Double_t Charge;
	//Double_t Mass;

public:
	ProtonDetectorBeamInfo();
	~ProtonDetectorBeamInfo();

	  inline Double_t GetEnergyEntrance() const { return energyEntrance; }

	  inline Double_t GetThetaEntrance() const { return thetaEntrance; }
	  inline Double_t GetPhiEntrance() const { return phiEntrance; }


	  inline Double_t GetXInitial() const { return xInitial; }
	  inline Double_t GetYInitial() const { return yInitial; }
	  inline Double_t GetZInitial() const { return zInitial; }

	  //inline Double_t GetCharge() const { return Charge; }
	  //inline Double_t GetMass() const { return Mass; }

	  inline void SetEnergyEntrance(Double_t e) { energyEntrance = e; }

	  inline void SetThetaEntrance(Double_t t) { thetaEntrance = t; }
	  inline void SetPhiEntrance(Double_t p) { phiEntrance = p; }

	  inline void SetXInitial(Double_t x) { xInitial = x; }
	  inline void SetYInitial(Double_t y) { yInitial = y; }
	  inline void SetZInitial(Double_t z) { zInitial = z; }

	  //inline void SetCharge(Double_t ch) { Charge=ch; }
	  //inline void SetMass(Double_t m)  { Mass=m; }

	  void print(void);

	  ClassDef(ProtonDetectorBeamInfo,1) //ROOT CINT

};

#endif /* PROTONDETECTORBEAMINFO_HH_ */

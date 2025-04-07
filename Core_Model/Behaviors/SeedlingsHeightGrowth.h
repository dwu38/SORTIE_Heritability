/*
 * SeedlingsHeightGrowth.h
 *
 *  Created on: Mar 21, 2018
 *      Author: ucmuser
 */

#ifndef BEHAVIORS_SEEDLINGSHEIGHTGROWTH_H_
#define BEHAVIORS_SEEDLINGSHEIGHTGROWTH_H_

#include "GrowthBase.h"

class clTree;
class clTreePopulation;
//class clGrid;

class clSeedlingsHeightGrowth: virtual public clGrowthBase {
//note: need the virtual keyword to avoid base class ambiguity.

public:

	/**
	 * Constructor.  Sets the namestring.
	 */
	clSeedlingsHeightGrowth(clSimManager *p_oSimManager);

	/**
	 * Destructor.  Frees memory.
	 */
	~clSeedlingsHeightGrowth();

	/**
	 * Calculates the new size class
	 *
	 * @param p_oTree Tree for which to calculate growth.
	 * @param p_oPop Tree population object, just in case it's needed.
	 * @param fDiameterGrowth Amount of diameter growth for this tree, in cm.
	 * @return Amount, in m, by which to increase the tree's height.
	 */
	int CalcSizeClassGrowthValue(clTree *p_oTree, clTreePopulation *p_oPop,
			float fDiameterGrowth);

	/**
	 * Does the setup for this behavior.  This reads in the parameters from the
	 * parameter file.
	 * @param p_oDoc Parsed parameter file.
	 */
	void DoShellSetup(xercesc::DOMDocument *p_oDoc);

protected:

	/**Intercept - sized number species*/
	float *mp_fInter;

	/**Xtra Small seedlings (first and second year) - sized number species*/
	float *mp_fXSm;

	/**Small seedlings (<10) - sized number species*/
	float *mp_fSm;

	/**Medium seedlings (10-50cm) - sized number species*/
	float *mp_fMd;

	/**large seedlings (50-140cm) - sized number species*/
	float *mp_fLg;

	/**Effect of the July maximum temperature - sized number species*/
	float *mp_fJMxC;

	/**Effect of the precipitation - sized number species*/
	float *mp_fPC;

	/**Effect of the basal area - sized number species*/
	float *mp_fBA;

	/**Neighborhood search radius.*/
	float m_fRadius;

	/**Minimum sapling height. For doing neighbor searches.*/
	float m_fMinSaplingHeight;

	/**Number of years per timestep - from sim manager*/
	float m_fNumberYearsPerTimestep;

	/**
	 * Makes sure that the neighbor
	 * search radius >= 0.
	 *
	 * @throws modelErr if the above conditions are not met.
	 */
	void ValidateData();

	/**
	 * Gets the total adult neighborhood basal area within the specified radius
	 * from a given point. Neighbors must have a DBH greater than the minimum. They
	 * also cannot be dead from a disturbance event; but any trees that have a
	 * dead code of "natural" are assumed to have died in the current time step
	 * mortality cycle and thus should be counted.
	 * @param fX X coordinate of point for which to calculate neighborhood basal
	 * area
	 * @param fY Y coordinate of point for which to calculate neighborhood basal
	 * area
	 * @param p_oPop Tree population.
	 * @returns Total adult basal area in square meters.
	 */
	float GetBAT(float &fX, float &fY, clTreePopulation *p_oPop);

};

#endif /* BEHAVIORS_SEEDLINGSHEIGHTGROWTH_H_ */

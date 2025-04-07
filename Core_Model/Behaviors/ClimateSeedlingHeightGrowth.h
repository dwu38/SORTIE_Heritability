/*
 * ClimateSeedlingHeightGrowth.h
 *
 *  Created on: Mar 14, 2018
 *      Author: ucmuser
 */

#ifndef BEHAVIORS_CLIMATESEEDLINGHEIGHTGROWTH_H_
#define BEHAVIORS_CLIMATESEEDLINGHEIGHTGROWTH_H_

#include "GrowthBase.h"
class clGrid;

/**
 * Growth model for the seedlings
 *
 * Compute a probability of switching from onesize category to another
 * p (probability of switching) depends on the size, the climate and the basal area
 */

class clClimateSeedlingHeightGrowth : virtual public clGrowthBase {
//note: need the virtual keyword to avoid base class ambiguity.

  public:

  /**
  * Constructor.  Sets the namestring.
  */
  clClimateSeedlingHeightGrowth(clSimManager *p_oSimManager);

  /**
  * Destructor.  Frees memory.
  */
  ~clClimateSeedlingHeightGrowth();

  /**
  * Calculates the amount of Height growth increase for a particular seedling
  */
  float CalcDiameterGrowthValue(clTree *p_oTree, clTreePopulation *p_oPop, float fHeightGrowth);

  /**
  * Does the setup for this behavior.  This reads in the parameters from the
  * parameter file, and retrieves the "Storm Light" grid if present.
  * @param p_oDoc Parsed parameter file.
  */
  void DoShellSetup(xercesc::DOMDocument *p_oDoc);

  /**
  * Captures the namestring passed to this behavior.  This is overridden from
  * clBehaviorBase so we can capture the namestring passed.  Since this class
  * can create multiple kinds of behaviors that function differently, this will
  * capture what kind of behavior this is supposed to be.
  *
  * @param sNameString Behavior's namestring.
  */
  void SetNameData(std::string sNameString);


  protected:

  /**Intercept - sized number species*/
  float *mp_fInter;

	/**Xtra Small seedlings (first and second year) - sized number species*/
	float *mp_fXSm; //added seedling size 0

  /**Small seedlings (10-25cm) - sized number species*/
  float *mp_fSm;

  /**Medium seedlings (25-50cm) - sized number species*/
  float *mp_fMd;

  /**large seedlings (50-140cm) - sized number species*/
  float *mp_fLg;

  /**Effect of the July maximum temperature - sized number species*/
  float *mp_fJMxC;

  /**Effect of the snow - sized number species*/
  float *mp_fPC;

  /**Effect of the basal area - sized number species*/
  float *mp_fBA;

  /**Neighborhood search radius.*/
  float m_fRadius;

  /**Minimum sapling height. For doing neighbor searches.*/
  float m_fMinSaplingHeight;

  /**Return code for the "BAT" grid data member.*/
  short int m_iBATCode;

  /**To help access the other arrays*/
  int *mp_iIndexes;

  /**Conversion factor to translate the results of the function to the
  * appropriate units per timestep*/
  float m_fYearsPerTimestep;

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
  float GetBAT (float &fX, float &fY, clTreePopulation *p_oPop);

};
//---------------------------------------------------------------------------

#endif /* BEHAVIORS_CLIMATESEEDLINGHEIGHTGROWTH_H_ */

/*
 * ClimateGridGrowth.h
 *
 *  Created on: Mar 14, 2017
 *      Author: ucmuser
 */

#ifndef BEHAVIORS_CLIMATEGRIDGROWTH_H_
#define BEHAVIORS_CLIMATEGRIDGROWTH_H_

//---------------------------------------------------------------------------
// inspired by DoubleMMRelGrowth
//---------------------------------------------------------------------------
//#if !defined(DoubleMMRelGrowth_H)
//  #define DoubleMMRelGrowth_H
//
//#include "MichMenGrowthBase.h"

#include "GrowthBase.h"

//class clTree;
//class clTreePopulation;
class clGrid;

/**
 * ClimateGridGrowth, version 1
 * Use the same function as ClimateDepGrowth, but to make climate depend on the altitude, use a grid instead of a constant
 */

class clClimateGridGrowth: virtual public clGrowthBase {
//note: need the virtual keyword to avoid base class ambiguity.

public:

	/**
	 * Constructor.  Sets the namestring.
	 *
	 * @param p_oSimManager Sim Manager object.
	 */
	clClimateGridGrowth(clSimManager *p_oSimManager);

	/**
	 * Destructor.  Frees memory.
	 */
	~clClimateGridGrowth();

	/**
	 * Applies growth as described in the equation above.
	 *
	 * @param p_oTree Tree for which to calculate growth.
	 * @param p_oPop Tree population.
	 * @param fHeightGrowth Amount of height growth, in m (ignored).
	 * @return Amount of growth increase, in cm.
	 */
	float CalcDiameterGrowthValue(clTree *p_oTree, clTreePopulation *p_oPop,
			float fHeightGrowth);

	/**
	 * Captures the behavior name passed from the parameter file.  This is useful
	 * since this class can produce a few different kinds of behaviors.
	 *
	 * @param sNameString Behavior name from parameter file.
	 */
	void SetNameData(std::string sNameString);

	/**
	 * Does setup.  Reads in values from the parameter file, and validates that
	 * all species/type combos use light (each must have "lgm" registered).
	 *
	 * @param p_oDoc DOM tree of parsed input file.
	 * @throws modelErr if any species/type combo to which this behavior is
	 * applied does not have a light behavior, or if there is no grids called
	 * "Resource1" and "Resource2".
	 */
	void DoShellSetup(xercesc::DOMDocument *p_oDoc);

protected:

	/**Grid containing the levels of the first resource.  This is not created
	 * by this behavior; it is expected to already be available.  It should be
	 * named "Resource1" and have one float data member called "Resource1".*/
	clGrid *mp_oPrecipGrid;

	/**Grid containing the levels of the second resource.  This is not created
	 * by this behavior; it is expected to already be available.  It should be
	 * named "Resource2" and have one float data member called "Resource2".*/
	clGrid *mp_oTemperatureGrid;

	/**DBH slope of growth equation - b - sized number of behavior
	 * species*/
	float *mp_fDBHSlope;

	/**Intercept of growth equation - a - sized number of behavior
	 * species*/
	float *mp_fCdIntercept;

	/**temperature slope of growth equation - c - sized number of behavior
	 * species*/
	float *mp_fTempSlope;

	/**precipitation slope of growth equation - d - sized number of behavior
	 * species*/
	float *mp_fPrecipSlope;

	/**To help access the other arrays*/
	int *mp_iIndexes;

	/**Conversion factor to translate the results of the function to the
	 * appropriate units per timestep*/
	float m_fYearsPerTimestep;

	/**The float data member code for the first resource grid.*/
	short int m_iPrecipCode;

	/**The float data member code for the second resource grid.*/
	short int m_iTemperatureCode;
};
//---------------------------------------------------------------------------

#endif /* BEHAVIORS_CLIMATEGRIDGROWTH_H_ */

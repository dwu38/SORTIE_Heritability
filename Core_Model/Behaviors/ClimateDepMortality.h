/*
 * ClimateDepMortality.h
 *
 *  Created on: Jul 7, 2017
 *      Author: ucmuser
 */

#ifndef BEHAVIORS_CLIMATEDEPMORTALITY_H_
#define BEHAVIORS_CLIMATEDEPMORTALITY_H_

#include "MortalityBase.h"

class clGrid;
class clTreePopulation;
//class clTree;

/***
 * ClimateDepMortality
 * This evaluate mortality acording to a logistic equation,
 * linking mortality to the actual DBH and climate variables (2 variables)
 * The equation used in this behavior is:
 * <center>p = exp(a + b * D + c * var1 + d * var2) / (1 + exp(a + b * D + c * var1 + d * var2))</center>
 * where
 * <ul>
 * <li>p = annual probability of survival</li>
 * <li>a = parameter</li>
 * <li>b = parameter</li>
 * <li>D = tree diameter, in cm; diam10 for seedlings, DBH for others</li>
 * <li>c= parameter</li>
 * <li>d= parameter</li>
 * <li>var1 = climatic variable 1 </li>
 * <li>var2 = climatic variable 2 </li>
 *
 * </ul>
 *
 * This class's namestring is "climatedep mortshell".  The parameter file
 * call string is "ClimateDepMortality".
 *
 */

class clClimateDepMortality: virtual public clMortalityBase {

public:

	/* Constructor */
	clClimateDepMortality(clSimManager *p_oSimManager);

	/* Destructor*/
	~clClimateDepMortality();

	/*Reads in values from the parameter file.*/
	void DoShellSetup(xercesc::DOMDocument *p_oDoc);

	/* Calculates mortality according to the logistic equation.*/
	deadCode DoMort(clTree *p_oTree, const float &fDBH,
			const short int &iSpecies);

	/* get the climate variables */
//	void PreMortCalcs( clTreePopulation *p_oPop );
protected:

	/* Tree population*/
	clTreePopulation *mp_oPop;

	/* Mortality equation "a" - sized number of behavior species*/
	float *mp_fInter;

	/* Mortality equation "b" - sized number of behavior species*/
	float *mp_fDiam;

	/* Mortality equation "c" - sized number of behavior species*/
	float *mp_fTemp;

	/* Mortality equation "d" - sized number of behavior species*/
	float *mp_fPrecip;

	/**Climate var 1 function.*/
//	float *mp_fClimVar1Function;
	/**Climate var 2 function.*/
//	float *mp_fClimVar2Function;
	/* To help access the other arrays*/
	int *mp_iIndexes;

	/*Conversion factor to translate the results of the function to the
	 * appropriate units per timestep*/
	float m_fYearsPerTimestep;

};

#endif /* BEHAVIORS_CLIMATEDEPMORTALITY_H_ */

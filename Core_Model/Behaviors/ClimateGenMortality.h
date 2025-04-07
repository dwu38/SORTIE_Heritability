/*
 * ClimateGenMortality.h
 *
 *  Created on: May 4, 2018
 *      Author: ucmuser
 */

#ifndef BEHAVIORS_CLIMATEGENMORTALITY_H_
#define BEHAVIORS_CLIMATEGENMORTALITY_H_

#include "MortalityBase.h"

class clGrid;
class clTreePopulation;
//class clTree;

/**
 * Change ClimateDepmortality to take the tree individual parameters ParamM1 and ParamM2 into account when computing mortality probability
 */

class clClimateGenMortality: virtual public clMortalityBase {

public:

	/* Constructor */
	clClimateGenMortality(clSimManager *p_oSimManager);

	/* Destructor*/
	~clClimateGenMortality();

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

#endif /* BEHAVIORS_CLIMATEGENMORTALITY_H_ */

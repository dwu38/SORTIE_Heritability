/*
 * LogisticMortality.h
 *
 *  Created on: Jul 7, 2016
 *      Author: ucmuser
 */

#ifndef LogisticMortalityH
#define LogisticMortalityH

//---------------------------------------------------------------------------
#include "MortalityBase.h"

class clGrid;
class clTreePopulation;
/**
 * Logistic Mortality - Version 1.0
 *
 * This evaluates mortality according to a logistic equation, with the
 * possibility of two sets of parameters for each species.
 *
 * The equation used in this behavior is:
 * <center>p = exp(a + b * D) / (1 + exp(a + b * D))</center>
 * where
 * <ul>
 * <li>p = annual probability of survival</li>
 * <li>a = parameter</li>
 * <li>b = parameter</li>
 * <li>D = tree diameter, in cm; diam10 for seedlings, DBH for others</li>
 * </ul>
 *
 * This class's namestring is "logistic mortshell".  The parameter file
 * call string is "LogisticMortality".
 *
 * Copyright 2005 Charles D. Canham.
 * @author Lora E. Murphy
 *
 * <br>Edit history:
 * <br>-----------------
 * <br>July 7, 2016 - Created (LEM)
 */
class clLogisticMortality: virtual public clMortalityBase {
//note: need the virtual keyword to avoid base class ambiguity.

public:

	/**
	 * Constructor.
	 *
	 * @param p_oSimManager Sim Manager object.
	 */
	clLogisticMortality(clSimManager *p_oSimManager);

	/**
	 * Destructor.
	 */
	~clLogisticMortality();

	/**
	 * Reads in values from the parameter file.
	 *
	 * @param p_oDoc DOM tree of parsed input file.
	 */
	void DoShellSetup(xercesc::DOMDocument *p_oDoc);

	/**
	 * Calculates mortality according to the logistic equation.
	 *
	 * @param p_oTree Tree being evaluated
	 * @param fDbh Tree's DBH
	 * @param iSpecies Species of the tree being evaluated
	 * @return natural if the tree is to die, notdead if it lives.
	 */
	deadCode DoMort(clTree *p_oTree, const float &fDbh,
			const short int &iSpecies);

protected:

	/**Tree population - for getting data codes*/
	clTreePopulation *mp_oPop;

	/**Mortality equation "b" - sized number of behavior species*/
	float *mp_fB;

	/**Mortality equation "a" - sized number of behavior species*/
	float *mp_fA;

	/**To help access the other arrays*/
	int *mp_iIndexes;

	/**Conversion factor to translate the results of the function to the
	 * appropriate units per timestep*/
	float m_fYearsPerTimestep;

};
//---------------------------------------------------------------------------

#endif /* BEHAVIORS_LOGISTICMORTALITY_H_ */

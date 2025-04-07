/*
 * ParamChange.h
 *
 *  Created on: Dec 6, 2017
 *      Author: ucmuser
 */

#ifndef BEHAVIORS_PARAMCHANGE_H_
#define BEHAVIORS_PARAMCHANGE_H_

#include "BehaviorBase.h"

/*
 * ParamChange behavior
 * This behavior change the value of the mean parameter that represent a species mean 'trait'
 * The mean parameter is the one changing (it is the mean of the population parameters)
 * The mean has to evolve as the population evolve
 * then individual parameters are sampled, based on the mean
 *
 * Inspired by ClimateChange
 */

class clParamChange: virtual public clBehaviorBase {

public:

	/* Constructor
	 * @param p_oSimManager Sim Manager object.
	 */
	clParamChange(clSimManager *p_oSimManager);

	/* Destuctor
	 * - not needed
	 */

	/**
	 * Reads in values from the parameter file.
	 * @param p_oDoc DOM tree of parsed input file.
	 */
	void GetData(xercesc::DOMDocument *p_oDoc);

	/*
	 * Update the parameter according to the parameters of all trees in the plot
	 */

	void Action();

protected:

	/** heritability */
	float m_fH;

};

#endif /* BEHAVIORS_PARAMCHANGE_H_ */

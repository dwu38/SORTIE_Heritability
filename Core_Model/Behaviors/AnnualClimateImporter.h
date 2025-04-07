/*
 * AnnualClimateImporter.h
 *
 *  Created on: Jul 26, 2017
 *      Author: ucmuser
 */

#ifndef BEHAVIORS_ANNUALCLIMATEIMPORTER_H_
#define BEHAVIORS_ANNUALCLIMATEIMPORTER_H_

#include "BehaviorBase.h"

/**
 * Inspired by the class ClimateImporter
 * Class to import annual climate
 */

class clAnnualClimateImporter: virtual public clBehaviorBase {

public:

	/*
	 * Constructor
	 */
	clAnnualClimateImporter(clSimManager *p_oSimManager);

	/*
	 * Destructor
	 */
	~clAnnualClimateImporter();

	/*
	 * reads in values from the parameter file
	 */
	void GetData(xercesc::DOMDocument *p_oDoc);

	/*
	 * Update the plot climate
	 */
	void Action();

protected:

	/*
	 * January min temperature. Array length is # timesteps.
	 */
	double *mp_fJanMin;

	/*
	 * Annual precipitation. Array length is # timesteps.
	 */
	double *mp_fPpt;

	/*
	 * July Max temperature. Array length is # timesteps.
	 */
	double *mp_fJulMax;

	/*
	 * read parameter file data
	 */
	void ReadParameterFileData(DOMElement *p_oElement, double *p_fJanMin,
			double *p_fPpt, double *p_fJulMax);

	void ReadAnnualData(xercesc::DOMElement *p_oParent, std::string sParentTag,
			std::string sSubTag, double *p_fVal);

};

#endif /* BEHAVIORS_ANNUALCLIMATEIMPORTER_H_ */

/*
 * LogisticMortality.cpp
 *
 *  Created on: Jul 7, 2016
 *      Author: ucmuser
 */

//---------------------------------------------------------------------------
#include "LogisticMortality.h"
#include "MortalityOrg.h"
#include "Grid.h"
#include "ModelMath.h"
#include "ParsingFunctions.h"
#include "TreePopulation.h"
#include "SimManager.h"
#include <math.h>

////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////
clLogisticMortality::clLogisticMortality(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clMortalityBase(
				p_oSimManager) {
	try {

		m_sNameString = "logistic mortshell";
		m_sXMLRoot = "LogisticMortality";

		mp_fB = NULL;
		mp_fA = NULL;
		mp_iIndexes = NULL;

		mp_oPop = NULL;
		m_fYearsPerTimestep = 0;
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clLogisticMortality::clLogisticMortality";
		throw(stcErr);
	}
}

////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////
clLogisticMortality::~clLogisticMortality() {
	delete[] mp_fB;
	delete[] mp_fA;
	delete[] mp_iIndexes;
}

////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
////////////////////////////////////////////////////////////////////////////
void clLogisticMortality::DoShellSetup(xercesc::DOMDocument *p_oDoc) {
	floatVal *p_fTempValues = NULL; //for getting species-specific values
	try {
		mp_oPop = (clTreePopulation*) mp_oSimManager->GetPopulationObject(
				"treepopulation");

		DOMElement *p_oElement = GetParentParametersElement(p_oDoc);
		short int iNumSpecies = mp_oPop->GetNumberOfSpecies(), i;

		//Get number of years per timestep
		m_fYearsPerTimestep = mp_oSimManager->GetNumberOfYearsPerTimestep();

		//Declare the arrays we'd like read
		mp_fA = new float[m_iNumBehaviorSpecies];
		mp_fB = new float[m_iNumBehaviorSpecies];
		mp_iIndexes = new int[iNumSpecies];

		//Set up the array indexes
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_iIndexes[mp_iWhatSpecies[i]] = i;

		//Declare the species-specific temp array and pre-load with the species that
		//this behavior affects
		p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			p_fTempValues[i].code = mp_iWhatSpecies[i];

		//mortality "a"
		FillSpeciesSpecificValue(p_oElement, "mo_logteA", "mo_laVal",
				p_fTempValues, m_iNumBehaviorSpecies, mp_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fA[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		// mortality "b"
		FillSpeciesSpecificValue(p_oElement, "mo_logteB", "mo_lbVal",
				p_fTempValues, m_iNumBehaviorSpecies, mp_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fB[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		delete[] p_fTempValues;
	} catch (modelErr &err) {
		delete[] p_fTempValues;
		throw(err);
	} catch (modelMsg &msg) {
		delete[] p_fTempValues;
		throw(msg);
	} //non-fatal error
	catch (...) {
		delete[] p_fTempValues;
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clLogisticMortality::DoShellSetup";
		throw(stcErr);
	}
}

///////////////////////////////////////////////////////////////////////////
// DoMort()
///////////////////////////////////////////////////////////////////////////
deadCode clLogisticMortality::DoMort(clTree *p_oTree, const float &fDiam,
		const short int &iSpecies) {

	float fSurvivalProb, fMortDiam;
	int iTp = p_oTree->GetType();

	//If the tree is a seedling, get the tree's diam10
	if (iTp == clTreePopulation::seedling && 0 == fDiam)
		p_oTree->GetValue(mp_oPop->GetDiam10Code(iSpecies, iTp), &fMortDiam);
	else
		fMortDiam = fDiam;

	//Calculate the function value
	fSurvivalProb = exp(
			mp_fA[mp_iIndexes[iSpecies]]
					+ mp_fB[mp_iIndexes[iSpecies]] * fMortDiam);

	fSurvivalProb = fSurvivalProb / (1 + fSurvivalProb);
	//Compound by number of years per timestep
	if (1 != m_fYearsPerTimestep)
		fSurvivalProb = pow(fSurvivalProb, m_fYearsPerTimestep);

	if (clModelMath::GetRand() < fSurvivalProb)
		return notdead;
	else
		return natural;
}


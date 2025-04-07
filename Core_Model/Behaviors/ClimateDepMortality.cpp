/*
 * ClimateDepMortality.cpp
 *
 *  Created on: Jul 7, 2017
 *      Author: ucmuser
 */

#include "ClimateDepMortality.h"
#include "MortalityOrg.h"
#include "Grid.h"
#include "ModelMath.h"
#include "ParsingFunctions.h"
#include "TreePopulation.h"
#include "SimManager.h"
#include "Plot.h"
#include <math.h>
#include <stdio.h>
#include <sstream>

// Constructor

clClimateDepMortality::clClimateDepMortality(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clMortalityBase(
				p_oSimManager) {
	try {

		m_sNameString = "climatedep mortshell";
		m_sXMLRoot = "ClimateDepMortality";

		mp_fInter = NULL;
		mp_fDiam = NULL;
		mp_fTemp = NULL;
		mp_fPrecip = NULL;

//		mp_fClimVar1Function = NULL;
//		mp_fClimVar2Function = NULL;

		mp_iIndexes = NULL;

		mp_oPop = NULL;
		m_fYearsPerTimestep = 0;
	}

	catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clClimateDepMortality::clClimateDepMortality";
		throw(stcErr);
	}
}

// Destructor

clClimateDepMortality::~clClimateDepMortality() {
	delete[] mp_fInter;
	delete[] mp_fDiam;
	delete[] mp_fTemp;
	delete[] mp_fPrecip;
	delete[] mp_iIndexes;
}

// DoShellSetup

void clClimateDepMortality::DoShellSetup(xercesc::DOMDocument *p_oDoc) {
	floatVal *p_fTempValues = NULL; // for getting species-specific values
	try {
		mp_oPop = (clTreePopulation*) mp_oSimManager->GetPopulationObject(
				"treepopulation");

		DOMElement *p_oElement = GetParentParametersElement(p_oDoc);
		short int iNumSpecies = mp_oPop->GetNumberOfSpecies(), i;

		// Get number of years per timestep
		m_fYearsPerTimestep = mp_oSimManager->GetNumberOfYearsPerTimestep();

		// Declare the arrays we'd like read
		mp_fInter = new float[m_iNumBehaviorSpecies];
		mp_fDiam = new float[m_iNumBehaviorSpecies];
		mp_fTemp = new float[m_iNumBehaviorSpecies];
		mp_fPrecip = new float[m_iNumBehaviorSpecies];
		mp_iIndexes = new int[iNumSpecies];

		// Set up the array indexes
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_iIndexes[mp_iWhatSpecies[i]] = i;

		// Declare the species-specific temp array and pre-load with the species that this behavior affects
		p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			p_fTempValues[i].code = mp_iWhatSpecies[i];

		// growth "a"
		FillSpeciesSpecificValue(p_oElement, "mo_climInter", "mo_clInterVal",
				p_fTempValues, m_iNumBehaviorSpecies, mp_oPop, true);
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fInter[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		// growth "b"
		FillSpeciesSpecificValue(p_oElement, "mo_climDiam", "mo_clDiamVal",
				p_fTempValues, m_iNumBehaviorSpecies, mp_oPop, true);
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fDiam[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		// growth "c"
		FillSpeciesSpecificValue(p_oElement, "mo_climTemp", "mo_clTempVal",
				p_fTempValues, m_iNumBehaviorSpecies, mp_oPop, true);
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fTemp[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		// growth "d"
		FillSpeciesSpecificValue(p_oElement, "mo_climPrecip", "mo_clPrecipVal",
				p_fTempValues, m_iNumBehaviorSpecies, mp_oPop, true);
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fPrecip[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		delete[] p_fTempValues;
	} catch (modelErr &err) {
		delete[] p_fTempValues;
		throw(err);
	} catch (modelMsg &msg) {
		delete[] p_fTempValues;
		throw(msg);
	} catch (...) {
		delete[] p_fTempValues;
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clClimateDepMortality::DoShellSetup";
		throw(stcErr);
	}
}

// DoMort

deadCode clClimateDepMortality::DoMort(clTree *p_oTree, const float &fDiam,
		const short int &iSpecies) {

	//clTreePopulation * p_oPop = ( clTreePopulation * ) mp_oSimManager->GetPopulationObject( "treepopulation" );
	clPlot *p_oPlot = mp_oSimManager->GetPlotObject();
	float fSurvivalProb, fMortDiam, fPlotClim1 = p_oPlot->GetMeanAnnualTemp(),
			fPlotClim2 = p_oPlot->GetMeanAnnualPrecip();
	int iTp = p_oTree->GetType();

	// if the tree is a seedling, get the tree's diam10
	if (iTp == clTreePopulation::seedling && 0 == fDiam)
		p_oTree->GetValue(mp_oPop->GetDiam10Code(iSpecies, iTp), &fMortDiam);
	else
		fMortDiam = fDiam;

	//Calculate the function value
	fSurvivalProb = exp(
			mp_fInter[mp_iIndexes[iSpecies]]
					+ mp_fDiam[mp_iIndexes[iSpecies]] * fMortDiam
					+ mp_fTemp[mp_iIndexes[iSpecies]] * fPlotClim1
					+ mp_fPrecip[mp_iIndexes[iSpecies]] * fPlotClim2);
	fSurvivalProb = fSurvivalProb / (1 + fSurvivalProb);

	// Compound by number of years per timestep
	if (1 != m_fYearsPerTimestep)
		fSurvivalProb = pow(fSurvivalProb, m_fYearsPerTimestep);

	if (clModelMath::GetRand() < fSurvivalProb)
		return notdead;
	else
		return natural;
}


/*
 * SeedlingsGrowth.cpp
 *
 *  Created on: Mar 14, 2018
 *      Author: ucmuser
 */
#include "SeedlingsGrowth.h"
#include "TreePopulation.h"
#include "SimManager.h"
#include "ParsingFunctions.h"
#include "Grid.h"
#include "Plot.h"
#include "ModelMath.h"
#include "GrowthOrg.h"
#include <sstream>

//////////////////////////////////////////////////////////////////////////////
// Constructor
/////////////////////////////////////////////////////////////////////////////*/
clSeedlingsGrowth::clSeedlingsGrowth(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clGrowthBase(
				p_oSimManager) {

	try {
		mp_fInter = NULL;
		mp_fSm = NULL;
		mp_fMd = NULL;
		mp_fLg = NULL;
		mp_fJMxC = NULL;
		mp_fPC = NULL;
		mp_fBA = NULL;

		// mp_oGrid = NULL;
		m_fRadius = 10;
		m_iBATCode = -1;
		m_fMinSaplingHeight = 0;

		mp_iIndexes = NULL;
		m_fYearsPerTimestep = 0;

		m_sNameString = "Seedlingsgrowthshell";
		m_sXMLRoot = "SeedlingsGrowth";
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clSeedlingsGrowth::clSeedlingsGrowth";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////////////////////////////////////*/
clSeedlingsGrowth::~clSeedlingsGrowth() {
	delete[] mp_fInter;
	delete[] mp_fSm;
	delete[] mp_fMd;
	delete[] mp_fLg;
	delete[] mp_fJMxC;
	delete[] mp_fPC;
	delete[] mp_fBA;
	delete[] mp_iIndexes;
}

//////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
/////////////////////////////////////////////////////////////////////////////*/
void clSeedlingsGrowth::DoShellSetup(DOMDocument *p_oDoc) {
	floatVal *p_fTempValues = NULL; //for getting species-specific values
	try {
		clTreePopulation *p_oPop =
				(clTreePopulation*) mp_oSimManager->GetPopulationObject(
						"treepopulation");
		DOMElement *p_oElement = GetParentParametersElement(p_oDoc);
		short int iNumSpecies = p_oPop->GetNumberOfSpecies(), i;

		//Get number of years per timestep
		m_fYearsPerTimestep = mp_oSimManager->GetNumberOfYearsPerTimestep();

		//Declare the arrays we'd like read
		mp_fInter = new float[m_iNumBehaviorSpecies];
		mp_fSm = new float[m_iNumBehaviorSpecies];
		mp_fMd = new float[m_iNumBehaviorSpecies];
		mp_fLg = new float[m_iNumBehaviorSpecies];
		mp_fJMxC = new float[m_iNumBehaviorSpecies];
		mp_fPC = new float[m_iNumBehaviorSpecies];
		mp_fBA = new float[m_iNumBehaviorSpecies];
		mp_iIndexes = new int[iNumSpecies];

		//Set up the array indexes
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_iIndexes[mp_iWhatSpecies[i]] = i;

		//Declare the species-specific temp array and pre-load with the species that
		//this behavior affects
		p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			p_fTempValues[i].code = mp_iWhatSpecies[i];

		// Intercept
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsIntercept",
				"gr_SIVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fInter[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		// Small seedlings parameter
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsSm", "gr_SSmVal",
				p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fSm[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		// Medium seedlings parameter
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsMd", "gr_SMdVal",
				p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fMd[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		// Large seedlings parameter
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsLg", "gr_SLgVal",
				p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fLg[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		// July max parameter
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsJMxC", "gr_SJMxCVal",
				p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fJMxC[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		// Precip max parameter
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsPC", "gr_SPCVal",
				p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fPC[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		// Basal area max parameter
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsBA", "gr_SBAVal",
				p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fBA[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		delete[] p_fTempValues;
	}

	catch (modelErr &err) {
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
		stcErr.sFunction = "clSeedlingsGrowth::DoShellSetup";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// SetNameData()
/////////////////////////////////////////////////////////////////////////////*/
void clSeedlingsGrowth::SetNameData(std::string sNameString) {
	try {

		//Check the string passed and set the flags accordingly
		if (sNameString.compare("SeedlingsGrowth") == 0) {
			m_iGrowthMethod = diameter_auto;
		} else if (sNameString.compare("SeedlingsGrowth diam only") == 0) {
			m_iGrowthMethod = diameter_only;
		} else {
			modelErr stcErr;
			stcErr.iErrorCode = BAD_DATA;
			std::stringstream s;
			s << "Unrecognized behavior name \"" << sNameString << "\".";
			stcErr.sFunction = "clSeedlingsGrowth::SetNameData";
			stcErr.sMoreInfo = s.str();
			throw(stcErr);
		}
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clSeedlingsGrowth::SetNameData";
		throw(stcErr);
	}
}

////////////////////////////////////////////////////////////////////////////
// GetBAT
////////////////////////////////////////////////////////////////////////////
float clSeedlingsGrowth::GetBAT(float &fX, float &fY,
		clTreePopulation *p_oPop) {
	clTreeSearch *p_oAllNeighbors; //neighborhood trees within crowding radius
	clTree *p_oNeighbor; //competing neighbor
	char cQuery[75]; //format search strings into this
	float fDbh, //neighbor's dbh
			fBAT = 0;
	int iIsDead; //whether a neighbor is dead
	short int iSpecies, iType, //species and type for neighbor
			iDeadCode; //neighbor's dead code

	//Get all trees taller than seedlings within the max crowding radius -
	//seedlings don't compete
	sprintf(cQuery, "%s%f%s%f%s%f%s%f", "distance=", m_fRadius, "FROM x=", fX,
			"y=", fY, "::height=", m_fMinSaplingHeight);
	p_oAllNeighbors = p_oPop->Find(cQuery);

	//Loop through and find all the adults
	p_oNeighbor = p_oAllNeighbors->NextTree();

	while (p_oNeighbor) {

		iType = p_oNeighbor->GetType();
		if (clTreePopulation::adult == iType) {
			iSpecies = p_oNeighbor->GetSpecies();
			//Get the neighbor's dbh
			p_oNeighbor->GetValue(p_oPop->GetDbhCode(iSpecies, iType), &fDbh);

			//Make sure the neighbor's not dead due to a disturbance event
			iDeadCode = p_oPop->GetIntDataCode("dead", iSpecies, iType);
			if (-1 != iDeadCode)
				p_oNeighbor->GetValue(iDeadCode, &iIsDead);
			else
				iIsDead = notdead;

			if (notdead == iIsDead || natural == iIsDead)
				fBAT += clModelMath::CalculateBasalArea(fDbh);

		}

		p_oNeighbor = p_oAllNeighbors->NextTree();
	}
	return fBAT;
}

//////////////////////////////////////////////////////////////////////////////
// CalcDiameterGrowthValue()
/////////////////////////////////////////////////////////////////////////////*/
float clSeedlingsGrowth::CalcDiameterGrowthValue(clTree *p_oTree,
		clTreePopulation *p_oPop, float fHeightGrowth) {

	clPlot *p_oPlot = mp_oSimManager->GetPlotObject();
	float fX, fY, fDiam, result, fTransProb, //Probability of transition from one size class to another
			fPlotTemp = p_oPlot->GetMeanAnnualTemp(), //Current plot temperature
			fPlotPrecip = p_oPlot->GetMeanAnnualPrecip(), //Current plot precipitation;
			fDBHEffect, // DBH effect
			fBAT; //Neighborhood basal area of adult trees
	int iSp = p_oTree->GetSpecies(), iTp = p_oTree->GetType();

	// Get the tree's diameter
	p_oTree->GetValue(mp_oGrowthOrg->GetDiamCode(iSp, iTp), &fDiam);

	// Get basal area of neighbors:
	p_oTree->GetValue(p_oPop->GetXCode(iSp, iTp), &fX);
	p_oTree->GetValue(p_oPop->GetYCode(iSp, iTp), &fY);
	fBAT = GetBAT(fX, fY, p_oPop);

	//Calculate the probability of changing size class

	// effect of diameter (at 10 cm)
	if (fDiam < 0.25) {
		fDBHEffect = mp_fSm[iSp];
	} else if (fDiam < 0.67) {
		fDBHEffect = mp_fMd[iSp];
	} else {
		fDBHEffect = mp_fLg[iSp];
	}

	fTransProb = 1
			/ (1
					+ exp(
							-(mp_fInter[iSp] + fDBHEffect
									+ mp_fJMxC[iSp] * fPlotTemp
									+ mp_fPC[iSp] * fPlotPrecip
									+ mp_fBA[iSp] * fBAT)));

	if (clModelMath::GetRand() >= fTransProb) {
		// stay in the same size class
		if (fDiam < 0.25) {
			// growth between 0 and 0.25-diam
			result = clModelMath::GetRand() * (0.25 - fDiam);
		} else if (fDiam < 0.67) {
			// growth between 0 and 0.67-diam
			result = clModelMath::GetRand() * (0.67 - fDiam);
		} else {
			// growth between 0 and 2.2-diam
			result = clModelMath::GetRand() * (2.2 - fDiam);
		}
	} else {
		// transition from one size class to another
		if (fDiam < 0.25) {
			// growth between 0.25-diam and 0.67-diam
			result = clModelMath::GetRand() * (0.67 - 0.25) + 0.25 - fDiam;
		} else if (fDiam < 0.67) {
			// growth between 0.67-diam and 2.2-diam
			result = clModelMath::GetRand() * (2.2 - 0.67) + 0.67 - fDiam;
		} else {
			// growth between 2.2-diam and 3-diam (? how big can the seedling become?)
			result = clModelMath::GetRand() * (3 - 2.2) + 2.2 - fDiam;
		}
	}

	return m_fYearsPerTimestep * result;

}

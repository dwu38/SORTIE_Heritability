/*
 * SeedlingsHeightGrowth.cpp
 *
 *  Created on: Mar 21, 2018
 *      Author: ucmuser
 */

#include "SeedlingsHeightGrowth.h"
#include "TreePopulation.h"
#include "SimManager.h"
#include "ParsingFunctions.h"
#include "Grid.h"
#include "GrowthOrg.h"
#include "Plot.h"
#include "ModelMath.h"
#include <math.h>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////
// Constructor
/////////////////////////////////////////////////////////////////////////////*/
clSeedlingsHeightGrowth::clSeedlingsHeightGrowth(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clGrowthBase(
				p_oSimManager) {

	try {
		mp_fInter = NULL;
		mp_fXSm = NULL;
		mp_fSm = NULL;
		mp_fMd = NULL;
		mp_fLg = NULL;
		mp_fJMxC = NULL;
		mp_fPC = NULL;
		mp_fBA = NULL;

		m_fRadius = 10;
		m_fMinSaplingHeight = 0;

		m_iGrowthMethod = height_only;
		m_fNumberYearsPerTimestep = 0;

		m_sNameString = "seedlingsheightgrowthshell";
		m_sXMLRoot = "SeedlingsHeightGrowth";

	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clSeedlingsHeightGrowth::clSeedlingsHeightGrowth";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////////////////////////////////////*/
clSeedlingsHeightGrowth::~clSeedlingsHeightGrowth() {

	delete[] mp_fInter;
	delete[] mp_fXSm;
	delete[] mp_fSm;
	delete[] mp_fMd;
	delete[] mp_fLg;
	delete[] mp_fJMxC;
	delete[] mp_fPC;
	delete[] mp_fBA;
}

//////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
/////////////////////////////////////////////////////////////////////////////*/
void clSeedlingsHeightGrowth::DoShellSetup(DOMDocument *p_oDoc) {
	floatVal *p_fTempValues = NULL; //for getting species-specific values
	try {
		clTreePopulation *p_oPop =
				(clTreePopulation*) mp_oSimManager->GetPopulationObject(
						"treepopulation");
		DOMElement *p_oElement = GetParentParametersElement(p_oDoc);
		short int iNumSpecies = p_oPop->GetNumberOfSpecies(), i;

		m_fMinSaplingHeight = 50;
		//Get the minimum sapling height
		for (i = 0; i < iNumSpecies; i++)
			if (p_oPop->GetMaxSeedlingHeight(i) < m_fMinSaplingHeight)
				m_fMinSaplingHeight = p_oPop->GetMaxSeedlingHeight(i);

		//Declare the arrays we'd like read
		mp_fInter = new float[iNumSpecies];
		mp_fXSm = new float[iNumSpecies];
		mp_fSm = new float[iNumSpecies];
		mp_fMd = new float[iNumSpecies];
		mp_fLg = new float[iNumSpecies];
		mp_fJMxC = new float[iNumSpecies];
		mp_fPC = new float[iNumSpecies];
		mp_fBA = new float[iNumSpecies];

		//Declare the species-specific temp array and pre-load with the species that
		//this behavior affects
		p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			p_fTempValues[i].code = mp_iWhatSpecies[i];

		//Intercept
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsHeightInter",
				"gr_shiVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);

		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fInter[p_fTempValues[i].code] = p_fTempValues[i].val;

		//Xtra Small seedlings
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsHeightXSm",
				"gr_shxsmVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);

		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fXSm[p_fTempValues[i].code] = p_fTempValues[i].val;

		//Small seedlings
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsHeightSm",
				"gr_shsmVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);

		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fSm[p_fTempValues[i].code] = p_fTempValues[i].val;

		//Medium seedlings
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsHeightMd",
				"gr_shmdVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);

		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fMd[p_fTempValues[i].code] = p_fTempValues[i].val;

		//Large seedlings
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsHeightLg",
				"gr_shlgVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);

		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fLg[p_fTempValues[i].code] = p_fTempValues[i].val;

		//July maximum temperature (current)
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsHeightJMxC",
				"gr_shjmxcVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);

		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fJMxC[p_fTempValues[i].code] = p_fTempValues[i].val;

		//Current precipitation
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsHeightPC",
				"gr_shpcVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);

		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fPC[p_fTempValues[i].code] = p_fTempValues[i].val;

		//Basal Area
		FillSpeciesSpecificValue(p_oElement, "gr_SeedlingsHeightBA",
				"gr_shbaVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);

		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fBA[p_fTempValues[i].code] = p_fTempValues[i].val;

		//Max crowding radius
		FillSingleValue(p_oElement, "gr_SeedlingsHeightR", &m_fRadius, true);

		m_fNumberYearsPerTimestep =
				mp_oSimManager->GetNumberOfYearsPerTimestep();

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
		stcErr.sFunction = "clSeedlingsHeightGrowth::GetNameData";
		throw(stcErr);
	}
}

////////////////////////////////////////////////////////////////////////////
// ValidateData
////////////////////////////////////////////////////////////////////////////
void clSeedlingsHeightGrowth::ValidateData() {
	try {

		//Make sure that the neighbor search radius is not less than zero
		if (0 > m_fRadius) {
			modelErr stcErr;
			stcErr.sMoreInfo = "Neighbor search radius cannot be less than 0.";
			throw(stcErr);
		}
	} catch (modelErr &err) {
		err.sFunction = "clSeedlingsHeightGrowth::ValidateData";
		err.iErrorCode = BAD_DATA;
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clSeedlingsHeightGrowth::ValidateData";
		throw(stcErr);
	}
}

////////////////////////////////////////////////////////////////////////////
// GetBAT
////////////////////////////////////////////////////////////////////////////
float clSeedlingsHeightGrowth::GetBAT(float &fX, float &fY,
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
// CalcHeightGrowthValue()
//////////////////////////////////////////////////////////////////////////////
int clSeedlingsHeightGrowth::CalcSizeClassGrowthValue(clTree *p_oTree,
		clTreePopulation *p_oPop, float fDiameterGrowth) {
	clPlot *p_oPlot = mp_oSimManager->GetPlotObject();
	float fX, fY, fBAT, //Neighborhood basal area of adult trees
//  	  fHeight, //tree's height value
//      fNewHeight, //tree's new height value
			fHEffect, // Effect of the size class
			fTransProb, //probability of transition from one size class to another
			fPlotJulMaxTemp = p_oPlot->GetJulMaxTemp(), //Current plot July max temp
			fPlotPrecip = p_oPlot->GetMeanAnnualPrecip(); //Current plot precipitation;
//     fAmountHeightIncrease; //amount by which the tree's height will
	int iSizeClass, // tree's size class
			iNewSizeClass; // tree's new size class
	//increase
	short int iSpecies = p_oTree->GetSpecies(), iType = p_oTree->GetType();

	//Get the size claa for this tree
	p_oTree->GetValue(p_oPop->GetSizeClassCode(iSpecies, iType), &iSizeClass);
	iNewSizeClass = iSizeClass;
//  fHeight *= 100.0; //transform to cm
//  fNewHeight = fHeight;

	// Get basal area of neighbors:
	p_oTree->GetValue(p_oPop->GetXCode(iSpecies, iType), &fX);
	p_oTree->GetValue(p_oPop->GetYCode(iSpecies, iType), &fY);
	fBAT = GetBAT(fX, fY, p_oPop);

	//Calculate the probability of changing size class

	// effect of height class
	if (iSizeClass == 0) {
		// XS
		fHEffect = mp_fXSm[iSpecies];
	}
	if (iSizeClass == 1) {
		// S
		fHEffect = mp_fSm[iSpecies];
	} else if (iSizeClass == 2) {
		// M
		fHEffect = mp_fMd[iSpecies];
	} else {
		// L
		fHEffect = mp_fLg[iSpecies];
	}

	fTransProb = 1
			/ (1
					+ exp(
							-(mp_fInter[iSpecies] + fHEffect
									+ mp_fJMxC[iSpecies] * fPlotJulMaxTemp / 10
									+ mp_fPC[iSpecies] * fPlotPrecip / 1000
									+ mp_fBA[iSpecies] * fBAT)));

	if (clModelMath::GetRand() < fTransProb)
		iNewSizeClass++;
//  {
//	  // stay in the same size class
//	  if (fHeight < 25){
//		  fNewHeight = 10 + clModelMath::GetRand() * (25 - 10);
//	  }
//	  else if (fHeight < 50){
//		  fNewHeight = 25 + clModelMath::GetRand() * (50 - 25);
//	  }
//	  else {
//		  fNewHeight = 50 + clModelMath::GetRand() * (140 - 50);
//	  }
//  }
//  else
//  {
//	  // transition from one size class to another
//	  if (fHeight < 25){
//		  // growth between 25-h and 50-h
//		  fNewHeight = 25 + clModelMath::GetRand() * (50 - 25);
//	  }
//	  else if (fHeight < 50){
//		  // growth between 50-h and 140-h
//		  fNewHeight = 50 + clModelMath::GetRand() * (140 - 50);
//	  }
//	  else {
//		  // growth between 140-h and 150-h(? how big can the seedling become?)
//		  fNewHeight = 140 + clModelMath::GetRand() * (150 - 140);
//	  }
//  }
//
//  fAmountHeightIncrease = fNewHeight - fHeight;
//  //Transform to m
//  fAmountHeightIncrease /= 100.0;
//
//  return m_fNumberYearsPerTimestep * fAmountHeightIncrease;
	return iNewSizeClass;
}

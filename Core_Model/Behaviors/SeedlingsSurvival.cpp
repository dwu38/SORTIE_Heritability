/*
 * SeedlingsSurvival.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: ucmuser
 */
//---------------------------------------------------------------------------
#include "SeedlingsSurvival.h"
#include "TreePopulation.h"
#include "SimManager.h"
#include "ParsingFunctions.h"
#include "Plot.h"
#include "Allometry.h"
#include <stdio.h>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////////////
clSeedlingsSurvival::clSeedlingsSurvival(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clMortalityBase(
				p_oSimManager) {
	try {
		//Set namestring
		m_sNameString = "Seedlingsmortshell";
		m_sXMLRoot = "SeedlingsSurvival";

		//Null out our pointers
		mp_fInter = NULL;
		mp_fSize0 = NULL;
		mp_fSize1 = NULL;
		mp_fSize2 = NULL;
		mp_fSize3 = NULL;
		mp_fJulyMaxTemp = NULL;
		mp_fPrecip = NULL;
		mp_fBA = NULL;
		mp_fBAt = NULL;
		mp_oPop = NULL;
		m_fYearsPerTimestep = 0;
		mp_iIndexes = NULL;

		m_fRadius = 10;
		m_fMinSaplingHeight = 0;

		//Version 1
		m_fVersionNumber = 1.0;
		m_fMinimumVersionNumber = 1.0;
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clSeedlingsSurvival::clSeedlingsSurvival";
		throw(stcErr);
	}
}

////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////
clSeedlingsSurvival::~clSeedlingsSurvival() {
	delete[] mp_fInter;
	delete[] mp_fSize0;
	delete[] mp_fSize1;
	delete[] mp_fSize2;
	delete[] mp_fSize3;
	delete[] mp_fJulyMaxTemp;
	delete[] mp_fPrecip;
	delete[] mp_fBA;
	delete[] mp_fBAt;
	delete[] mp_iIndexes;
}

////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
////////////////////////////////////////////////////////////////////////////
void clSeedlingsSurvival::DoShellSetup(xercesc::DOMDocument *p_oDoc) {
	floatVal *p_fTempValues = NULL; //for getting species-specific values
	try {
		mp_oPop = (clTreePopulation*) mp_oSimManager->GetPopulationObject(
				"treepopulation");

		DOMElement *p_oElement = GetParentParametersElement(p_oDoc);
		short int iNumSpecies = mp_oPop->GetNumberOfSpecies(), i;

		//Get number of years per timestep
		m_fYearsPerTimestep = mp_oSimManager->GetNumberOfYearsPerTimestep();

		//Declare the arrays we'd like read
		mp_fInter = new float[m_iNumBehaviorSpecies];
		mp_fSize0 = new float[m_iNumBehaviorSpecies];
		mp_fSize1 = new float[m_iNumBehaviorSpecies];
		mp_fSize2 = new float[m_iNumBehaviorSpecies];
		mp_fSize3 = new float[m_iNumBehaviorSpecies];
		mp_fJulyMaxTemp = new float[m_iNumBehaviorSpecies];
		mp_fPrecip = new float[m_iNumBehaviorSpecies];
		mp_fBA = new float[m_iNumBehaviorSpecies];
		mp_fBAt = new float[m_iNumBehaviorSpecies];
		mp_iIndexes = new int[iNumSpecies];

		//Set up the array indexes
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_iIndexes[mp_iWhatSpecies[i]] = i;

		//Declare the species-specific temp array and pre-load with the species that
		//this behavior affects
		p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			p_fTempValues[i].code = mp_iWhatSpecies[i];

		//mortality Intercept
		FillSpeciesSpecificValue(p_oElement, "mo_SeedlingsInter",
				"mo_SeedlingsInterVal", p_fTempValues, m_iNumBehaviorSpecies,
				mp_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fInter[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		// size 0 parameter
		FillSpeciesSpecificValue(p_oElement, "mo_SeedlingsS0",
				"mo_SeedlingsS0Val", p_fTempValues, m_iNumBehaviorSpecies,
				mp_oPop, true);
		//Transfer to the appropriate array buckets
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fSize0[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		// size 1 parameter
		FillSpeciesSpecificValue(p_oElement, "mo_SeedlingsS1",
				"mo_SeedlingsS1Val", p_fTempValues, m_iNumBehaviorSpecies,
				mp_oPop, true);
		//Transfer to the appropriate array buckets
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fSize1[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		// size 2 parameter
		FillSpeciesSpecificValue(p_oElement, "mo_SeedlingsS2",
				"mo_SeedlingsS2Val", p_fTempValues, m_iNumBehaviorSpecies,
				mp_oPop, true);
		//Transfer to the appropriate array buckets
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fSize2[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		// size 3 parameter
		FillSpeciesSpecificValue(p_oElement, "mo_SeedlingsS3",
				"mo_SeedlingsS3Val", p_fTempValues, m_iNumBehaviorSpecies,
				mp_oPop, true);
		//Transfer to the appropriate array buckets
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fSize3[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		// average temperature
		FillSpeciesSpecificValue(p_oElement, "mo_SeedlingsJulyMaxTemp",
				"mo_SeedlingsJMTVal", p_fTempValues, m_iNumBehaviorSpecies,
				mp_oPop, true);
		//Transfer to the appropriate array buckets
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fJulyMaxTemp[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		// water deficit
		FillSpeciesSpecificValue(p_oElement, "mo_SeedlingsPrecip",
				"mo_SeedlingsPVal", p_fTempValues, m_iNumBehaviorSpecies,
				mp_oPop, true);
		//Transfer to the appropriate array buckets
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fPrecip[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		// basal area
		FillSpeciesSpecificValue(p_oElement, "mo_SeedlingsBA",
				"mo_SeedlingsBAVal", p_fTempValues, m_iNumBehaviorSpecies,
				mp_oPop, true);
		//Transfer to the appropriate array buckets
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fBA[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		// basal area threshold
		FillSpeciesSpecificValue(p_oElement, "mo_SeedlingsBAt",
				"mo_SeedlingsBAtVal", p_fTempValues, m_iNumBehaviorSpecies,
				mp_oPop, true);
		//Transfer to the appropriate array buckets
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fBAt[mp_iIndexes[p_fTempValues[i].code]] = p_fTempValues[i].val;

		//Max crowding radius
		FillSingleValue(p_oElement, "mo_SeedlingsRadius", &m_fRadius, true);

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
		stcErr.sFunction = "clSeedlingsSurvival::DoShellSetup";
		throw(stcErr);
	}
}

////////////////////////////////////////////////////////////////////////////
// DoMort
////////////////////////////////////////////////////////////////////////////
deadCode clSeedlingsSurvival::DoMort(clTree *p_oTree, const float &fDbh,
		const short int &iSpecies) {
	clPlot *p_oPlot = mp_oSimManager->GetPlotObject();
	float fX, fY, fSurvivalProb, fMortDiam, fJulyMaxTemp =
			p_oPlot->GetJulMaxTemp(), //Current plot temperature!!
	fPrecip = p_oPlot->GetMeanAnnualPrecip(), //Current plot precipitation
	fBAT, //Neighborhood basal area of trees
			fheightEffect;
	int iTp = p_oTree->GetType(), iSizeClass;
	short int iType = p_oTree->GetType();

	//the tree is a seedling, get the tree's diam10
	if (iTp == clTreePopulation::seedling && 0 == fDbh)
		p_oTree->GetValue(mp_oPop->GetDiam10Code(iSpecies, iTp), &fMortDiam);
	else
		fMortDiam = fDbh;

	// Get tree's location
	p_oTree->GetValue(mp_oPop->GetXCode(iSpecies, iType), &fX);
	p_oTree->GetValue(mp_oPop->GetYCode(iSpecies, iType), &fY);

	//Get the basal area of neighbors
	fBAT = GetBAT(fX, fY, mp_oPop);

	// Get tree's size class and choose the height effect
	p_oTree->GetValue(mp_oPop->GetSizeClassCode(iSpecies, iType), &iSizeClass);
	if (iSizeClass <= 0) {
		fheightEffect = mp_fSize0[iSpecies];
	}
	if (iSizeClass == 1) {
		fheightEffect = mp_fSize1[iSpecies];
	} else if (iSizeClass == 2) {
		fheightEffect = mp_fSize2[iSpecies];
	} else if (iSizeClass == 3) {
		fheightEffect = mp_fSize3[iSpecies];
	}

	//Calculate survival probability for this tree
	fSurvivalProb = mp_fInter[iSpecies] + fheightEffect
			+ mp_fJulyMaxTemp[iSpecies] * fJulyMaxTemp / 10
			+ mp_fPrecip[iSpecies] * fPrecip / 1000 + mp_fBA[iSpecies] * fBAT;
	fSurvivalProb = 1 / (1 + exp(-fSurvivalProb));
	//if the basal area is above threshold, the seedling will automatically die
	if (fBAT > mp_fBAt[iSpecies]) {
		fSurvivalProb = 0;
	}
	//Compound by number of years per timestep
	if (1 != m_fYearsPerTimestep)
		fSurvivalProb = pow(fSurvivalProb, m_fYearsPerTimestep);
	if (clModelMath::GetRand() < fSurvivalProb)
		return notdead;
	else
		return natural;
}

////////////////////////////////////////////////////////////////////////////
// GetBAT
////////////////////////////////////////////////////////////////////////////
float clSeedlingsSurvival::GetBAT(float &fX, float &fY,
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

////////////////////////////////////////////////////////////////////////////
// PreMortCalcs
////////////////////////////////////////////////////////////////////////////
//void clSeedlingsSurvival::PreMortCalcs( clTreePopulation * p_oPop )
//{
//  float fValue;
//  short int i, iX, iY,
//      iNumTotalSpecies = p_oPop->GetNumberOfSpecies(),
//      iNumXCells = mp_oGrid->GetNumberXCells(),
//      iNumYCells = mp_oGrid->GetNumberYCells();
//
//  //Reset the values in the survival grid to -1
//  fValue = -1;
//  for (iX = 0; iX < iNumXCells; iX++)
//    for (iY = 0; iY < iNumYCells; iY++) {
//      mp_oGrid->SetValueOfCell(iX, iY, m_iBATCode, fValue);
//      for (i = 0; i < iNumTotalSpecies; i++)
//          mp_oGrid->SetValueOfCell(iX, iY, mp_iGridSurvivalCodes[i], fValue);
//    }
//
//}

/////////////////////////////////////////////////////////////////////////////
// SetupGrid()
/////////////////////////////////////////////////////////////////////////////
//void clSeedlingsSurvival::SetupGrid(clTreePopulation *p_oPop)
//{
//  std::stringstream sLabel;
//  short int iNumTotalSpecies = p_oPop->GetNumberOfSpecies(),
//            i; //loop counter
//
//  //Declare the arrays for our grid codes
//  mp_iGridSurvivalCodes = new short int[iNumTotalSpecies];
//
//  //Check to see if this grid has already been created
//  mp_oGrid = mp_oSimManager->GetGridObject("Seedlings Survival");
//
//  if (NULL == mp_oGrid) {
//
//    //Create the grid with five float data members for each species
//    mp_oGrid = mp_oSimManager->CreateGrid( "Seedlings Survival",
//                          0,                    //number of ints
//                          iNumTotalSpecies + 1, //number of floats
//                          0,                    //number of chars
//                          0);                   //number of bools
//
//    //Register the data member called "survival_x"
//    for (i = 0; i < iNumTotalSpecies; i++) {
//      sLabel << "survival_" << i;
//      mp_iGridSurvivalCodes[i] = mp_oGrid->RegisterFloat(sLabel.str());
//      sLabel.str("");
//    }
//
//    //Register the data member called "BAT"
//    m_iBATCode = mp_oGrid->RegisterFloat("BAT");
//  }
//  else {
//    //Grid already exists - get the codes
//    //Get the data member called "survival_x"
//    for (i = 0; i < iNumTotalSpecies; i++) {
//      sLabel << "survival_" << i;
//      mp_iGridSurvivalCodes[i] = mp_oGrid->GetFloatDataCode(sLabel.str());
//      if (-1 == mp_iGridSurvivalCodes[i]) {
//        modelErr stcErr;
//        stcErr.sFunction = "clSeedlingsSurvival::SetupGrid";
//        std::stringstream s;
//        s << "Couldn't find the \"" << sLabel.str()
//          << "\" member of the \"Seedlings Survival\" grid.";
//        stcErr.iErrorCode = BAD_DATA;
//        throw( stcErr );
//      }
//
//      m_iBATCode = mp_oGrid->GetFloatDataCode( "BAT" );
//      if (-1 == m_iBATCode) {
//        modelErr stcErr;
//        stcErr.sFunction = "clSeedlingsSurvival::SetupGrid" ;
//        stcErr.sMoreInfo = "Couldn't find the \"BAT\" member of the \"Seedlings Survival\" grid.";
//        stcErr.iErrorCode = BAD_DATA;
//        throw( stcErr );
//      }
//    }
//  }
//
//}
//


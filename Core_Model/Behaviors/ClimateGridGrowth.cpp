/*
 * ClimateGridGrowth.cpp
 *
 *  Created on: Mar 14, 2017
 *      Author: ucmuser
 */

#include "ClimateGridGrowth.h"
#include "TreePopulation.h"
#include "SimManager.h"
#include "ParsingFunctions.h"
#include "Grid.h"
#include "Plot.h"
#include "ModelMath.h"
#include "GrowthOrg.h"
#include <sstream>

//---------------------------------------------------------------------------
// DoubleMMRelGrowth.cpp
//---------------------------------------------------------------------------
//#include "DoubleMMRelGrowth.h"
//#include "GrowthOrg.h"
//#include "Tree.h"
//#include "TreePopulation.h"
//#include "SimManager.h"
//#include "Grid.h"
//#include "ParsingFunctions.h"
//#include <math.h>
//#include <stdio.h>
//#include <sstream>

//////////////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////////////
clClimateGridGrowth::clClimateGridGrowth(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clGrowthBase(
				p_oSimManager) {

	try {
		mp_fDBHSlope = NULL;
		mp_fTempSlope = NULL;
		mp_fPrecipSlope = NULL;
		mp_fCdIntercept = NULL;
		mp_iIndexes = NULL;
		m_fYearsPerTimestep = 0;

		m_sNameString = "ClimateGridGrowthshell";
		m_sXMLRoot = "ClimateGridGrowth";

		mp_oPrecipGrid = NULL;
		mp_oTemperatureGrid = NULL;
		m_iPrecipCode = -1;
		m_iTemperatureCode = -1;
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clClimateGridGrowth::clClimateGridGrowth";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////////////
clClimateGridGrowth::~clClimateGridGrowth() {
//  delete[] mp_fResourceInfluence;
	delete[] mp_fDBHSlope;
	delete[] mp_fTempSlope;
	delete[] mp_fPrecipSlope;
	delete[] mp_fCdIntercept;
	delete[] mp_iIndexes;
}

//////////////////////////////////////////////////////////////////////////////
// SetNameData()
//////////////////////////////////////////////////////////////////////////////
void clClimateGridGrowth::SetNameData(std::string sNameString) {
	try {

		//Check the string passed and set the flags accordingly
		if (sNameString.compare("ClimateGridGrowth") == 0) {
			m_iGrowthMethod = diameter_auto;
		} else if (sNameString.compare("ClimateGridGrowth diam only") == 0) {
			m_iGrowthMethod = diameter_only;
		}

		else {
			modelErr stcErr;
			stcErr.iErrorCode = BAD_DATA;
			std::stringstream s;
			s << "Unrecognized behavior name \"" << sNameString << "\".";
			stcErr.sFunction = "clClimateGridGrowth::SetNameData";
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
		stcErr.sFunction = "clClimateGridGrowth::SetNameData";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
//////////////////////////////////////////////////////////////////////////////
void clClimateGridGrowth::DoShellSetup(DOMDocument *p_oDoc) {
	floatVal *p_fTempValues = NULL;
	try {
		clTreePopulation *p_oPop =
				(clTreePopulation*) mp_oSimManager->GetPopulationObject(
						"treepopulation");
		DOMElement *p_oElement = GetParentParametersElement(p_oDoc);

//	    short int iNumSpecies = mp_oGrowthOrg->GetNumberOfSpecies(), i;
		short int iNumSpecies = p_oPop->GetNumberOfSpecies(), i;

		//Get number of years per timestep
		m_fYearsPerTimestep = mp_oSimManager->GetNumberOfYearsPerTimestep();

		//Declare the arrays we'd like read
		mp_fDBHSlope = new float[m_iNumBehaviorSpecies];
		mp_fTempSlope = new float[m_iNumBehaviorSpecies];
		mp_fPrecipSlope = new float[m_iNumBehaviorSpecies];
		mp_fCdIntercept = new float[m_iNumBehaviorSpecies];
		mp_iIndexes = new int[iNumSpecies];

		//Set up the array indexes
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_iIndexes[mp_iWhatSpecies[i]] = i;

		//Declare the species-specific temp array and pre-load with the species that
		//this behavior affects
		p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			p_fTempValues[i].code = mp_iWhatSpecies[i];

		//Read the base variables with the base class function
		// GetParameterFileData( p_oDoc );

		//growth intercept
		FillSpeciesSpecificValue(p_oElement, "gr_ClimateGridIntercept",
				"gr_CdiVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fCdIntercept[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		//DBH growth slope
		FillSpeciesSpecificValue(p_oElement, "gr_ClimateGridDBHSlope",
				"gr_CdDBHsVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fDBHSlope[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		//Temperature growth slope
		FillSpeciesSpecificValue(p_oElement, "gr_ClimateGridTempSlope",
				"gr_CdTempsVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fTempSlope[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		//Precipitation growth slope
		FillSpeciesSpecificValue(p_oElement, "gr_ClimateGridPrecipSlope",
				"gr_CdPrecipsVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fPrecipSlope[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		//Make sure all species/type combos have "Light" registered
		// no need for light in this behaviour
//    for ( i = 0; i < m_iNumSpeciesTypeCombos; i++ )
//    {
//      if ( -1 == mp_oGrowthOrg->GetLightCode( mp_whatSpeciesTypeCombos[i].iSpecies, mp_whatSpeciesTypeCombos[i].iType ) )
//      {
//        modelErr stcErr;
//        stcErr.sFunction = "clDoubleMMRelGrowth::DoShellSetup" ;
//        std::stringstream s;
//        s << "Type/species combo species="
//          << mp_whatSpeciesTypeCombos[i].iSpecies << " type="
//          << mp_whatSpeciesTypeCombos[i].iType
//          << " does not have a required light behavior.";
//        stcErr.sMoreInfo = s.str();
//        stcErr.iErrorCode = BAD_DATA;
//        throw( stcErr );
//      }
//    }

		//Get "Precip" grid
		mp_oPrecipGrid = mp_oSimManager->GetGridObject("Precip");
		if ( NULL == mp_oPrecipGrid) {
			modelErr stcErr;
			stcErr.sFunction = "clClimateGridGrowth::DoShellSetup";
			stcErr.sMoreInfo = "Can't find required grid object \"Precip\".";
			stcErr.iErrorCode = CANT_FIND_OBJECT;
			throw(stcErr);
		}

		m_iPrecipCode = mp_oPrecipGrid->GetFloatDataCode("Precip");
		if (-1 == m_iPrecipCode) {
			modelErr stcErr;
			stcErr.sFunction = "clClimateGridGrowth::DoShellSetup";
			stcErr.sMoreInfo = "Grid object \"Precip\" is set up incorrectly.";
			stcErr.iErrorCode = CANT_FIND_OBJECT;
			throw(stcErr);
		}

		//Get Temperature grid
		mp_oTemperatureGrid = mp_oSimManager->GetGridObject("Temperature");
		if ( NULL == mp_oTemperatureGrid) {
			modelErr stcErr;
			stcErr.sFunction = "clClimateGridGrowth::DoShellSetup";
			stcErr.sMoreInfo =
					"Can't find required grid object \"Temperature\".";
			stcErr.iErrorCode = CANT_FIND_OBJECT;
			throw(stcErr);
		}

		m_iTemperatureCode = mp_oTemperatureGrid->GetFloatDataCode(
				"Temperature");
		if (-1 == m_iTemperatureCode) {
			modelErr stcErr;
			stcErr.sFunction = "clClimateGridGrowth::DoShellSetup";
			stcErr.sMoreInfo =
					"Grid object \"Temperature\" is set up incorrectly.";
			stcErr.iErrorCode = CANT_FIND_OBJECT;
			throw(stcErr);
		}

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
		stcErr.sFunction = "clClimateGridGrowth::DoShellSetup";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// CalcDiameterGrowthValue()
//////////////////////////////////////////////////////////////////////////////
float clClimateGridGrowth::CalcDiameterGrowthValue(clTree *p_oTree,
		clTreePopulation *p_oPop, float fHeightGrowth) {

	float result, //amount of annual relative growth
			fDiam, //tree's diameter value
			fPrecip, //amount of Precipitation
			fTemperature, //amount of Temperature
			fX, fY; //tree's coordinates

	short int iSpecies = p_oTree->GetSpecies(), iType = p_oTree->GetType();
//      iLightCode = mp_oGrowthOrg->GetLightCode( iSpecies, iType );

	//***************************************
	// Calculate the amount of growth
	//***************************************

	//Get the gli value for this tree
//  p_oTree->GetValue( iLightCode, & fGli );

	//Get the proper resource level from the tree's coordinates
	p_oTree->GetValue(p_oPop->GetXCode(iSpecies, iType), &fX);
	p_oTree->GetValue(p_oPop->GetYCode(iSpecies, iType), &fY);
	mp_oPrecipGrid->GetValueAtPoint(fX, fY, m_iPrecipCode, &fPrecip);
	mp_oTemperatureGrid->GetValueAtPoint(fX, fY, m_iTemperatureCode,
			&fTemperature);

	//Get the appropriate diameter for this tree
	p_oTree->GetValue(mp_oGrowthOrg->GetDiamCode(iSpecies, iType), &fDiam);

	//Get the value of the Michaelis-Menton function
//  fAnnualRelativeGrowth = mp_fAsympDiamGrowth[iSpecies] + (mp_fResourceInfluence[iSpecies] * fResource);
//  fAnnualRelativeGrowth = ( fAnnualRelativeGrowth * fGli ) / ( (fAnnualRelativeGrowth / mp_fSlopeDiamGrowthResponse[iSpecies]) + fGli );

	//Compound the relative growth over the number of years/time step
//  fCompoundRelativeGrowth = pow( 1.0 + fAnnualRelativeGrowth, m_fNumberYearsPerTimestep );

	//Calculate amount of diameter increase based on the tree's diameter
//  fAmountDiamIncrease = fDiam * ( fCompoundRelativeGrowth - 1.0 );

	result = fPrecip * mp_fPrecipSlope[mp_iIndexes[iSpecies]]
			+ fTemperature * mp_fTempSlope[mp_iIndexes[iSpecies]]
			+ fDiam * mp_fDBHSlope[mp_iIndexes[iSpecies]]
			+ mp_fCdIntercept[mp_iIndexes[iSpecies]];

	return m_fYearsPerTimestep * result;
}

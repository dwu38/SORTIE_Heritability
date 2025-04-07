//---------------------------------------------------------------------------
// ClimateDepGrowth.cpp
//---------------------------------------------------------------------------
#include "ClimateDepGrowth.h"
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
clClimateDepGrowth::clClimateDepGrowth(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clGrowthBase(
				p_oSimManager) {

	try {
		mp_fDBHSlope = NULL;
		mp_fTempSlope = NULL;
		mp_fPrecipSlope = NULL;
		mp_fCdIntercept = NULL;
		mp_iIndexes = NULL;
		m_fYearsPerTimestep = 0;

		m_sNameString = "ClimateDepgrowthshell";
		m_sXMLRoot = "ClimateDepGrowth";
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clClimateDepGrowth::clClimateDepGrowth";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////////////////////////////////////*/
clClimateDepGrowth::~clClimateDepGrowth() {
	delete[] mp_fDBHSlope;
	delete[] mp_fTempSlope;
	delete[] mp_fPrecipSlope;
	delete[] mp_fCdIntercept;
	delete[] mp_iIndexes;
}

//////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
/////////////////////////////////////////////////////////////////////////////*/
void clClimateDepGrowth::DoShellSetup(DOMDocument *p_oDoc) {
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

		//growth intercept
		FillSpeciesSpecificValue(p_oElement, "gr_ClimateDepIntercept",
				"gr_CdiVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fCdIntercept[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		//DBH growth slope
		FillSpeciesSpecificValue(p_oElement, "gr_ClimateDepDBHSlope",
				"gr_CdDBHsVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fDBHSlope[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		//Temperature growth slope
		FillSpeciesSpecificValue(p_oElement, "gr_ClimateDepTempSlope",
				"gr_CdTempsVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fTempSlope[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		//Precipitation growth slope
		FillSpeciesSpecificValue(p_oElement, "gr_ClimateDepPrecipSlope",
				"gr_CdPrecipsVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
				true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fPrecipSlope[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

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
		stcErr.sFunction = "clClimateDepGrowth::DoShellSetup";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// SetNameData()
/////////////////////////////////////////////////////////////////////////////*/
void clClimateDepGrowth::SetNameData(std::string sNameString) {
	try {

		//Check the string passed and set the flags accordingly
		if (sNameString.compare("ClimateDepGrowth") == 0) {
			m_iGrowthMethod = diameter_auto;
		} else if (sNameString.compare("ClimateDepGrowth diam only") == 0) {
			m_iGrowthMethod = diameter_only;
		} else {
			modelErr stcErr;
			stcErr.iErrorCode = BAD_DATA;
			std::stringstream s;
			s << "Unrecognized behavior name \"" << sNameString << "\".";
			stcErr.sFunction = "clClimateDepGrowth::SetNameData";
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
		stcErr.sFunction = "clClimateDepGrowth::SetNameData";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// CalcDiameterGrowthValue()
/////////////////////////////////////////////////////////////////////////////*/
float clClimateDepGrowth::CalcDiameterGrowthValue(clTree *p_oTree,
		clTreePopulation *p_oPop, float fHeightGrowth) {

	clPlot *p_oPlot = mp_oSimManager->GetPlotObject();
	float fDiam, result, fPlotTemp = p_oPlot->GetMeanAnnualTemp(), //Current plot temperature
	fPlotPrecip = p_oPlot->GetMeanAnnualPrecip(); //Current plot precipitation;
	int iSp = p_oTree->GetSpecies(), iTp = p_oTree->GetType();

	//Get the tree's diameter
	p_oTree->GetValue(mp_oGrowthOrg->GetDiamCode(iSp, iTp), &fDiam);

	//Calculate the function value
	result = fPlotPrecip * mp_fPrecipSlope[mp_iIndexes[iSp]]
			+ fPlotTemp * mp_fTempSlope[mp_iIndexes[iSp]]
			+ fDiam * mp_fDBHSlope[mp_iIndexes[iSp]]
			+ mp_fCdIntercept[mp_iIndexes[iSp]];
	return m_fYearsPerTimestep * result;

}


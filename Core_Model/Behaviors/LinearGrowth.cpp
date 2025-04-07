//---------------------------------------------------------------------------
// LinearGrowth.cpp
//---------------------------------------------------------------------------
#include "LinearGrowth.h"
#include "TreePopulation.h"
#include "SimManager.h"
#include "ParsingFunctions.h"
#include "Grid.h"
#include "ModelMath.h"
#include "GrowthOrg.h"
#include <sstream>

//////////////////////////////////////////////////////////////////////////////
// Constructor
/////////////////////////////////////////////////////////////////////////////*/
clLinearGrowth::clLinearGrowth(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clGrowthBase(
				p_oSimManager) {

	try {
		mp_fSlope = NULL;
		mp_fIntercept = NULL;
		mp_iIndexes = NULL;
		m_fYearsPerTimestep = 0;

		m_sNameString = "lineargrowthshell";
		m_sXMLRoot = "LinearGrowth";
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clLinearGrowth::clLinearGrowth";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////////////////////////////////////*/
clLinearGrowth::~clLinearGrowth() {
	delete[] mp_fSlope;
	delete[] mp_fIntercept;
	delete[] mp_iIndexes;
}

//////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
/////////////////////////////////////////////////////////////////////////////*/
void clLinearGrowth::DoShellSetup(DOMDocument *p_oDoc) {
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
		mp_fSlope = new float[m_iNumBehaviorSpecies];
		mp_fIntercept = new float[m_iNumBehaviorSpecies];
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
		FillSpeciesSpecificValue(p_oElement, "gr_linearIntercept", "gr_liVal",
				p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fIntercept[mp_iIndexes[p_fTempValues[i].code]] =
					p_fTempValues[i].val;

		//growth slope
		FillSpeciesSpecificValue(p_oElement, "gr_linearSlope", "gr_lsVal",
				p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);
		//Transfer values to our permanent array
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fSlope[mp_iIndexes[p_fTempValues[i].code]] =
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
		stcErr.sFunction = "clLinearGrowth::DoShellSetup";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// SetNameData()
/////////////////////////////////////////////////////////////////////////////*/
void clLinearGrowth::SetNameData(std::string sNameString) {
	try {

		//Check the string passed and set the flags accordingly
		if (sNameString.compare("LinearGrowth") == 0) {
			m_iGrowthMethod = diameter_auto;
		} else if (sNameString.compare("LinearGrowth diam only") == 0) {
			m_iGrowthMethod = diameter_only;
		} else {
			modelErr stcErr;
			stcErr.iErrorCode = BAD_DATA;
			std::stringstream s;
			s << "Unrecognized behavior name \"" << sNameString << "\".";
			stcErr.sFunction = "clLinearGrowth::SetNameData";
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
		stcErr.sFunction = "clLinearGrowth::SetNameData";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// CalcDiameterGrowthValue()
/////////////////////////////////////////////////////////////////////////////*/
float clLinearGrowth::CalcDiameterGrowthValue(clTree *p_oTree,
		clTreePopulation *p_oPop, float fHeightGrowth) {

	float fDiam;
	int iSp = p_oTree->GetSpecies(), iTp = p_oTree->GetType();

	//Get the tree's diameter
	p_oTree->GetValue(mp_oGrowthOrg->GetDiamCode(iSp, iTp), &fDiam);

	//Calculate the function value
	return m_fYearsPerTimestep
			* clModelMath::CalcPointValue(fDiam, mp_fSlope[mp_iIndexes[iSp]],
					mp_fIntercept[mp_iIndexes[iSp]]);
}

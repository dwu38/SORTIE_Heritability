//---------------------------------------------------------------------------
// LognormalGrowth.cpp
//---------------------------------------------------------------------------
#include "LognormalGrowth.h"
#include "TreePopulation.h"
#include "SimManager.h"
#include "ParsingFunctions.h"
#include "GrowthOrg.h"
#include <math.h>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////
// Constructor
/////////////////////////////////////////////////////////////////////////////*/
clLognormalGrowth::clLognormalGrowth(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clGrowthBase(
				p_oSimManager) {

	try {
		mp_fIncAtDiam36 = NULL;
		mp_fShape = NULL;
		mp_fShade = NULL;
		mp_iIndexes = NULL;

		m_sNameString = "lognormalgrowthshell";
		m_sXMLRoot = "LognormalGrowth";

		m_fConversionFactor = 0;
		m_fNumberOfYearsPerTimestep = 0;
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clLognormalGrowth::clLognormalGrowth";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////////////////////////////////////*/
clLognormalGrowth::~clLognormalGrowth() {
	delete[] mp_fIncAtDiam36;
	delete[] mp_fShape;
	delete[] mp_fShade;
	delete[] mp_iIndexes;
}

//////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
/////////////////////////////////////////////////////////////////////////////*/
void clLognormalGrowth::DoShellSetup(DOMDocument *p_oDoc) {
	try {
		clPopulationBase *p_oTemp = mp_oSimManager->GetPopulationObject(
				"treepopulation");
		clTreePopulation *p_oPop = (clTreePopulation*) p_oTemp;
		DOMElement *p_oElement = p_oDoc->getDocumentElement();
		floatVal *p_fTempValues; //for getting species-specific values
		short int iNumSpecies = p_oPop->GetNumberOfSpecies(), i;

		//Conversion factor
		if (height_only == m_iGrowthMethod) {
			//The conversion factor needs to convert from cm to m
			m_fConversionFactor = 0.01;
		} else {
			//The conversion factor needs to convert from mm to cm and from radial to diameter
			m_fConversionFactor = 0.2;
		}

		m_fNumberOfYearsPerTimestep =
				mp_oSimManager->GetNumberOfYearsPerTimestep();

		//Make the list of indexes
		mp_iIndexes = new short int[iNumSpecies];
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_iIndexes[mp_iWhatSpecies[i]] = i;

		//Declare the arrays we'd like read
		mp_fIncAtDiam36 = new float[m_iNumBehaviorSpecies];
		mp_fShape = new float[m_iNumBehaviorSpecies];
		mp_fShade = new float[m_iNumBehaviorSpecies];

		//Declare the species-specific temp array and pre-load with the species that
		//this behavior affects
		p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			p_fTempValues[i].code = mp_iWhatSpecies[i];

		if (height_only == m_iGrowthMethod) {
			//Growth increment at diameter 36
			FillSpeciesSpecificValue(p_oElement,
					"gr_heightLognormalIncAtDiam36", "gr_hliad36Val",
					p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);

			//Transfer values to our permanent array
			for (i = 0; i < m_iNumBehaviorSpecies; i++)
				mp_fIncAtDiam36[mp_iIndexes[p_fTempValues[i].code]] =
						p_fTempValues[i].val;

			//Shape parameter
			FillSpeciesSpecificValue(p_oElement, "gr_heightLognormalShapeParam",
					"gr_hlspVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
					true);

			//Transfer values to our permanent array
			for (i = 0; i < m_iNumBehaviorSpecies; i++)
				mp_fShape[mp_iIndexes[p_fTempValues[i].code]] =
						p_fTempValues[i].val;

			//Effects of shade parameter
			FillSpeciesSpecificValue(p_oElement,
					"gr_heightLognormalEffectOfShade", "gr_hleosVal",
					p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);

			//Transfer values to our permanent array
			for (i = 0; i < m_iNumBehaviorSpecies; i++)
				mp_fShade[mp_iIndexes[p_fTempValues[i].code]] =
						p_fTempValues[i].val;
		} else {
			//Growth increment at diameter 36
			FillSpeciesSpecificValue(p_oElement, "gr_diamLognormalIncAtDiam36",
					"gr_dliad36Val", p_fTempValues, m_iNumBehaviorSpecies,
					p_oPop, true);

			//Transfer values to our permanent array
			for (i = 0; i < m_iNumBehaviorSpecies; i++)
				mp_fIncAtDiam36[mp_iIndexes[p_fTempValues[i].code]] =
						p_fTempValues[i].val;

			//Shape parameter
			FillSpeciesSpecificValue(p_oElement, "gr_diamLognormalShapeParam",
					"gr_dlspVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
					true);

			//Transfer values to our permanent array
			for (i = 0; i < m_iNumBehaviorSpecies; i++)
				mp_fShape[mp_iIndexes[p_fTempValues[i].code]] =
						p_fTempValues[i].val;

			//Effects of shade parameter
			FillSpeciesSpecificValue(p_oElement,
					"gr_diamLognormalEffectOfShade", "gr_dleosVal",
					p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);

			//Transfer values to our permanent array
			for (i = 0; i < m_iNumBehaviorSpecies; i++)
				mp_fShade[mp_iIndexes[p_fTempValues[i].code]] =
						p_fTempValues[i].val;
		}

		//Make sure all species/type combos have "Light" registered
		for (i = 0; i < m_iNumSpeciesTypeCombos; i++) {
			if (-1
					== mp_oGrowthOrg->GetLightCode(
							mp_whatSpeciesTypeCombos[i].iSpecies,
							mp_whatSpeciesTypeCombos[i].iType)) {
				modelErr stcErr;
				stcErr.sFunction = "clLognormalGrowth::DoShellSetup";
				std::stringstream s;
				s << "Type/species combo species="
						<< mp_whatSpeciesTypeCombos[i].iSpecies << " type="
						<< mp_whatSpeciesTypeCombos[i].iType
						<< " does not have a required light behavior.";
				stcErr.sMoreInfo = s.str();
				stcErr.iErrorCode = BAD_DATA;
				throw(stcErr);

			}
		}

		delete[] p_fTempValues;
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clLognormalGrowth::DoShellSetup";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// SetNameData()
/////////////////////////////////////////////////////////////////////////////*/
void clLognormalGrowth::SetNameData(std::string sNameString) {
	try {

		//Check the string passed and set the flags accordingly
		if (sNameString.compare("LognormalGrowth") == 0) {
			m_iGrowthMethod = diameter_auto;
		} else if (sNameString.compare("LognormalGrowth diam only") == 0) {
			m_iGrowthMethod = diameter_only;
		} else if (sNameString.compare("LognormalGrowth height only") == 0) {
			m_iGrowthMethod = height_only;
			m_bGoLast = true;
		} else {
			modelErr stcErr;
			stcErr.iErrorCode = BAD_DATA;
			std::stringstream s;
			s << "Unrecognized behavior name \"" << sNameString << "\".";
			stcErr.sFunction = "clLognormalGrowth::SetNameData";
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
		stcErr.sFunction = "clLognormalGrowth::SetNameData";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// CalculateFunctionValue()
/////////////////////////////////////////////////////////////////////////////*/
float clLognormalGrowth::CalculateFunctionValue(int iSpecies, float fGLI,
		float fDiam) {
	short int iIndex = mp_iIndexes[iSpecies];

	//Calculate the function value
	return mp_fIncAtDiam36[iIndex]
			* exp(-0.5 * pow(log(fDiam / 36) / mp_fShape[iIndex], 2))
			* pow(fGLI / 100, mp_fShade[iIndex]);
}

//////////////////////////////////////////////////////////////////////////////
// CalcHeightGrowthValue()
/////////////////////////////////////////////////////////////////////////////*/
float clLognormalGrowth::CalcHeightGrowthValue(clTree *p_oTree,
		clTreePopulation *p_oPop, float fDiameterGrowth) {

	float fGLI, fDiam, fDiamIncrement, fHeightIncrement = 0;
	int iSpecies = p_oTree->GetSpecies(), iType = p_oTree->GetType(),
			iLightCode = mp_oGrowthOrg->GetLightCode(iSpecies, iType),
			iDiamCode = mp_oGrowthOrg->GetDiamCode(iSpecies, iType), i;

	//Get the GLI and the diameter
	p_oTree->GetValue(iLightCode, &fGLI);
	p_oTree->GetValue(iDiamCode, &fDiam);

	//Calculate the diameter increment
	fDiamIncrement = fDiameterGrowth / m_fNumberOfYearsPerTimestep;

	//Increment the height
	for (i = 0; i < m_fNumberOfYearsPerTimestep; i++) {
		fHeightIncrement += CalculateFunctionValue(iSpecies, fGLI, fDiam);
		fDiam += fDiamIncrement;
	}

	return fHeightIncrement * m_fConversionFactor;
}

//////////////////////////////////////////////////////////////////////////////
// CalcDiameterGrowthValue()
/////////////////////////////////////////////////////////////////////////////*/
float clLognormalGrowth::CalcDiameterGrowthValue(clTree *p_oTree,
		clTreePopulation *p_oPop, float fHeightGrowth) {
	float fGLI, fDiam, fDiamIncrement = 0, fTotalDiamIncrement = 0;
	int iSpecies = p_oTree->GetSpecies(), iType = p_oTree->GetType(),
			iLightCode = mp_oGrowthOrg->GetLightCode(iSpecies, iType),
			iDiamCode = mp_oGrowthOrg->GetDiamCode(iSpecies, iType), i;

	//Get the GLI and the diameter
	p_oTree->GetValue(iLightCode, &fGLI);
	p_oTree->GetValue(iDiamCode, &fDiam);

	//Increment the height
	for (i = 0; i < m_fNumberOfYearsPerTimestep; i++) {
		fDiamIncrement = CalculateFunctionValue(iSpecies, fGLI, fDiam)
				* m_fConversionFactor;
		fTotalDiamIncrement += fDiamIncrement;
		fDiam += fDiamIncrement;
	}

	return fTotalDiamIncrement;
}

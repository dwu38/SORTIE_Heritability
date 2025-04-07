/*
 * JanMinPrecipGrowth.cpp
 *
 *  Created on: May 14, 2018
 *      Author: ucmuser
 */

//---------------------------------------------------------------------------
#include <JanMinPrecipGrowth.h>
#include "TreePopulation.h"
#include "SimManager.h"
#include "ParsingFunctions.h"
#include "Plot.h"
#include "GrowthOrg.h"
#include "Allometry.h"

#include "NCI/NCITermBase.h"

#include <stdio.h>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////////////
clJanMinPrecipGrowth::clJanMinPrecipGrowth(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager), clGrowthBase(
				p_oSimManager), clNCIBehaviorBase() {

	//Set namestring
	m_sNameString = "JanMinPrecipgrowthshell";
	m_sXMLRoot = "JanMinPrecipGrowth";

	//Version 1
	m_fVersionNumber = 1.0;
	m_fMinimumVersionNumber = 1.0;

	//Null out our pointers
	mp_iGrowthCodes = NULL;
	mp_fIntercept = NULL;
	mp_fDBHSlope = NULL;
	mp_fNCISlope = NULL;
	mp_fJanMinSlope = NULL;
	mp_fPrecipSlope = NULL;

	mp_fRandParameter = NULL;

	m_iStochasticGrowthMethod = deterministic_pdf;
	Adjust = NULL;

	m_iNumTotalSpecies = 0;

}

////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////
clJanMinPrecipGrowth::~clJanMinPrecipGrowth() {

	int i;
	if (mp_iGrowthCodes) {
		for (i = 0; i < m_iNumTotalSpecies; i++) {
			delete[] mp_iGrowthCodes[i];
		}
		delete[] mp_iGrowthCodes;
	}
	delete[] mp_fIntercept;
	delete[] mp_fDBHSlope;
	delete[] mp_fNCISlope;
	delete[] mp_fJanMinSlope;
	delete[] mp_fPrecipSlope;
	delete[] mp_fRandParameter;
}

////////////////////////////////////////////////////////////////////////////
// GetTreeMemberCodes()
////////////////////////////////////////////////////////////////////////////
void clJanMinPrecipGrowth::GetTreeMemberCodes() {
	clTreePopulation *p_oPop =
			(clTreePopulation*) mp_oSimManager->GetPopulationObject(
					"treepopulation");
	int iNumTypes = p_oPop->GetNumberOfTypes(), i, j;

	//Get codes for growth
	mp_iGrowthCodes = new short int*[m_iNumTotalSpecies];
	for (i = 0; i < m_iNumTotalSpecies; i++) {
		mp_iGrowthCodes[i] = new short int[iNumTypes];
		for (j = 0; j < iNumTypes; j++) {
			mp_iGrowthCodes[i][j] = -1;
		}
	}

	for (i = 0; i < m_iNumSpeciesTypeCombos; i++) {
		//Get the code from growth org
		mp_iGrowthCodes[mp_whatSpeciesTypeCombos[i].iSpecies][mp_whatSpeciesTypeCombos[i].iType] =
				mp_oGrowthOrg->GetGrowthCode(
						mp_whatSpeciesTypeCombos[i].iSpecies,
						mp_whatSpeciesTypeCombos[i].iType);
	}
}

////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
////////////////////////////////////////////////////////////////////////////
void clJanMinPrecipGrowth::DoShellSetup(xercesc::DOMDocument *p_oDoc) {
	clTreePopulation *p_oPop =
			(clTreePopulation*) mp_oSimManager->GetPopulationObject(
					"treepopulation");

	m_iNumTotalSpecies = p_oPop->GetNumberOfSpecies();

	DOMElement *p_oElement = GetParentParametersElement(p_oDoc);
	floatVal *p_fTempValues; //for getting species-specific values
	int iTemp;
	short int i; //loop counters

	//Set up our floatVal array that will extract values only for the species
	//assigned to this behavior
	p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
	for (i = 0; i < m_iNumBehaviorSpecies; i++)
		p_fTempValues[i].code = mp_iWhatSpecies[i];

	//Intercept
	mp_fIntercept = new float[m_iNumTotalSpecies];
	FillSpeciesSpecificValue(p_oElement, "gr_JanMinPrecipIntercept",
			"gr_JanMinPrecipiVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
			true);
	//Transfer to the appropriate array buckets
	for (i = 0; i < m_iNumBehaviorSpecies; i++)
		mp_fIntercept[p_fTempValues[i].code] = p_fTempValues[i].val;

	//DBH slope
	mp_fDBHSlope = new float[m_iNumTotalSpecies];
	FillSpeciesSpecificValue(p_oElement, "gr_JanMinPrecipDBHSlope",
			"gr_JanMinPrecipdVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
			true);
	//Transfer to the appropriate array buckets
	for (i = 0; i < m_iNumBehaviorSpecies; i++)
		mp_fDBHSlope[p_fTempValues[i].code] = p_fTempValues[i].val;

	// NCI slope
	mp_fNCISlope = new float[m_iNumTotalSpecies];
	FillSpeciesSpecificValue(p_oElement, "gr_JanMinPrecipNCISlope",
			"gr_JanMinPrecipnVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
			true);
	//Transfer to the appropriate array buckets
	for (i = 0; i < m_iNumBehaviorSpecies; i++)
		mp_fNCISlope[p_fTempValues[i].code] = p_fTempValues[i].val;

	// Temperature slope
	mp_fJanMinSlope = new float[m_iNumTotalSpecies];
	FillSpeciesSpecificValue(p_oElement, "gr_JanMinPrecipJanMinSlope",
			"gr_JanMinPrecipjmnVal", p_fTempValues, m_iNumBehaviorSpecies,
			p_oPop, true);
	//Transfer to the appropriate array buckets
	for (i = 0; i < m_iNumBehaviorSpecies; i++)
		mp_fJanMinSlope[p_fTempValues[i].code] = p_fTempValues[i].val;

	// Water availability slope
	mp_fPrecipSlope = new float[m_iNumTotalSpecies];
	FillSpeciesSpecificValue(p_oElement, "gr_JanMinPrecipPrecipSlope",
			"gr_JanMinPrecippVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop,
			true);
	//Transfer to the appropriate array buckets
	for (i = 0; i < m_iNumBehaviorSpecies; i++)
		mp_fPrecipSlope[p_fTempValues[i].code] = p_fTempValues[i].val;

	//Get the type of growth stochasticity desired
	FillSingleValue(p_oElement, "gr_stochGrowthMethod", &iTemp, true);

	//Make sure the value is valid
	if (deterministic_pdf == iTemp) {

		m_iStochasticGrowthMethod = deterministic_pdf;
		Adjust = &clJanMinPrecipGrowth::DeterministicAdjust;

	} else if (lognormal_pdf == iTemp) {

		m_iStochasticGrowthMethod = lognormal_pdf;
		Adjust = &clJanMinPrecipGrowth::LognormalAdjust;

	} else if (normal_pdf == iTemp) {

		m_iStochasticGrowthMethod = normal_pdf;
		Adjust = &clJanMinPrecipGrowth::NormalAdjust;

	} else {
		modelErr stcErr;
		stcErr.sFunction = "clJanMinPrecipGrowth::DoShellSetup";
		std::stringstream s;
		s << "Unrecognized value for gr_stochGrowthMethod: " << iTemp;
		stcErr.sMoreInfo = s.str();
		stcErr.iErrorCode = BAD_DATA;
		throw(stcErr);
	}

	//If lognormal or normal, populate the parameters array
	if (lognormal_pdf == m_iStochasticGrowthMethod
			|| normal_pdf == m_iStochasticGrowthMethod) {

		mp_fRandParameter = new float[m_iNumTotalSpecies];
		FillSpeciesSpecificValue(p_oElement, "gr_standardDeviation", "gr_sdVal",
				p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true);

		//Transfer to the appropriate array buckets
		for (i = 0; i < m_iNumBehaviorSpecies; i++)
			mp_fRandParameter[p_fTempValues[i].code] = p_fTempValues[i].val;

	}

	delete[] p_fTempValues;

	m_sQuery = FormatSpeciesTypeQueryString();
	ReadParameterFile(p_oElement, p_oPop, this, true);
	GetTreeMemberCodes();

}

//////////////////////////////////////////////////////////////////////////////
// SetNameData()
//////////////////////////////////////////////////////////////////////////////
void clJanMinPrecipGrowth::SetNameData(std::string sNameString) {
	try {

		//Check the string passed and set the flags accordingly
		if (sNameString.compare("JanMinPrecipGrowth") == 0) {
			m_iGrowthMethod = diameter_auto;
		} else if (sNameString.compare("JanMinPrecipGrowth diam only") == 0) {
			m_iGrowthMethod = diameter_only;
		} else {
			modelErr stcErr;
			stcErr.iErrorCode = BAD_DATA;
			std::stringstream s;
			s << "Unrecognized behavior name \"" << sNameString << "\".";
			stcErr.sMoreInfo = s.str();
			stcErr.sFunction = "clJanMinPrecipGrowth::SetNameData";
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
		stcErr.sFunction = "clJanMinPrecipGrowth::SetNameData";
		throw(stcErr);
	}
}

////////////////////////////////////////////////////////////////////////////
// CalcGrowthValue
////////////////////////////////////////////////////////////////////////////
float clJanMinPrecipGrowth::CalcDiameterGrowthValue(clTree *p_oTree,
		clTreePopulation *p_oPop, float fHeightGrowth) {
	float fAmountDiamIncrease; //amount diameter increase

	//Get the tree's growth - it's already calculated
	p_oTree->GetValue(
			mp_iGrowthCodes[p_oTree->GetSpecies()][p_oTree->GetType()],
			&fAmountDiamIncrease);

	return fAmountDiamIncrease;
}

////////////////////////////////////////////////////////////////////////////
// PreGrowthCalcs
////////////////////////////////////////////////////////////////////////////
void clJanMinPrecipGrowth::PreGrowthCalcs(clTreePopulation *p_oPop) {
//  float *p_fTempEffect,
//        *p_fWaterEffect;
	try {
#ifdef NCI_WRITER
    using namespace std;
    int iTS = mp_oSimManager->GetCurrentTimestep();
    char cFilename[100];
    sprintf(cFilename, "%s%d%s", "GrowthNCI", iTS, ".txt");
    fstream out( cFilename, ios::trunc | ios::out );
    out << "Timestep\tSpecies\tDBH\tNCI\tGrowth\n";
#endif

		clTreeSearch *p_oNCITrees; //trees that this growth behavior applies to
		clAllometry *p_oAllom = p_oPop->GetAllometryObject();
		clPlot *p_oPlot = mp_oSimManager->GetPlotObject();
		clTree *p_oTree; //a single tree we're working with
		clNCITermBase::ncivals nci;
		float
				fDiam, //tree's diameter
				fX, fY, fNumberYearsPerTimestep =
						mp_oSimManager->GetNumberOfYearsPerTimestep(),
				fAmountDiamIncrease, //amount diameter increase
				fTempDiamIncrease, //amount diameter increase - intermediate
				fPlotTemp = p_oPlot->GetJanMinTemp(), //Temperature parameter
				fPlotWater = p_oPlot->GetMeanAnnualPrecip(); //water availability parameter
		int iIsDead;
		short int iSpecies, iType, //type and species of a tree
				i, //loop counter
				iDeadCode; //tree's dead code

		//Pre-calcs for other effects
		mp_oNCITerm->PreCalcs(p_oPop);

		p_oNCITrees = p_oPop->Find(m_sQuery);

		//************************************
		// Loop through and to calculate growth for each tree
		//************************************
		p_oTree = p_oNCITrees->NextTree();
		while (p_oTree) {
			iSpecies = p_oTree->GetSpecies();
			iType = p_oTree->GetType();

			if (-1 != mp_iGrowthCodes[iSpecies][iType]) {

				//Make sure tree's not dead
				iDeadCode = p_oPop->GetIntDataCode("dead",
						p_oTree->GetSpecies(), p_oTree->GetType());
				if (-1 != iDeadCode) {
					p_oTree->GetValue(iDeadCode, &iIsDead);
				} else
					iIsDead = notdead;

				if (notdead == iIsDead) {

					if (iType == clTreePopulation::seedling) {
						p_oTree->GetValue(
								p_oPop->GetDiam10Code(iSpecies, iType), &fDiam);
					} else {
						p_oTree->GetValue(p_oPop->GetDbhCode(iSpecies, iType),
								&fDiam);
					}

					//First calculate the pieces that have no DBH component and thus will
					//not change in our loop

					//Get NCI
					p_oTree->GetValue(p_oPop->GetXCode(iSpecies, iType), &fX);
					p_oTree->GetValue(p_oPop->GetYCode(iSpecies, iType), &fY);
					nci = mp_oNCITerm->CalculateNCITerm(p_oTree, p_oPop,
							p_oPlot, fX, fY, iSpecies);

					//To correctly compound growth over the number of years per timestep,
					//we have to loop over the number of years, re-calculating the parts
					//with DBH and incrementing the DBH each time
					fAmountDiamIncrease = 0;
					for (i = 0; i < fNumberYearsPerTimestep; i++) {

						//Calculate actual growth in cm/yr
						fTempDiamIncrease = mp_fIntercept[iSpecies]
								+ mp_fDBHSlope[iSpecies] * fDiam
								+ mp_fNCISlope[iSpecies] * nci.fNCI1
								+ mp_fJanMinSlope[iSpecies] * fPlotTemp
								+ mp_fPrecipSlope[iSpecies] * fPlotWater;
						fTempDiamIncrease =
								(fTempDiamIncrease > 0) ? fTempDiamIncrease : 0;

						//Add it to the running total of diameter increase
						fAmountDiamIncrease += fTempDiamIncrease;

						//Increase the DBH for the next loop.  If this is a sapling,
						//convert to a dbh value from what would be a diam10 increase
						// why?
						//if ( clTreePopulation::sapling == iType )
						//{
						//  fDiam += p_oAllom->ConvertDiam10ToDbh( fTempDiamIncrease, iSpecies );
						//}
						//else
						//{
						fDiam += fTempDiamIncrease;
						//}
					}

					//Adjust stochastically (if deterministic, we'll get the same
					//number back
					fAmountDiamIncrease = (*this.*Adjust)(fAmountDiamIncrease,
							iSpecies);

					//Assign the growth back to "Growth" and hold it
					p_oTree->SetValue(mp_iGrowthCodes[iSpecies][iType],
							fAmountDiamIncrease);

#ifdef NCI_WRITER
            p_oTree->GetValue( p_oPop->GetDbhCode( iSpecies, iType ), & fDbh );
            out << iTS << "\t" << p_oTree->GetSpecies() << "\t" << fDbh
                << "\t" << fNCI << "\t" << fAmountDiamIncrease << "\n";
#endif

				} //end of if (bIsNCITree)
			}

			p_oTree = p_oNCITrees->NextTree();
		}

//    delete[] p_fTempEffect;
//    delete[] p_fWaterEffect;

#ifdef NCI_WRITER
    out.close();
#endif

	} catch (modelErr &err) {
//    delete[] p_fTempEffect;
//    delete[] p_fWaterEffect;
		throw(err);
	} catch (...) {
//    delete[] p_fTempEffect;
		//   delete[] p_fWaterEffect;
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clJanMinPrecipGrowth::PreCalcGrowth";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////
// DeterministicAdjust()
/////////////////////////////////////////////////////////////////////////////
float clJanMinPrecipGrowth::DeterministicAdjust(float fNumber, int iSpecies) {
	return fNumber;
}

//////////////////////////////////////////////////////////////////////////////
// NormalAdjust()
/////////////////////////////////////////////////////////////////////////////
float clJanMinPrecipGrowth::NormalAdjust(float fNumber, int iSpecies) {
	float fReturn = fNumber
			+ clModelMath::NormalRandomDraw(mp_fRandParameter[iSpecies]);
	if (fReturn < 0)
		fReturn = 0;
	return fReturn;
}

//////////////////////////////////////////////////////////////////////////////
// LognormalAdjust()
/////////////////////////////////////////////////////////////////////////////
float clJanMinPrecipGrowth::LognormalAdjust(float fNumber, int iSpecies) {
	float fReturn = clModelMath::LognormalRandomDraw(fNumber,
			mp_fRandParameter[iSpecies]);
	if (fReturn < 0)
		fReturn = 0;
	return fReturn;
}


/*
 * ParamChange.cpp
 *
 *  Created on: Dec 6, 2017
 *      Author: ucmuser
 */

#include "ParamChange.h"
#include "SimManager.h"
#include "TreePopulation.h"
#include "ParsingFunctions.h"
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////
clParamChange::clParamChange(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager) {
	try {
		m_sNameString = "ParamChange";
		m_sXMLRoot = "ParamChange";

		//Allowed file types
		m_iNumAllowedTypes = 2;
		mp_iAllowedFileTypes = new int[m_iNumAllowedTypes];
		mp_iAllowedFileTypes[0] = parfile;
		mp_iAllowedFileTypes[1] = detailed_output;

		//Versions
		m_fVersionNumber = 1;
		m_fMinimumVersionNumber = 1;

		m_fH = 1;
	}

	catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clParamChange::clParamChange";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////////////
// GetData
//////////////////////////////////////////////////////////////////////////////////////
void clParamChange::GetData(DOMDocument *p_oDoc) {
	try {
		DOMElement *p_oElement = GetParentParametersElement(p_oDoc);

		FillSingleValue(p_oElement, "sc_ParamChangeH", &m_fH, true);
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clParamChange::GetData";
		throw(stcErr);
	}
}

//////////////////////////////////////////////////////////////////////////////////////
// Action
//////////////////////////////////////////////////////////////////////////////////////
void clParamChange::Action() {
	try {
		clTreeSearch *p_oAllTrees; //search object for getting all trees
		clTree *p_oTree; //for working with a single tree
		clTreePopulation *p_oPop =
				(clTreePopulation*) mp_oSimManager->GetPopulationObject(
						"treepopulation");
		int iDead, i;
//	iNumSpecies = p_oPop->GetNumberOfSpecies();
		int iNbSp; //number of trees from each species
		float fParamG1;
		float fMeanP1; // mean of parameter 1
		short unsigned int iTp; // type of a given tree
		short int iSp, //species of a given tree
				iDeadCode; //dead code for a tree
		std::stringstream sQueryTemp;
		std::string sQuery; // the Query for the find function

		// for each species
		for (i = 0; i < m_iNumBehaviorSpecies; i++) {

			iNbSp = 0;

			// Ask the population to find all adult trees from the species
			sQueryTemp << "species=" << mp_iWhatSpecies[i] << "::type="
					<< clTreePopulation::adult;
			sQuery = sQueryTemp.str();
			p_oAllTrees = p_oPop->Find(sQuery);

			//Go through the trees one at a time
			p_oTree = p_oAllTrees->NextTree();
			while (p_oTree) {
				//Cache tree species
				iSp = p_oTree->GetSpecies();

				//Make sure this tree is not dead from a previous disturbance
				iDeadCode = p_oPop->GetIntDataCode("dead", iSp,
						p_oTree->GetType());
				if (-1 != iDeadCode) {
					p_oTree->GetValue(iDeadCode, &iDead);
					if (iDead > notdead)
						goto nextTree;
				}

				//Make sure the tree is an adult
				iTp = p_oTree->GetType();
				if (3 == iTp) {
					//Get the tree's diameter
					p_oTree->GetValue(p_oPop->GetParamG1Code(iSp, iTp),
							&fParamG1);
					iNbSp++;
					fMeanP1 += fParamG1;
				}

				nextTree: p_oTree = p_oAllTrees->NextTree();
			} //end of while (p_oTree)

			// compute the mean by dividing by the number of adults
			fMeanP1 = fMeanP1 / iNbSp;

			// change the parameter in the population
			p_oPop->SetMeanParamG1(iSp, fMeanP1);

			//reset the stringstream
			sQueryTemp.str("");
			sQueryTemp.clear();
		}
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} //non-fatal error
	catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clParamChange::Action";
		throw(stcErr);
	}
}


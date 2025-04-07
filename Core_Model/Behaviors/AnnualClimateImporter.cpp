/*
 * AnnualClimateImporter.cpp
 *
 *  Created on: Jul 26, 2017
 *      Author: ucmuser
 */

#include "AnnualClimateImporter.h"
#include "SimManager.h"
#include "TreePopulation.h"
#include "ParsingFunctions.h"
#include "Plot.h"
#include <sstream>

// Constructor

clAnnualClimateImporter::clAnnualClimateImporter(clSimManager *p_oSimManager) :
		clWorkerBase(p_oSimManager), clBehaviorBase(p_oSimManager) {
	try {
		m_sNameString = "AnnualClimateImporter";
		m_sXMLRoot = "AnnualClimateImporter";

		m_iNumAllowedTypes = 2;
		mp_iAllowedFileTypes = new int[m_iNumAllowedTypes];
		mp_iAllowedFileTypes[0] = parfile;
		mp_iAllowedFileTypes[1] = detailed_output;

		m_fVersionNumber = 1;
		m_fMinimumVersionNumber = 1;

		mp_fPpt = NULL;
		mp_fJanMin = NULL;
		mp_fJulMax = NULL;
	} catch (modelErr &err) {
		throw(err);
	} catch (modelMsg &msg) {
		throw(msg);
	} catch (...) {
		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clAnnualClimateImporter::clAnnualClimateImporter";
		throw(stcErr);
	}
}

// Destructor

clAnnualClimateImporter::~clAnnualClimateImporter() {
	delete[] mp_fJanMin;
	delete[] mp_fPpt;
	delete[] mp_fJulMax;
}

// GetData

void clAnnualClimateImporter::GetData(DOMDocument *p_oDoc) {
	DOMElement *p_oElement = GetParentParametersElement(p_oDoc);
	int iNumTimesteps = mp_oSimManager->GetNumberOfTimesteps(), i;

	double *p_fJanMin = new double[iNumTimesteps], *p_fPpt =
			new double[iNumTimesteps], *p_fJulMax = new double[iNumTimesteps];

	try {
		// Get parameter file data
		ReadParameterFileData(p_oElement, p_fJanMin, p_fPpt, p_fJulMax);

		if (p_fPpt < 0) {
			modelErr stcErr;
			stcErr.iErrorCode = BAD_DATA;
			stcErr.sFunction = "clAnnualClimateImporter::ReadParameterFileData";
			stcErr.sMoreInfo = "Precipiation cannot be negative.";
			throw(stcErr);
		}

		mp_fJanMin = new double[iNumTimesteps];
		mp_fPpt = new double[iNumTimesteps];
		mp_fJulMax = new double[iNumTimesteps];

		for (i = 0; i < iNumTimesteps; i++) {
			mp_fJanMin[i] = p_fJanMin[i];
			mp_fPpt[i] = p_fPpt[i];
			mp_fJulMax[i] = p_fJulMax[i];
		}

		// Set the initial conditions values according to the first timestep
		clPlot *p_oPlot = mp_oSimManager->GetPlotObject();
		p_oPlot->SetJanMinTemp(mp_fJanMin[0]);
		p_oPlot->SetMeanAnnualPrecip(mp_fPpt[0]);
		p_oPlot->SetJulMaxTemp(mp_fJulMax[0]);

		delete[] p_fJanMin;
		delete[] p_fPpt;
		delete[] p_fJulMax;
	} catch (modelErr &err) {
		delete[] p_fJanMin;
		delete[] p_fPpt;
		delete[] p_fJulMax;
		throw(err);

	} catch (...) {
		delete[] p_fJanMin;
		delete[] p_fPpt;
		delete[] p_fJulMax;

		modelErr stcErr;
		stcErr.iErrorCode = UNKNOWN;
		stcErr.sFunction = "clDisturbance::CutTrees";
		throw(stcErr);
	}
}

// Action
void clAnnualClimateImporter::Action() {
	clPlot *p_oPlot = mp_oSimManager->GetPlotObject();
	int iTs = mp_oSimManager->GetCurrentTimestep() - 1;

	p_oPlot->SetJanMinTemp(mp_fJanMin[iTs]);
	p_oPlot->SetMeanAnnualPrecip(mp_fPpt[iTs]);
	p_oPlot->SetJulMaxTemp(mp_fJulMax[iTs]);
}

// ReadParameterFileData
void clAnnualClimateImporter::ReadParameterFileData(DOMElement *p_oElement,
		double *mp_fJanMin, double *mp_fPpt, double *mp_fJulMax) {

//	double fVal;
//	int i, j, iNumTimesteps = mp_oSimManager->GetNumberOfTimesteps();

// January minimu temperature for each year
	ReadAnnualData(p_oElement, "sc_aciJanMin", "sc_aciJanMinVal", mp_fJanMin);

	// Precipitation for each year
	ReadAnnualData(p_oElement, "sc_aciPpt", "sc_aciPptVal", mp_fPpt);

	// July max temperature for each year
	ReadAnnualData(p_oElement, "sc_aciJulMax", "sc_aciJulMaxVal", mp_fJulMax);

	// check our data
	if (mp_fPpt < 0) {
		modelErr stcErr;
		stcErr.iErrorCode = BAD_DATA;
		stcErr.sFunction = "clAnnualClimateImporter::ReadParameterFileData";
		stcErr.sMoreInfo = "Annual precipitation is negative.";
		throw(stcErr);
	}
}

///////////////////////////////////////////////////////////////////////////////
// ReadAnnualData
///////////////////////////////////////////////////////////////////////////////
void clAnnualClimateImporter::ReadAnnualData(DOMElement *p_oParent,
		std::string sParentTag, std::string sSubTag, double *p_fVal) {

	DOMNodeList *p_oNodeList;
	DOMElement *p_oElement, *p_oChildElement;
	DOMNode *p_oDocNode;
	XMLCh *sVal;
	char *cName;
	double fVal;
	int iNumChildren, iNumNodes, i, iCode;

	//Find the parent element
	sVal = XMLString::transcode(sParentTag.c_str());
	p_oNodeList = p_oParent->getElementsByTagName(sVal);
	XMLString::release(&sVal);
	iNumNodes = p_oNodeList->getLength();
	if (0 == iNumNodes) {
		modelErr stcErr;
		stcErr.sFunction = "AnnualClimateImporter::ReadAnnualData";
		stcErr.iErrorCode = DATA_MISSING;
		stcErr.sMoreInfo = sParentTag;
		throw(stcErr);
	}
	p_oDocNode = p_oNodeList->item(0); //get first element returned by the list

	//Now get the list of all subelements with the actual data
	p_oElement = (DOMElement*) p_oDocNode; //cast to DOM_Element
	sVal = XMLString::transcode(sSubTag.c_str());
	p_oNodeList = p_oElement->getElementsByTagName(sVal);
	XMLString::release(&sVal);
	iNumChildren = p_oNodeList->getLength();
	if (iNumChildren != mp_oSimManager->GetNumberOfTimesteps()) {
		modelErr stcErr;
		stcErr.iErrorCode = BAD_DATA;
		stcErr.sFunction = "clAnnualClimateImporter::ReadParameterFileData";
		stcErr.sMoreInfo = "Not enough timesteps of data.";
		throw(stcErr);
	}

	for (i = 0; i < iNumChildren; i++) {
		p_oDocNode = p_oNodeList->item(i);
		p_oChildElement = (DOMElement*) p_oDocNode;
		sVal = XMLString::transcode("ts");
		cName = XMLString::transcode(
				p_oChildElement->getAttributeNode(sVal)->getNodeValue());
		XMLString::release(&sVal);
		iCode = atoi(cName);
		if (-1 == iCode) { //throw an error
			modelErr stcErr;
			stcErr.iErrorCode = BAD_DATA;
			stcErr.sFunction = "clAnnualClimateImporter::ReadParameterFileData";
			stcErr.sMoreInfo =
					"Unrecognized timestep in annual climate importer.";
			throw(stcErr);
		}
		cName = XMLString::transcode(
				p_oChildElement->getChildNodes()->item(0)->getNodeValue());
		fVal = strtod(cName, NULL);
		if (fVal != 0 || cName[0] == '0') {
			p_fVal[(iCode - 1)] = fVal;
		}
	}
}

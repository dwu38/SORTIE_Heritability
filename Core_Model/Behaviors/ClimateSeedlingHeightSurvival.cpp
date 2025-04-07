/*
 * ClimateSeedlingHeightSurvival.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: ucmuser
 *
 *      Added seedling size 0 - DW
 *
 */
//---------------------------------------------------------------------------
#include "ClimateSeedlingHeightSurvival.h"
#include "TreePopulation.h"
#include "SimManager.h"
#include "ParsingFunctions.h"
#include "Plot.h"
#include "Allometry.h"
#include "Grid.h"
#include <stdio.h>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////////////
clClimateSeedlingHeightSurvival::clClimateSeedlingHeightSurvival( clSimManager * p_oSimManager ) :
  clWorkerBase( p_oSimManager ), clBehaviorBase( p_oSimManager ),
  clMortalityBase( p_oSimManager )
{
  //Set namestring
  m_sNameString = "Seedlingsmortshell";
  m_sXMLRoot = "ClimateSeedlingHeightSurvival";

  //Null out our pointers
  mp_fInter = NULL;
	mp_fSize0 = NULL; //added seedling size 0
  mp_fSize1 = NULL;
  mp_fSize2 = NULL;
  mp_fSize3 = NULL;
  mp_fJulyMaxTemp = NULL;
  mp_fPrecip = NULL;
  mp_fBA = NULL;
  mp_fBAt = NULL;
  mp_iGridSurvivalCodes = NULL;

  mp_oGrid = NULL;
  m_fRadius = 10;
  m_iBATCode = -1;
  m_fMinSaplingHeight = 0;

  //Version 1
  m_fVersionNumber = 1.0;
  m_fMinimumVersionNumber = 1.0;
}


////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////
clClimateSeedlingHeightSurvival::~clClimateSeedlingHeightSurvival()
{
  delete[] mp_fInter;
	delete[] mp_fSize0; //added seedling size 0
  delete[] mp_fSize1;
  delete[] mp_fSize2;
  delete[] mp_fSize3;
  delete[] mp_fJulyMaxTemp;
  delete[] mp_fPrecip;
  delete[] mp_fBA;
  delete[] mp_fBAt;
  delete[] mp_iGridSurvivalCodes;
}


////////////////////////////////////////////////////////////////////////////
// ReadParameterFile()
////////////////////////////////////////////////////////////////////////////
void clClimateSeedlingHeightSurvival::ReadParameterFile( xercesc::DOMDocument * p_oDoc, clTreePopulation *p_oPop )
{
  try
  {
    DOMElement * p_oElement = GetParentParametersElement(p_oDoc);
    floatVal * p_fTempValues; //for getting species-specific values
    short int i, //loop counter
      iNumTotalSpecies = p_oPop->GetNumberOfSpecies();

    m_fMinSaplingHeight = 50;
    //Get the minimum sapling height
    for ( i = 0; i < iNumTotalSpecies; i++ )
      if ( p_oPop->GetMaxSeedlingHeight( i ) < m_fMinSaplingHeight )
        m_fMinSaplingHeight = p_oPop->GetMaxSeedlingHeight( i );

    //The rest are sized number of species to which this behavior applies
    mp_fInter = new float[iNumTotalSpecies];
	mp_fSize0 = new float[iNumTotalSpecies]; //added seedling size 0
    mp_fSize1 = new float[iNumTotalSpecies];
    mp_fSize2 = new float[iNumTotalSpecies];
    mp_fSize3 = new float[iNumTotalSpecies];
    mp_fJulyMaxTemp = new float[iNumTotalSpecies];
    mp_fPrecip = new float[iNumTotalSpecies];
    mp_fBA = new float[iNumTotalSpecies];
    mp_fBAt = new float[iNumTotalSpecies];

    //Set up our floatVal array that will extract values only for the species
    //assigned to this behavior
    p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
      p_fTempValues[i].code = mp_iWhatSpecies[i];

    //Fill the variables

    // intercept parameter
    FillSpeciesSpecificValue( p_oElement, "mo_SeedlingsInter", "mo_SeedlingsInterVal",
        p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
    //Transfer to the appropriate array buckets
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
    	mp_fInter[p_fTempValues[i].code] = p_fTempValues[i].val;

    // size 0 parameter
    FillSpeciesSpecificValue( p_oElement, "mo_SeedlingsS0", "mo_SeedlingsS0Val",
        p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true ); //added seedling size 0
    //Transfer to the appropriate array buckets
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
    	mp_fSize0[p_fTempValues[i].code] = p_fTempValues[i].val;

    // size 1 parameter
    FillSpeciesSpecificValue( p_oElement, "mo_SeedlingsS1", "mo_SeedlingsS1Val",
        p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
    //Transfer to the appropriate array buckets
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
    	mp_fSize1[p_fTempValues[i].code] = p_fTempValues[i].val;

    // size 2 parameter
    FillSpeciesSpecificValue( p_oElement, "mo_SeedlingsS2", "mo_SeedlingsS2Val",
        p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
    //Transfer to the appropriate array buckets
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
    	mp_fSize2[p_fTempValues[i].code] = p_fTempValues[i].val;

    // size 3 parameter
    FillSpeciesSpecificValue( p_oElement, "mo_SeedlingsS3", "mo_SeedlingsS3Val",
        p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
    //Transfer to the appropriate array buckets
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
    	mp_fSize3[p_fTempValues[i].code] = p_fTempValues[i].val;

    // average temperature
    FillSpeciesSpecificValue( p_oElement, "mo_SeedlingsJulyMaxTemp", "mo_SeedlingsJMTVal",
        p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
    //Transfer to the appropriate array buckets
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
      mp_fJulyMaxTemp[p_fTempValues[i].code] = p_fTempValues[i].val;

    // water deficit
    FillSpeciesSpecificValue( p_oElement, "mo_SeedlingsPrecip", "mo_SeedlingsPVal",
        p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
    //Transfer to the appropriate array buckets
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
      mp_fPrecip[p_fTempValues[i].code] = p_fTempValues[i].val;

    // basal area
    FillSpeciesSpecificValue( p_oElement, "mo_SeedlingsBA", "mo_SeedlingsBAVal",
        p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
    //Transfer to the appropriate array buckets
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
      mp_fBA[p_fTempValues[i].code] = p_fTempValues[i].val;

    // basal area threshold
    FillSpeciesSpecificValue( p_oElement, "mo_SeedlingsBAt", "mo_SeedlingsBAtVal",
        p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
    //Transfer to the appropriate array buckets
    for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
      mp_fBAt[p_fTempValues[i].code] = p_fTempValues[i].val;

    //Max crowding radius
    FillSingleValue( p_oElement, "mo_SeedlingsRadius", &m_fRadius, true );

    delete[] p_fTempValues;
  }
  catch ( modelErr & err )
  {
    throw( err );
  }
  catch ( modelMsg & msg )
  {
    throw( msg );
  } //non-fatal error
  catch ( ... )
  {
    modelErr stcErr;
    stcErr.iErrorCode = UNKNOWN;
    stcErr.sFunction = "clClimateSeedlingHeightSurvival::ReadParameterFile" ;
    throw( stcErr );
  }
}

////////////////////////////////////////////////////////////////////////////
// ValidateData
////////////////////////////////////////////////////////////////////////////
void clClimateSeedlingHeightSurvival::ValidateData()
{
  try
  {

    //Make sure that the neighbor search radius is not less than zero
    if ( 0 > m_fRadius )
    {
      modelErr stcErr;
      stcErr.sMoreInfo = "Neighbor search radius cannot be less than 0.";
      throw( stcErr );
    }
  }
  catch ( modelErr & err )
  {
    err.sFunction = "clClimateSeedlingHeightSurvival::ValidateData";
    err.iErrorCode = BAD_DATA;
    throw( err );
  }
  catch ( modelMsg & msg )
  {
    throw( msg );
  } //non-fatal error
  catch ( ... )
  {
    modelErr stcErr;
    stcErr.iErrorCode = UNKNOWN;
    stcErr.sFunction = "clClimateSeedlingHeightSurvival::ValidateData" ;
    throw( stcErr );
  }
}

////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
////////////////////////////////////////////////////////////////////////////
void clClimateSeedlingHeightSurvival::DoShellSetup( xercesc::DOMDocument * p_oDoc )
{
  clTreePopulation * p_oPop = ( clTreePopulation * ) mp_oSimManager->GetPopulationObject( "treepopulation" );
  ReadParameterFile( p_oDoc, p_oPop );
  ValidateData();
  //SetupGrid(p_oPop);
  mp_fInter = new float[p_oPop->GetNumberOfSpecies()];
	mp_fSize0 = new float[p_oPop->GetNumberOfSpecies()]; //added seedling size 0
  mp_fSize1 = new float[p_oPop->GetNumberOfSpecies()];
  mp_fSize2 = new float[p_oPop->GetNumberOfSpecies()];
  mp_fSize3 = new float[p_oPop->GetNumberOfSpecies()];
  mp_fJulyMaxTemp = new float[p_oPop->GetNumberOfSpecies()];
  mp_fPrecip = new float[p_oPop->GetNumberOfSpecies()];
  mp_fBA = new float[p_oPop->GetNumberOfSpecies()];
  mp_fBAt = new float[p_oPop->GetNumberOfSpecies()];
}

////////////////////////////////////////////////////////////////////////////
// DoMort
////////////////////////////////////////////////////////////////////////////
deadCode clClimateSeedlingHeightSurvival::DoMort( clTree * p_oTree, const float & fDbh, const short int & iSpecies )
{
  clTreePopulation * p_oPop = ( clTreePopulation * ) mp_oSimManager->GetPopulationObject( "treepopulation" );
  clPlot *p_oPlot = mp_oSimManager->GetPlotObject();
  float fX, fY,
        fNumberYearsPerTimestep = mp_oSimManager->GetNumberOfYearsPerTimestep(),
        fSurvivalProb, //Tree's annual survival probability
        fBAT, //Neighborhood basal area of adult trees
  	  	ftheta, // linear function that has to be pass to the logit link
		fJulyMaxTemp = p_oPlot->GetJulMaxTemp(), //Current plot temperature!!
        fPrecip = p_oPlot->GetMeanAnnualPrecip(),//Current plot precipitation
  	  	//fHeight, // height of the current tree
  	  	fheightEffect; // height effect
  short int iType = p_oTree->GetType();
 //           iX, iY;
  int iSizeClass; // size class of the current tree
  bool bIsDead;

//  //Has survival been calculated for this tree's cell?
  p_oTree->GetValue( p_oPop->GetXCode( iSpecies, iType ), & fX );
  p_oTree->GetValue( p_oPop->GetYCode( iSpecies, iType ), & fY );
//  mp_oGrid->GetValueAtPoint( fX, fY, mp_iGridSurvivalCodes[iSpecies], &fSurvivalProb);
//  if ( -1 == fSurvivalProb ) {

    //Get the basal area of neighbors, if it has not already been gotten
//    mp_oGrid->GetCellOfPoint(fX, fY, &iX, &iY);
//    mp_oGrid->GetPointOfCell( iX, iY, & fX, & fY );
//    mp_oGrid->GetValueOfCell(iX, iY, m_iBATCode, &fBAT);
//    if (fBAT < 0) {
      fBAT = GetBAT(fX, fY, p_oPop);
//      mp_oGrid->SetValueOfCell(iX, iY, m_iBATCode, fBAT);
//    }

    // Get tree's size class and choose the height effect
    p_oTree->GetValue( p_oPop->GetSizeClassCode( iSpecies, iType ), & iSizeClass );
	if (iSizeClass <= 0) {
		fheightEffect = mp_fSize0[iSpecies]; //added seedling size 0
	}
    if (iSizeClass == 1){
    	fheightEffect = mp_fSize1[iSpecies];
    }
    else if (iSizeClass == 2){
    	fheightEffect = mp_fSize2[iSpecies];
    }
    else if (iSizeClass == 3){
    	fheightEffect = mp_fSize3[iSpecies];
    }

    //Calculate survival probability for this tree
    ftheta = mp_fInter[iSpecies] + fheightEffect + mp_fJulyMaxTemp[iSpecies] * fJulyMaxTemp/10 + mp_fPrecip[iSpecies] * fPrecip/1000 + mp_fBA[iSpecies] * fBAT;
    fSurvivalProb = 1/(1+exp(-ftheta));
    // if the basal area is above threshold, the seedling will automatically die
    if (fBAT > mp_fBAt[iSpecies])
    	fSurvivalProb = 0;
    fSurvivalProb = pow( fSurvivalProb, fNumberYearsPerTimestep );
//    mp_oGrid->SetValueOfCell(iX, iY, mp_iGridSurvivalCodes[iSpecies], fSurvivalProb);
//  }

  bIsDead = clModelMath::GetRand() >= fSurvivalProb;

  if (bIsDead) return natural;
  else return notdead;
}



////////////////////////////////////////////////////////////////////////////
// GetBAT
////////////////////////////////////////////////////////////////////////////
float clClimateSeedlingHeightSurvival::GetBAT(float &fX, float &fY, clTreePopulation *p_oPop)
{
  clTreeSearch * p_oAllNeighbors; //neighborhood trees within crowding radius
  clTree * p_oNeighbor; //competing neighbor
  char cQuery[75]; //format search strings into this
  float fDbh, //neighbor's dbh
        fBAT = 0;
  int iIsDead; //whether a neighbor is dead
  short int iSpecies, iType, //species and type for neighbor
            iDeadCode; //neighbor's dead code

  //Get all trees taller than seedlings within the max crowding radius -
  //seedlings don't compete
  sprintf( cQuery, "%s%f%s%f%s%f%s%f", "distance=", m_fRadius, "FROM x=", fX,
       "y=", fY, "::height=", m_fMinSaplingHeight );
  p_oAllNeighbors = p_oPop->Find( cQuery );

  //Loop through and find all the adults
  p_oNeighbor = p_oAllNeighbors->NextTree();

  while ( p_oNeighbor ) {

    iType = p_oNeighbor->GetType();
    if ( clTreePopulation::adult == iType) {
      iSpecies = p_oNeighbor->GetSpecies();
      //Get the neighbor's dbh
      p_oNeighbor->GetValue( p_oPop->GetDbhCode( iSpecies, iType ), & fDbh );

      //Make sure the neighbor's not dead due to a disturbance event
      iDeadCode = p_oPop->GetIntDataCode( "dead", iSpecies, iType );
      if ( -1 != iDeadCode )  p_oNeighbor->GetValue( iDeadCode, & iIsDead );
      else iIsDead = notdead;

      if ( notdead == iIsDead || natural == iIsDead)
        fBAT += clModelMath::CalculateBasalArea(fDbh);

    }

    p_oNeighbor = p_oAllNeighbors->NextTree();
  }
  return fBAT;
}

////////////////////////////////////////////////////////////////////////////
// PreMortCalcs
////////////////////////////////////////////////////////////////////////////
//void clClimateSeedlingHeightSurvival::PreMortCalcs( clTreePopulation * p_oPop )
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
//void clClimateSeedlingHeightSurvival::SetupGrid(clTreePopulation *p_oPop)
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
//        stcErr.sFunction = "clClimateSeedlingHeightSurvival::SetupGrid";
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
//        stcErr.sFunction = "clClimateSeedlingHeightSurvival::SetupGrid" ;
//        stcErr.sMoreInfo = "Couldn't find the \"BAT\" member of the \"Seedlings Survival\" grid.";
//        stcErr.iErrorCode = BAD_DATA;
//        throw( stcErr );
//      }
//    }
//  }
//
//}
//




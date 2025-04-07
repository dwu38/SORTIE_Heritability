/*
 * ClimateNCIMortality.cpp
 *
 *  Created on: May 14, 2018
 *      Author: ucmuser
 */

//---------------------------------------------------------------------------
#include <ClimateNCIMortality.h>
#include "TreePopulation.h"
#include "SimManager.h"
#include "ParsingFunctions.h"
#include "Plot.h"
#include "GrowthOrg.h"
#include "Allometry.h"

//#include "NCI/CrowdingEffectBase.h"
//#include "NCI/DamageEffectBase.h"
#include "NCI/NCITermBase.h"
//#include "NCI/ShadingEffectBase.h"
//#include "NCI/SizeEffectBase.h"
//#include "NCI/TemperatureEffectBase.h"
//#include "NCI/PrecipitationEffectBase.h"
//#include "NCI/InfectionEffectBase.h"
//#include "NCI/NitrogenEffectBase.h"

#include <stdio.h>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////////////
clClimateNCIMortality::clClimateNCIMortality( clSimManager * p_oSimManager ) :
clWorkerBase( p_oSimManager ), clBehaviorBase( p_oSimManager ), clMortalityBase( p_oSimManager ), clNCIBehaviorBase()
{
  try
  {
    //Set namestring
    m_sNameString = "JulMaxPrecip mortshell";
    m_sXMLRoot = "ClimateNCIMortality";

    //Version 1
    m_fVersionNumber = 1.0;
    m_fMinimumVersionNumber = 1.0;
    m_fMaxSurvivalPeriod = 1;

    //Null out our pointers
    mp_iDeadCodes = NULL;
    //mp_fMaxPotentialValue = NULL;

    m_iNumTotalSpecies = 0;

	mp_fInter = NULL;
	mp_fDiam = NULL;
	mp_fDiam2 = NULL;
	mp_fNCI = NULL;
	mp_fJulMax = NULL;
	mp_fPrecip= NULL;

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
    stcErr.sFunction = "clClimateNCIMortality::clClimateNCIMortality" ;
    throw( stcErr );
  }
}

////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////
clClimateNCIMortality::~clClimateNCIMortality() {

  int i;
  if (mp_iDeadCodes) {
    for ( i = 0; i < m_iNumTotalSpecies; i++ ) {
      delete[] mp_iDeadCodes[i];
    }
    delete[] mp_iDeadCodes;
  }
  //delete[] mp_fMaxPotentialValue;

  delete[] mp_fInter;
  delete[] mp_fDiam;
  delete[] mp_fDiam2;
  delete[] mp_fNCI;
  delete[] mp_fJulMax;
  delete[] mp_fPrecip;
}


////////////////////////////////////////////////////////////////////////////
// GetTreeMemberCodes()
////////////////////////////////////////////////////////////////////////////
void clClimateNCIMortality::GetTreeMemberCodes() {
  clTreePopulation *p_oPop = (clTreePopulation*) mp_oSimManager->GetPopulationObject("treepopulation");
  int iNumTypes = p_oPop->GetNumberOfTypes(),
      i, j;

  //Declare the dead codes array
  mp_iDeadCodes = new short int * [m_iNumTotalSpecies];
  for ( i = 0; i < m_iNumTotalSpecies; i++ ) {
    mp_iDeadCodes[i] = new short int[iNumTypes];
    for (j = 0; j < iNumTypes; j++) {
      mp_iDeadCodes[i][j] = -1;
    }
  }
  for ( i = 0; i < m_iNumSpeciesTypeCombos; i++ ) {
    mp_iDeadCodes[mp_whatSpeciesTypeCombos[i].iSpecies] [mp_whatSpeciesTypeCombos[i].iType] =
        p_oPop->GetIntDataCode( "dead", mp_whatSpeciesTypeCombos[i].iSpecies, mp_whatSpeciesTypeCombos[i].iType );

    if ( -1 == mp_iDeadCodes[mp_whatSpeciesTypeCombos[i].iSpecies]
                             [mp_whatSpeciesTypeCombos[i].iType] )
    {
      modelErr stcErr;
      stcErr.sFunction = "clClimateNCIMortality::GetTreeMemberCodes" ;
      stcErr.sMoreInfo = "Can't find the \"dead\" data member.";
      stcErr.iErrorCode = BAD_DATA;
      throw( stcErr );
    }
  }
}

////////////////////////////////////////////////////////////////////////////
// DoShellSetup()
////////////////////////////////////////////////////////////////////////////
void clClimateNCIMortality::DoShellSetup(xercesc::DOMDocument * p_oDoc) {
  clTreePopulation * p_oPop = ( clTreePopulation * ) mp_oSimManager->GetPopulationObject( "treepopulation" );

  m_iNumTotalSpecies = p_oPop->GetNumberOfSpecies();

  DOMElement * p_oElement = GetParentParametersElement(p_oDoc);
  floatVal * p_fTempValues; //for getting species-specific values
  short int i; //loop counters


  //If any of the types is seedling, error out
  for ( i = 0; i < m_iNumSpeciesTypeCombos; i++ )
    if ( clTreePopulation::sapling != mp_whatSpeciesTypeCombos[i].iType
        && clTreePopulation::adult != mp_whatSpeciesTypeCombos[i].iType )
    {
      modelErr stcErr;
      stcErr.iErrorCode = BAD_DATA;
      stcErr.sFunction = "clClimateNCIMortality::DoShellSetup" ;
      stcErr.sMoreInfo = "This behavior can only be applied to saplings and adults.";
      throw( stcErr );
    }

	// Declare the arrays we'd like read
	mp_fInter = new float[m_iNumBehaviorSpecies];
	mp_fDiam = new float[m_iNumBehaviorSpecies];
	mp_fDiam2 = new float[m_iNumBehaviorSpecies];
	mp_fNCI = new float[m_iNumBehaviorSpecies];
	mp_fJulMax = new float[m_iNumBehaviorSpecies];
	mp_fPrecip= new float[m_iNumBehaviorSpecies];

  //Set up our floatVal array that will extract values only for the species
  //assigned to this behavior
  p_fTempValues = new floatVal[m_iNumBehaviorSpecies];
  for ( i = 0; i < m_iNumBehaviorSpecies; i++ )
    p_fTempValues[i].code = mp_iWhatSpecies[i];

	// parameter a
	FillSpeciesSpecificValue( p_oElement, "mo_JulMaxPrecipInter",
		"mo_JulMaxPrecipInterVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
	for (i = 0; i < m_iNumBehaviorSpecies; i++ )
		mp_fInter[p_fTempValues[i].code] = p_fTempValues[i].val;

	// parameter b
	FillSpeciesSpecificValue( p_oElement, "mo_JulMaxPrecipDiam",
		"mo_JulMaxPrecipDiamVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
	for (i = 0; i < m_iNumBehaviorSpecies; i++ )
		mp_fDiam[p_fTempValues[i].code] = p_fTempValues[i].val;

	// parameter c
	FillSpeciesSpecificValue( p_oElement, "mo_JulMaxPrecipDiam2",
		"mo_JulMaxPrecipDiam2Val", p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
	for (i = 0; i < m_iNumBehaviorSpecies; i++ )
		mp_fDiam2[p_fTempValues[i].code] = p_fTempValues[i].val;

	// parameter d
	FillSpeciesSpecificValue( p_oElement, "mo_JulMaxPrecipNCI",
		"mo_JulMaxPrecipNCIVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
	for (i = 0; i < m_iNumBehaviorSpecies; i++ )
		mp_fNCI[p_fTempValues[i].code] = p_fTempValues[i].val;

	// parameter e
	FillSpeciesSpecificValue( p_oElement, "mo_JulMaxPrecipJulMax",
		"mo_JulMaxPrecipJulMaxVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
	for (i = 0; i < m_iNumBehaviorSpecies; i++ )
		mp_fJulMax[p_fTempValues[i].code] = p_fTempValues[i].val;

	// parameter f
	FillSpeciesSpecificValue( p_oElement, "mo_JulMaxPrecipPrecip",
		"mo_JulMaxPrecipPrecipVal", p_fTempValues, m_iNumBehaviorSpecies, p_oPop, true );
	for (i = 0; i < m_iNumBehaviorSpecies; i++ )
		mp_fPrecip[p_fTempValues[i].code] = p_fTempValues[i].val;


  delete[] p_fTempValues;

  m_sQuery = FormatSpeciesTypeQueryString();
  ReadParameterFile(p_oElement, p_oPop, this, true);
  GetTreeMemberCodes();
}

////////////////////////////////////////////////////////////////////////////
// DoMort
////////////////////////////////////////////////////////////////////////////
deadCode clClimateNCIMortality::DoMort( clTree * p_oTree, const float & fDbh, const short int & iSpecies )
{
  int iDead;
  p_oTree->GetValue(mp_iDeadCodes[iSpecies][p_oTree->GetType()], &iDead);
  return (deadCode)iDead;
}

////////////////////////////////////////////////////////////////////////////
// PreMortCalcs
////////////////////////////////////////////////////////////////////////////
void clClimateNCIMortality::PreMortCalcs( clTreePopulation * p_oPop )
{
  try
  {

#ifdef NCI_WRITER
    using namespace std;
    float fX, fY;
    int iTS = mp_oSimManager->GetCurrentTimestep();
    char cFilename[100];
    sprintf(cFilename, "%s%d%s", "MortJulMaxPrecip", iTS, ".txt");
    fstream out( cFilename, ios::trunc | ios::out );
    out << "Timestep\tSpecies\tDBH\tNCI\tSurv Prob\tDead?\n";
#endif

    clTreeSearch * p_oNCITrees; //trees that this growth behavior applies to
    clPlot * p_oPlot = mp_oSimManager->GetPlotObject();
    clTree * p_oTree; //a single tree we're working with
    clNCITermBase::ncivals nci;
    float fDbh, //tree's dbh
    fX, fY,
    fNumberYearsPerTimestep = mp_oSimManager->GetNumberOfYearsPerTimestep(),
    fSurvivalProb, //amount diameter increase - intermediate
	fPlotJulMax = p_oPlot->GetJulMaxTemp(), //July maximum temperature
	fPlotPrecip = p_oPlot->GetMeanAnnualPrecip(); //Current plot water deficit

    int iDead;
    short int iSpecies, iType; //type and species of a tree
    bool bIsDead;

    //Pre-calcs for other effects
    mp_oNCITerm->PreCalcs(p_oPop);

    p_oNCITrees = p_oPop->Find( m_sQuery );
    p_oTree = p_oNCITrees->NextTree();
    while ( p_oTree )
    {
      iSpecies = p_oTree->GetSpecies();
      iType = p_oTree->GetType();

      if ( mp_bUsesThisMortality[iSpecies][iType] )
      {

        //Make sure tree's not dead
        p_oTree->GetValue( mp_iDeadCodes[iSpecies][iType], & iDead );

        if ( iDead <= natural )
        {

          p_oTree->GetValue( p_oPop->GetDbhCode( iSpecies, iType ), & fDbh );

          //Get NCI
          p_oTree->GetValue( p_oPop->GetXCode( iSpecies, iType ), & fX );
          p_oTree->GetValue( p_oPop->GetYCode( iSpecies, iType ), & fY );
          nci = mp_oNCITerm->CalculateNCITerm(p_oTree, p_oPop, p_oPlot, fX, fY, iSpecies);

          //Determine whether tree survives, using the parameters for a mortality model
          // (that's why it is not the logit expression, which was used to compute a mortality model)
          fSurvivalProb = exp(mp_fInter[iSpecies] + mp_fDiam[iSpecies] * fDbh + mp_fDiam2[iSpecies] * fDbh * fDbh +
        		  mp_fNCI[iSpecies] * nci.fNCI1 +
				  mp_fJulMax[iSpecies] * fPlotJulMax +
				  mp_fPrecip[iSpecies] * fPlotPrecip);
          fSurvivalProb = 1/(1+fSurvivalProb);
          //Get annual survival
          fSurvivalProb = pow( fSurvivalProb, 1/m_fMaxSurvivalPeriod);
          //Get timestep survival
          fSurvivalProb = pow( fSurvivalProb, fNumberYearsPerTimestep );

          bIsDead = clModelMath::GetRand() >= fSurvivalProb;

          //Assign the value back to the dead data member
          if (bIsDead) iDead = natural;
          else iDead = notdead;
          p_oTree->SetValue( mp_iDeadCodes[iSpecies][iType], iDead );

#ifdef NCI_WRITER
          //if (8 == p_oTree->GetSpecies()) {
          out << iTS << "\t"  << p_oTree->GetSpecies() << "\t" << fDbh
              << "\t" << fNCI  << "\t"
              << fSurvivalProb << "\t" << bIsDead << "\n";
          //}
#endif

        } //end of if ( notdead == iDead )
      }

      p_oTree = p_oNCITrees->NextTree();
    }

#ifdef NCI_WRITER
    out.close();
#endif

    }
  catch ( modelMsg & msg )
  {
    throw( msg );
  } //non-fatal error
  catch ( ... )
  {
    modelErr stcErr;
    stcErr.iErrorCode = UNKNOWN;
    stcErr.sFunction = "clClimateNCIMortality::PreMortGrowth" ;
    throw( stcErr );
  }
}




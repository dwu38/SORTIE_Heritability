/*
 * ClimateNCIMortality.h
 *
 *  Created on: May 14, 2018
 *      Author: ucmuser
 */

#ifndef BEHAVIORS_CLIMATENCIMORTALITY_H_
#define BEHAVIORS_CLIMATENCIMORTALITY_H_


#include "MortalityBase.h"
#include "NCI/NCIBehaviorBase.h"

class clGrid;
class clTreePopulation;
class clTree;

/***
 * ClimateNCIMortality
 * This evaluate mortality according to a logistic equation,
 * linking mortality to the actual DBH , competition (NCI), JulMax and Precip
 * The equation used in this behavior is:
* <center>p = logit-1(a + b DBH + c DBH2 + d NCI + e JulMax + f Precip)</center>
* where
* <ul>
* <li>p = annual Probability of mortality</li>
* <li>a,b,c,d,e,f = parameter</li>
* <li>DBH = tree diameter, in cm; diam10 for seedlings, DBH for others</li>
* <li>JulMax = actual evapotranspiration </li>
* <li>Precip = climate water deficit </li>
* <li>NCI is the sum of neighbors DBH/dist</li>
*
* </ul>
*
* This class's namestring is "JulMaxPrecip mortshell".  The parameter file
* call string is "ClimateNCIMortality".
*
 */

class clClimateNCIMortality : virtual public clMortalityBase, clNCIBehaviorBase {

public:

	/*
	 * Constructor
	 * @param p_oSimManager Sim Manager object.
	 */
	clClimateNCIMortality(clSimManager *p_oSimManager);

	/*
	 * Destructor
	 */
	~clClimateNCIMortality();

	  /**
	  * Determines mortality for a tree.
	  *
	  * @param p_oTree Tree being evaluated.
	  * @param fDbh DBH of tree being evaluated.
	  * @param iSpecies Species of tree being evaluated.
	  * @return natural if the tree is to die, notdead if it lives.
	  */
	  deadCode DoMort (clTree *p_oTree, const float &fDbh, const short int &iSpecies);

	  /**
	  * Does setup.
	  * <ol>
	  * <li>ReadParameterFile() is called to read the parameter file's data.</li>
	  * <li>GetTreeMemberCodes() is called to get tree data return codes.</li>
	  * <li>FormatQueryString() is called.</li>
	  * </ol>
	  *
	  * @param p_oDoc DOM tree of parsed input tree.
	  */
	  void DoShellSetup(xercesc::DOMDocument *p_oDoc);

	  /**
	  * Performs calculations before any trees have been killed.  This finds all
	  * trees to which this behavior applies and performs their NCI calculations.
	  * Then, having done all that work, this function goes ahead and assesses the
	  * tree's mortality.  Whether it lives or dies is then stashed in the
	  * "NCI Mort" bool tree data member.
	  *
	  * @param p_oPop Tree population.
	  */
	  void PreMortCalcs( clTreePopulation *p_oPop );

	  protected:

		/* Mortality equation "a" - sized number of behavior species*/
		float *mp_fInter;

		/* Mortality equation "b" - sized number of behavior species*/
		float *mp_fDiam;

		/* Mortality equation "c" - sized number of behavior species*/
		float *mp_fDiam2;

		/* Mortality equation "d" - sized number of behavior species*/
		float *mp_fNCI;

		/* Mortality equation "e" - sized number of behavior species*/
		float *mp_fJulMax;

		/* Mortality equation "f" - sized number of behavior species*/
		float *mp_fPrecip;


	  /**Maximum survival value. Array sized number of species.*/
	  //float *mp_fMaxPotentialValue;

	  /** The length of the time period of the max survival, if needed for
	   * adjustment of survival rates. For instance, if the max survival is for
	   * a 5-year time period, then this value is 5, and the 5th root is taken
	   * of the final survival rate to arrive at the yearly value. 1 indicates
	   * that the max rate is yearly already.*/
	  float m_fMaxSurvivalPeriod;

	  /**Return codes for the "dead" tree int data member variable.  Array size
	   * is number of species by number of tree types (even if not every species
	   * and type is represented).*/
	  short int **mp_iDeadCodes;

	  /**Total number of species - for the destructor */
	  short int m_iNumTotalSpecies;

	  /**For finding trees*/
	  std::string m_sQuery;

	  /**
	  * Gets the return codes for needed tree data members.
	  * @throws modelErr if a light code comes back -1 for any species which uses
	  * the shading effect.
	  */
	  void GetTreeMemberCodes();

};

#endif /* BEHAVIORS_CLIMATENCIMORTALITY_H_ */

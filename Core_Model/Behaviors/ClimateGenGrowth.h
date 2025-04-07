/*
 * ClimateGenGrowth.h
 *
 *  Created on: Oct 3, 2017
 *      Author: ucmuser
 */

#ifndef BEHAVIORS_CLIMATEGENGROWTH_H_
#define BEHAVIORS_CLIMATEGENGROWTH_H_

//---------------------------------------------------------------------------
#include "GrowthBase.h"
#include "NCI/NCIBehaviorBase.h"


class clTree;
class clTreePopulation;

/**
* Change ClimateNCIGrowth to take the tree individual parameters ParamG1 and ParamG2 into account when computing growth
*/

class clClimateGenGrowth : virtual public clGrowthBase, clNCIBehaviorBase {
//note: need the virtual keyword to avoid base class ambiguity.

  public:

  /**
  * Constructor.
  * @param p_oSimManager Sim Manager object.
  */
	clClimateGenGrowth(clSimManager *p_oSimManager);

  /**
  * Destructor.  Frees memory.
  */
  ~clClimateGenGrowth();

  /**
  * Returns the value in the tree's float data member that holds the value
  * that was calculated by PreGrowthCalcs().
  *
  * @param p_oTree Tree to which to apply growth.
  * @param p_oPop Tree population object.
  * @param fHeightGrowth Amount of height growth, in m (ignored).
  * @return Amount of diameter growth, in cm.
  */
  float CalcDiameterGrowthValue(clTree *p_oTree, clTreePopulation *p_oPop, float fHeightGrowth);

  /**
  * Calculates growth for all NCI trees. Climate effects, which do not depend
  * on individual trees, are called once; effects which do not depend on a
  * target tree's size are called once per tree; and effects with a size
  * component are called once per year, looping over the years in a timestep.
  * The values are stashed in the "Growth" tree float data member for later
  * application.

  * This must be called first of any growth stuff, since it uses other trees'
  * DBHs to calculate NCI, and these must be before growth has been applied.
  *
  * Growth per timestep is calculated by looping over the number of years
  * per timestep and incrementing the DBH.
  *
  * @param p_oPop Tree population object.
  */
  void PreGrowthCalcs( clTreePopulation *p_oPop );

  /**
  * Does setup.
  * <ol>
  * <li>ReadParameterFile() is called to read the parameter file's data.</li>
  * <li>GetTreeMemberCodes() is called to get tree data return codes.</li>
  * <li>FormatQueryString() is called.</li>
  * </ol>
  *
  * @param p_oDoc DOM tree of parsed input tree.
  * @throws modelErr if all values of max growth are not greater than 0.
  */
  void DoShellSetup(xercesc::DOMDocument *p_oDoc);

  /**
  * Captures the namestring passed to this behavior. This is overridden from
  * clBehaviorBase so we can capture the namestring passed. Since this class
  * can create multiple kinds of behaviors that function differently, this will
  * capture what kind of behavior this is supposed to be.
  *
  * @param sNameString Behavior's namestring.
  */
  void SetNameData(std::string sNameString);

  protected:

  /**Intercept. Array sized number of species.*/
  float *mp_fIntercept;

  /** DBH slope. Array sized number of species.*/
  float *mp_fDBHSlope;

  /** NCI slope. Array sized number of species.*/
  float *mp_fNCISlope;

  /** Temperature slope. Array sized number of species.*/
  float *mp_fTempSlope;

  /** Precip slope. Array sized number of species.*/
  float *mp_fPrecipSlope;

  /**Standard deviation if normal or lognormal distribution is desired. One for
   * each species.*/
  float *mp_fRandParameter;

  /**Holds return data codes for the "Growth" tree data member. Array size is
   * number of species by number of types.*/
  short int **mp_iGrowthCodes;

  /**Total number of species - for the destructor */
  short int m_iNumTotalSpecies;

  /**For finding trees*/
  std::string m_sQuery;

  /**What stochastic growth distribution applies to this run*/
  pdf m_iStochasticGrowthMethod;

  /**
  * Gets the return codes for needed tree data members.
  * @throws modelErr if a code comes back -1 for any species/type combo to
  * which this behavior is applied.
  */
  void GetTreeMemberCodes();

  /**
  * Performs a deterministic adjustment of growth.
  * @param fNumber Growth to adjust.
  * @param iSpecies Species.
  * @return Final growth.
  */
  float DeterministicAdjust(float fNumber, int iSpecies);

  /**
  * Performs a stochastic adjustment of growth according to a normal
  * distribution.
  * @param fNumber Growth to adjust.
  * @param iSpecies Species.
  * @return Final growth.
  */
  float NormalAdjust(float fNumber, int iSpecies);

  /**
  * Performs a stochastic adjustment of growth according to a lognormal
  * distribution.
  * @param fNumber Growth to adjust.
  * @param iSpecies Species.
  * @return Final growth.
  */
  float LognormalAdjust(float fNumber, int iSpecies);

  /**
  * Sets the Adjust function pointer according to the value of
  * m_iStochasticGrowthMethod.
  */
  void SetFunctionPointer();

  /**Function pointer for the appropriate Adjust function*/
  float (clClimateGenGrowth::*Adjust)(float fNumber, int iSpecies);

};
//---------------------------------------------------------------------------


#endif /* BEHAVIORS_CLIMATEGENGROWTH_H_ */

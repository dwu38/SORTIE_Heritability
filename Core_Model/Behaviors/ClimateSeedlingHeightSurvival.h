/*
 * ClimateSeedlingHeightSurvival.h
 *
 *  Created on: Nov 29, 2016
 *      Author: ucmuser
 *
 *      Added seedling size 0 - DW
 *
 */

#ifndef BEHAVIORS_CLIMATESEEDLINGHEIGHTSURVIVAL_H_
#define BEHAVIORS_CLIMATESEEDLINGHEIGHTSURVIVAL_H_

/*
 * SeedlingsMortality.h
 *
 *  Created on: Nov 29, 2016
 *      Author: ucmuser
 */

//---------------------------------------------------------------------------
#include "MortalityBase.h"
  /**
   * SeedlingsMortality using the model developed by Emily Moran:
   * survival probability p = intercept + size1 + size2 + size3 + JulyMaxTemp + Precip + BA
   * intercept : parameter
   * size : parameter, depends on size class: 10-25, 25-50, or 50-140 cm
   * JulyMaxTemp : July maximum temperature * parameter
   * Precip : Precipitation  * parameter
   * BA : total adult basal area within 10 m * parameter
  */

class clTree;
class clTreePopulation;
class clGrid;


class clClimateSeedlingHeightSurvival : virtual public clMortalityBase {
//note: need the virtual keyword to avoid base class ambiguity.

  public:

  /**
  * Constructor.
  * @param p_oSimManager Sim Manager object.
  */
   clClimateSeedlingHeightSurvival(clSimManager *p_oSimManager);

  /**
  * Destructor.
  */
  ~clClimateSeedlingHeightSurvival();

  /**
  * Determines mortality for a tree. This starts by checking the grid cell of
  * the tree to see if survival rate has already been calculated. If so, use it.
  * If not, calculate it for this species for this grid cell. If not already
  * determined, get the adult neighborhood basal area for the cell and calculate
  * the probability of survival using the equation above. Use the random number
  * generator to decide life or death; return the result.
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
  * <li>ValidateData() is called to validate the data.</li>
  * <li>SetupGrid() is called to create the survival grid.</li>
  * </ol>
  *
  * @param p_oDoc DOM tree of parsed input tree.
  */
  void DoShellSetup(xercesc::DOMDocument *p_oDoc);

  /**
  * Does beginning of the time step setup. This clears all grid values.
  * @param p_oPop Tree population.
  */
  //void PreMortCalcs( clTreePopulation *p_oPop );

  protected:

  /**
   * Grid holding survival rate for each species. The grid name is
   * "Seedlings". It has X+1 float data
   * members, where X = the total number of species. The data member names
   * are "survival_x", for the time step survival rate (where "x" is the
   * species number), and "BAT" for the total adult neighborhood basal area.
   */
  clGrid* mp_oGrid;

  /** intercept parameter - sized number of species.*/
  float *mp_fInter;
	/** size 0 parameter - sized number of species.*/
	float *mp_fSize0;
  /** size 1 parameter - sized number of species.*/
  float *mp_fSize1;
  /** size 2 parameter - sized number of species.*/
  float *mp_fSize2;
  /** size 3 parameter - sized number of species.*/
  float *mp_fSize3;

  /** Average Temperature parameter - sized number of species.*/
  float *mp_fJulyMaxTemp;

  /** climate water deficit parameter - sized number of species.*/
  float *mp_fPrecip;

  /** basal area parameter - sized number of species.*/
  float *mp_fBA;

  /** basal area threshold parameter - sized number of species.*/
  float *mp_fBAt;

  /** Holds data member codes for the "survival_x" data members of the
   * "Climate Comp Dep Neighborhood Survival" grid. Array size is total
   * # species.*/
  short int *mp_iGridSurvivalCodes;

  /**Neighborhood search radius.*/
  float m_fRadius;

  /**Minimum sapling height. For doing neighbor searches.*/
  float m_fMinSaplingHeight;

  /**Return code for the "BAT" grid data member.*/
  short int m_iBATCode;

  /**
  * Makes sure N and D for each species does not = 0 and that the neighbor
  * search radius >= 0.
  *
  * @throws modelErr if the above conditions are not met.
  */
  void ValidateData();

    /**
  * Reads data from the parameter file.
  * @param p_oDoc DOM tree of parsed input tree.
  * @param p_oPop Tree population.
  */
  void ReadParameterFile( xercesc::DOMDocument *p_oDoc, clTreePopulation *p_oPop );

  /**
  * Gets the total adult neighborhood basal area within the specified radius
  * from a given point. Neighbors must have a DBH greater than the minimum. They
  * also cannot be dead from a disturbance event; but any trees that have a
  * dead code of "natural" are assumed to have died in the current time step
  * mortality cycle and thus should be counted.
  * @param fX X coordinate of point for which to calculate neighborhood basal
  * area
  * @param fY Y coordinate of point for which to calculate neighborhood basal
  * area
  * @param p_oPop Tree population.
  * @returns Total adult basal area in square meters.
  */
  float GetBAT (float &fX, float &fY, clTreePopulation *p_oPop);

  /**
  * Sets up the "Climate Comp Dep Neighborhood Survival" grid. This
  * ignores any maps.
  * @param p_oPop Tree population object.
  */
  //void SetupGrid( clTreePopulation *p_oPop );
};

//---------------------------------------------------------------------------

#endif /* BEHAVIORS_CLIMATESEEDLINGHEIGHTSURVIVAL_H_ */

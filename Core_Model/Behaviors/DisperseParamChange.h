/*
 * DisperseParamChange.h
 *
 *  Created on: 10 janv. 2018
 *      Author: melaine
 */

#ifndef BEHAVIORS_DISPERSEPARAMCHANGE_H_
#define BEHAVIORS_DISPERSEPARAMCHANGE_H_

#include "DisperseBase.h"
#include "DataTypes.h"



class clTreePopulation;
class clPlot;
class clTree;

/*
 * DisperseParamChange
 *
 * This class creates and disperses seeds (based on masting spatial disperse)
 * as seeds are produced, the MeanParamG1 of each species
 * is updated with a weighted mean corresponding to the
 * number of seeds produced by each tree
 */

class clDisperseParamChange : virtual public clDisperseBase {

public:

	/*
	 * Constructor
	 */
	clDisperseParamChange(clSimManager *p_oSimManager);

	/*
	 * Destructor
	 */
	~clDisperseParamChange();

	  /**
	   * Get whether a given species masted this timestep.
	   * @param iSp Species to check
	   * @return mast or nonmast.
	   */
	  mastEvent GetMastEvent(int iSp) {
	    if (mp_iIndexes[iSp] > -1) return mp_iEvent[mp_iIndexes[iSp]];
	    else return nonmast;};

	  /**
	   * Gets the number of years since the last mast for a species.
	   * @param iSp Species.
	   * @return Number of years since the last mast.
	   */
	  int GetTimestepsSinceLastMast(int iSp) {return mp_iTimestepsSinceLastMast[iSp];};

	  /**
	   * Gets the value of the masting cumulative distribution function at a
	   * particular point.
	   * @param sp Species.
	   * @param yr Year.
	   * @return Value of masting CDF.
	   */
	  float GetMastCDF(int sp, int yr) {return mp_fMastCDF[sp][yr];};

	  /**
	    * Get max timesteps.
	    * @return Max timesteps.
	    */
	   int GetMaxTimesteps() {return m_iMaxTimesteps;};


	 /**
	  * Captures the behavior name passed from the parameter file.  This is useful
	  * since this class can produce a few different kinds of behaviors.
	  * @param sNameString Behavior name from parameter file.
	  */
	  //void SetNameData(std::string sNameString);


protected:

	   /**Cumulative probability array for seed dispersal. Array size is mast or
	    * non-mast by # species by max distance. The 0th bucket will never be
	    * accessed.*/
	   float ***mp_fSeedCDF;

	   /**STR mean (or STR value if deterministic). The array is 2D - the first
	    * index is masting or non-masting. The second index is species. This value
	    * comes from the parameter file.*/
	   float **mp_fStrMean;

	   /**STR standard deviation. The array is 2D - the first index is masting or
	    * non-masting. The second index is species. This value comes from the
	    * parameter file.*/
	   float **mp_fStrStdDev;

	   /**Beta parameter. The array is 2D - the first index is masting or
	    * non-masting. The second index is species. This value comes from the
	    * parameter file.*/
	   float **mp_fBeta;

	   /**Masting cumulative distribution function. Array size is number of species
	    * by m_iMaxYears.*/
	   float **mp_fMastCDF;

	   /** Fraction participating in disperse. Array is 2D. First index is masting
	    * or non-masting. The second index is species. This value comes from the
	    * parameter file.*/
	   float **mp_fFractionParticipating;

	   /**Whether this behavior is used by a species/type combo.  First array index
	    * is species, second is type.*/
	   bool **mp_bIsUsed;

	   /** Fecundity, if it is possible to pre-calculate it.*/
	   float *mp_fFecundity;

	   /**Array of species with each one's dbh for reproduction*/
	   float *mp_fDbhForReproduction;

	   /** Number of years since the last mast for each species */
	   int *mp_iTimestepsSinceLastMast;

	   /**This will speed access to the other arrays by storing each species' array
	    * index so the other arrays only have to be as big as the number of unique
	    * species for this behavior.*/
	   short int *mp_iIndexes;

	   /**Group affiliation of the species. Any species with the same group number
	    * will always mast together. Array size is species.*/
	   short int *mp_iGroup;

	   /** What PDF is used to draw the STR. Only deterministic, normal, and
	    * lognormal are supported.*/
	   // 4/1/2018: Zero inflated Poisson added
	   pdf *mp_iWhatPDFForSTR;

	   /**Which function is used by species*/
	   function *mp_iWhatFunction;

	   /** Whether to draw STR once per species (true) or once per tree (false).
	    * This is from the parameter file. */
	   bool *mp_bDrawSTRPerSpecies;

	   /**Which event is occurring in the current timestep for each species*/
	   enum mastEvent *mp_iEvent;

	   /**Query to perform to search for trees*/
	   char *m_cQuery;

	   /**Define a type for pointers to functions of the GetNumberOfSeeds type*/
	   typedef float (clDisperseParamChange::*Ptr2GetNumberOfSeeds)(const float &, const short int &);

	   /**Function pointer array for the appropriate function for calculating
	   * the number of seeds.  Array size is number of species to which this behavior
	   * applies.*/
	   Ptr2GetNumberOfSeeds* mp_GetSeeds;

	   /**Number of years per timestep*/
	   float m_fNumYearsPerTimestep;

	   /**Maximum distance, in meters, a seed can disperse - which is the maximum
	    * dimension of the grid with max of 1000 m*/
	   int m_iMaxDistance;

	   /** Maximum timesteps before a masting event occurs after a first masting
	    * event. This is calculated to be the value at which all species have a prob
	    * of 0.9999, or the number of timesteps in the run, whichever is smaller.*/
	   int m_iMaxTimesteps;

	   /**
	   * Does setup. This sets all values in mp_iTimestepsSinceLastMast to 1. Then
	   * it calls:
	   * <ol>
	   * <li>GetParameterFileData</li>
	   * <li>CalcMastCDF</li>
	   * <li>CalcSeedCDF</li>
	   * <li>FormatQueryString</li>
	   * <li>PopulateUsedTable</li>
	   * <li>SetGetSeedsFunctionPointers</li>
	   * </ol>
	   *
	   * @param p_oDoc DOM tree of parsed input file.
	   */
	   void DoShellSetup(xercesc::DOMDocument *p_oDoc);

	   /**
	   * Calculates the masting cumulative distribution function. This reads in the
	   * appropriate parameter values, finds the maximum number of years between
	   * mast events, then creates the CDF array. The parameters are thrown away
	   * because they are no longer needed.
	   * @param p_oDoc DOM tree of parsed input file.
	   */
	   void CalcMastCDF(xercesc::DOMDocument *p_oDoc);

	   /**
	   * Calculates the cumulative distribution functions for seed dispersal. This
	   * reads in the appropriate parameter values, finds the maximum seed dispersal
	   * distance, then creates the CDF array. The parameters are thrown away
	   * because they are no longer needed.
	   * @param p_oDoc DOM tree of parsed input file.
	   * @param p_oPop Tree population object.
	   * @throws stcErr if a weibull theta value is not less than 50 (to prevent pow
	   * overflows)
	   */
	   void CalcSeedCDF(xercesc::DOMDocument *p_oDoc, clTreePopulation *p_oPop);

	   /**
	    * Formats the string in m_cQuery. This value will be used in AddSeeds() in
	    * order to get the trees to act on.
	    * @param p_oPop Tree Population object
	    */
	   void FormatQueryString(clTreePopulation *p_oPop);

	   /**
	   * Performs dispersal of seeds for one tree. The number of seeds is calculated
	   * by using the pointer in mp_GetSeeds. Each seed is given a random azimuth
	   * direction from the parent. Then each seed is given a random distance from
	   * the parent that conforms to the chosen probability distribution function.
	   * This is done by comparing a random value to successive values in the
	   * cumulative probability array until the first array bucket that has a greater
	   * value than the random number. Once the seed has an azimuth direction and
	   * a distance, it is added to the species total in the appropriate grid cell.
	   * @param p_oTree Tree for which to perform dispersal.
	   * @param fDbh DBH of the tree, in cm.
	   * @param p_oPlot Plot object
	   * @param p_oPop Tree Population object
	   */
	   float DisperseOneParentSeeds(clTree * p_oTree, clTreePopulation * p_oPop, clPlot * p_oPlot, float fDbh );

	   /**
	   * Extracts needed parameter file data. (Some parameters are extracted by
	   * other setup functions and thrown away - this gets all other, permanent
	   * parameters.)
	   * @param p_oDoc Parsed parameter file document.
	   * @throws stcErr if:
	   * <ul>
	   * <li>The function codes are not valid enums</li>
	   * <li>A beta value is greater than 25 (to prevent pow overflows)</li>
	   * <li>A value in mp_fFractionParticipating is not between 0 and 1</li>
	   * </ul>
	   * @param p_oPop Tree population object.
	   */
	   void GetParameterFileData(xercesc::DOMDocument *p_oDoc, clTreePopulation * p_oPop);

	   /**
	   * Declares and populates the mp_bIsUsed array.
	   * @param p_oPop Tree Population object
	   */
	   void PopulateUsedTable(clTreePopulation *p_oPop);

	   /**
	    * This figures out how each species's number of seeds is calculated. It sets
	    * the appropriate function pointers in mp_GetSeeds.
	    */
	   void SetGetSeedsFunctionPointers();

	   /**
	   * Calculates the normalized probability distribution for a function
	   * for one species.
	   * @param p_fProbArray The array into which to put the normalized values.
	   * @param iMaxDistance The maximum distance out to which to calculate the
	   * function - which must equal the size of p_fProbArray.
	   * @param iFunction Function flag.
	   * @param fDispersalX0 Dispersal or X0 of the species in question, depending
	   * on the function
	   * @param fThetaXb Theta or Xb of the species in question, depending
	   * on the function
	   */
	   void CalculateProbabilityDistribution( float * p_fProbArray,
	   const int &iMaxDistance, const function &iFunction, const float &fDispersalX0,
	   const float &fThetaXb);

	   /**
	   * Performs disperse. First, DecideMast() is called to determine which species
	   * mast this timestep. Then, for any species with STRs that are drawn per
	   * species, fecundity is calculated. Then this gets the group of trees to
	   * which this behavior applies. For each, a random number is compared to the
	   * appropriate value in mp_fFractionParticipating to determine whether this
	   * tree will produce seeds. That tree is passed to DisperseOneParentSeeds()
	   * for seed creation and distribution.
	   */
	   void AddSeeds();

	   /**
	    * Decides whether masting events occur for each species. This takes the time
	    * since last mast and compares a random number to the value for that time in
	    * the mp_fMastCDF array. It then uses that decision to set the appropriate
	    * event flag in mp_iEvent. Then the counter value for each species in
	    * m_iTimestepsSinceLastMast is set appropriately; incremented if no mast
	    * occurred, or set to 0 if it did.
	    *
	    * If species have group affiliations, the mast decision is made for the
	    * first species. Then all other species in that group get the same value.
	    */
	   void DecideMast();

	   /**
	    * Gets the number of seeds for a tree when there does not need to be an
	    * STR draw.
	    * @htmlonly seeds = fecundity * DBH<sup>&beta;</sup>, @endhtmlonly
	    * Fecundity should have already been calculated.
	    * @param fDbh DBH of the tree, in cm
	    * @param iSp Tree species number
	    * @return Number of seeds.
	    */
	   float GetNumberOfSeedsNoDraw(const float &fDbh, const short int &iSp);

	   /**
	    * Gets the number of seeds for a tree when the tree needs a normal STR draw.
	    * @htmlonly seeds = STR * (DBH / 30)<sup>&beta;</sup> @endhtmlonly
	    * @param fDbh DBH of the tree, in cm
	    * @param iSp Tree species number
	    * @return Number of seeds.
	    */
	   float GetNumberOfSeedsDrawNormal(const float &fDbh, const short int &iSp);

	   /**
	    * Gets the number of seeds for a tree when the tree needs a lognormal STR
	    * draw.
	    * @htmlonly seeds = STR * (DBH / 30)<sup>&beta;</sup> @endhtmlonly
	    * @param fDbh DBH of the tree, in cm
	    * @param iSp Tree species number
	    * @return Number of seeds.
	    */
	   float GetNumberOfSeedsDrawLognormal(const float &fDbh, const short int &iSp);

	   /**
	    * Gets the number of seeds for a tree when the tree needs a zero-inflated Poisson STR
	    * draw.
	    * @htmlonly seeds = STR <sup>&beta;</sup> @endhtmlonly
	    * @param fDbh DBH of the tree, in cm
	    * @param iSp Tree species number
	    * @return Number of seeds.
	    */
	   float GetNumberOfSeedsDrawZIP(const float &fDbh, const short int &iSp);

	 };
	 //---------------------------------------------------------------------------
	 #endif

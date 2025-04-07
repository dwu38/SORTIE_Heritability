//---------------------------------------------------------------------------
// LinearGrowth
//---------------------------------------------------------------------------
#if !defined(LinearGrowth_H)
#define LinearGrowth_H

#include "GrowthBase.h"

class clGrid;

/**
 * Increments growth according to a simple linear equation.  This behavior
 * can be used to create a growth increment with no automatic height adjustment
 * or a growth increment with automatic height adjustment.
 *
 * The equation used in this behavior is:
 * <center>Y = a + b * diam</center>
 * where
 * <ul>
 * <li>Y = diameter growth in cm/year</li>
 * <li>a = growth intercept</li>
 * <li>b = growth slope</li>
 * <li>diam = diam10 for seedlings and saplings, or DBH for adults</li>
 * </ul>
 *
 *
 * The name string is "lineargrowthshell".  In the parameter file:  For
 * diameter growth with no automatic height adjustment, call
 * "LinearGrowth diam only".  For diameter growth with automatic height
 * adjustment, call "LinearGrowth".
 *
 * Copyright 2011 Charles D. Canham.
 * @author Lora E. Murphy
 *
 * <br>Edit history:
 * <br>-----------------
 * <br>October 20, 2011 - Wiped the slate clean for SORTIE 7.0 (LEM)
 *
 * 7/26/2016 Created from LinearBilevelGrowth MAK
 *
 */
class clLinearGrowth: virtual public clGrowthBase {
//note: need the virtual keyword to avoid base class ambiguity.

public:

	/**
	 * Constructor.  Sets the namestring.
	 */
	clLinearGrowth(clSimManager *p_oSimManager);

	/**
	 * Destructor.  Frees memory.
	 */
	~clLinearGrowth();

	/**
	 * Calculates the amount of diameter growth increase for a particular tree
	 * using the linear growth equation.
	 *
	 * @param p_oTree Tree for which to calculate growth.
	 * @param p_oPop Tree population object, just in case it's needed.
	 * @param fHeightGrowth Amount of height growth, in m (ignored).
	 * @return Amount, in cm, by which to increase the tree's diameter.
	 */
	float CalcDiameterGrowthValue(clTree *p_oTree, clTreePopulation *p_oPop,
			float fHeightGrowth);

	/**
	 * Does the setup for this behavior.  This reads in the parameters from the
	 * parameter file, and retrieves the "Storm Light" grid if present.
	 * @param p_oDoc Parsed parameter file.
	 * @throws modelErr if the high-light growth threshold is not between 0 and
	 * 100.
	 */
	void DoShellSetup(xercesc::DOMDocument *p_oDoc);

	/**
	 * Captures the namestring passed to this behavior.  This is overridden from
	 * clBehaviorBase so we can capture the namestring passed.  Since this class
	 * can create multiple kinds of behaviors that function differently, this will
	 * capture what kind of behavior this is supposed to be.
	 *
	 * @param sNameString Behavior's namestring.
	 */
	void SetNameData(std::string sNameString);

protected:

	/**Slope of growth equation - b - sized number of behavior
	 * species*/
	float *mp_fSlope;

	/**Intercept of growth equation - a - sized number of behavior
	 * species*/
	float *mp_fIntercept;

	/**To help access the other arrays*/
	int *mp_iIndexes;

	/**Conversion factor to translate the results of the function to the
	 * appropriate units per timestep*/
	float m_fYearsPerTimestep;
};
//---------------------------------------------------------------------------
#endif // SimpleLinearGrowth_H

//---------------------------------------------------------------------------
// SizeDepLogisticGrowth
//---------------------------------------------------------------------------
#if !defined(SizeDepLogisticGrowth_H)
#define SizeDepLogisticGrowth_H

#include "GrowthBase.h"

/**
 * Increments growth according to a size dependent logistic equation.  This can
 * be used to create a growth increment with no automatic height adjustment, a
 * growth increment with automatic height adjustment, or a height increment.
 *
 * The equation used in this behavior is:
 * <center>Y = (a + b * diam) / (1 + e<sup>c - d*GLI</sup>)</center>
 * where
 * <ul>
 * <li>Y = either radial growth in mm/year or height growth in cm/yr
 * <li>a = growth intercept
 * <li>b = growth slope
 * <li>c = shape parameter 1
 * <li>d = shape parameter 2
 * <li>diam = diameter of the tree at which to apply growth, in cm
 * <li>GLI is global light index, from 0 to 100
 * </ul>
 *
 * This behavior compounds this result for multi-year timesteps.  In the case
 * of diameter growth, a copy of the diameter value is incremented X times (with
 * X being the number of years per timestep) and the final increment is the sum
 * of all the interim increments.  In the case of height growth, the amount of
 * diameter growth is divided by X and this value used to increment height X
 * times, with the nth diameter (n = 1 to X) being the original diameter plus
 * (n * diameter growth / X).
 *
 * The name string is "sizedeplogisticgrowthshell". In the parameter file: For
 * diameter growth with no automatic height adjustment, call
 * "SizeDependentLogisticGrowth diam only". For diameter growth with automatic
 * height adjustment, call "SizeDependentLogisticGrowth".  For height growth,
 * call "SizeDependentLogisticGrowth height".
 *
 * Copyright 2011 Charles D. Canham.
 * @author Lora E. Murphy
 *
 * <br>Edit history:
 * <br>-----------------
 * <br>October 20, 2011 - Wiped the slate clean for SORTIE 7.0 (LEM)
 */
class clSizeDepLogisticGrowth: virtual public clGrowthBase {
//note: need the virtual keyword to avoid base class ambiguity.

public:

	/**
	 * Constructor.  Sets the namestring.
	 */
	clSizeDepLogisticGrowth(clSimManager *p_oSimManager);

	/**
	 * Destructor.  Frees memory.
	 */
	~clSizeDepLogisticGrowth();

	/**
	 * Calculates the amount of height growth increase for a particular tree using
	 * the size dependent logistic growth equation.  The function value is
	 * calculated for each year of the timestep by using the original diameter
	 * plus fDiameterGrowth divided by the years per timestep times the loop
	 * counter for the year being compounded.  The final increment is the sum of
	 * all the intermediate increments.
	 *
	 * @param p_oTree Tree for which to calculate growth.
	 * @param p_oPop Tree population object, just in case it's needed.
	 * @param fDiameterGrowth Amount of diameter growth for this tree, in cm.
	 * @return Amount, in m, by which to increase the tree's height.
	 */
	float CalcHeightGrowthValue(clTree *p_oTree, clTreePopulation *p_oPop,
			float fDiameterGrowth);

	/**
	 * Calculates the amount of diameter growth increase for a particular tree
	 * using the size dependent logistic growth equation.  The function value is
	 * calculated for each year of the timestep by repeatedly incrementing a copy
	 * of the diameter (the tree's actual diameter remains untouched, as it is
	 * supposed to).  The final increment is the sum of all the intermediate
	 * increments.
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
	 * parameter file, and validates that all species/type combos use light (each
	 * must have "Light" registered).
	 * @param p_oDoc Parsed parameter file.
	 * @throws modelErr if any species/type combo to which this behavior is
	 * applied does not have a light behavior.
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

	/**Growth slope - b - sized number of behavior species*/
	float *mp_fSlope;

	/**Growth intercept - a - sized number of behavior species*/
	float *mp_fIntercept;

	/**Shape parameter 1 - c - sized number of behavior species*/
	float *mp_fShape1;

	/**Shape parameter 2 - d - sized number of behavior species*/
	float *mp_fShape2;

	/**For accessing the other arrays*/
	short int *mp_iIndexes;

	/**Conversion factor to translate the results of the function to the
	 * appropriate units, depending on the type of growth behavior this is*/
	float m_fConversionFactor;

	/**Number of years per timestep*/
	float m_fNumberOfYearsPerTimestep;

	/**
	 * Calculates the value of the size dependent logistic growth function for one
	 * increment.  The meaning of what is returned depends on the type of growth
	 * the behavior is doing.
	 * @param iSpecies Tree species.
	 * @param fGLI The GLI value.
	 * @param fDiam The diameter for which to calculate this increment.
	 * @return The value of the size dependent logistic growth function.  Units
	 * depend on the type of growth that this is.
	 */
	inline float CalculateFunctionValue(int iSpecies, float fGLI, float fDiam);
};
//---------------------------------------------------------------------------
#endif // SizeDepLogisticGrowth_H

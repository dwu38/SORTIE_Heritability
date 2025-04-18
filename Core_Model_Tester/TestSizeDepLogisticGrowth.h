//---------------------------------------------------------------------------
// TestSizeDepLogisticGrowth
//---------------------------------------------------------------------------
#if !defined(TestSizeDepLogisticGrowth_H)
#define TestSizeDepLogisticGrowth_H

/**
 * Tests the clSizeDepLogisticGrowth class. This tests it in all versions - as a
 * diameter-only incrementer with no height updating, as a DiameterIncrementer
 * with auto-height updating, and as a HeightIncrementer.
 */

/**
 * Writes a parameter file to direct testing. Timestep length is 1 year.
 * This contains two species.
 * Behaviors:
 * <ul>
 * <li>size dependent logistic growth diam only - applied to species 0
 * seedlings and species 1 seedlings
 * <li>size dependent logistic growth height only - applied to species 0
 * saplings and species 1 seedlings
 * <li>constant radial growth - applied to species 0 saplings
 * <li>HeightIncrementer - applied to species 0 seedlings
 * <li> size dependent logistic growth - applied to species 1 saplings
 * </ul>
 *
 * @return Filename written.
 */
const char* WriteSizeDepLogisticGrowthXMLFile1();

/**
 * Writes a parameter file to direct testing. Timestep length is 3 years.
 *  This contains two species.
 * Behaviors:
 * <ul>
 * <li>size dependent logistic growth diam only - applied to species 0
 * seedlings and species 1 seedlings
 * <li>size dependent logistic growth height only - applied to species 0
 * saplings and species 1 seedlings
 * <li>constant radial growth - applied to species 0 saplings
 * <li>HeightIncrementer - applied to species 0 seedlings
 * <li> size dependent logistic growth - applied to species 1 saplings
 * </ul>
 *
 * @return Filename written.
 */
const char* WriteSizeDepLogisticGrowthXMLFile2();
//---------------------------------------------------------------------------
#endif // TestSizeDepLogisticGrowth_H

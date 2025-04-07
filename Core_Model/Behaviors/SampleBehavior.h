#ifndef SampleBehaviorH
#define SampleBehaviorH

#include "BehaviorBase.h"

/**
 The slash and two asterisks, closed by an asterisk and a single slash, denotes a comment for automatic documentation by Doxygen. The Doxygen web site has a manual, and there are examples throughout this header file and in existing SORTIE code for how to comment in Doxygen format. (See the SORTIE help topic for more on SORTIE documentation.) Doxygen's C++ output for SORTIE is here, if you want to see what it looks like. Put this comment immediately before your class declaration.

 Now to document what this behavior actually does. This behavior calculates biomass for adult trees according to the following equation:
 B = a * DBH d

 where:

 B = biomass in metric tons
 a = slope of biomass equation
 DBH = tree's DBH, in cm
 d = exponent of biomass equation

 Biomass is totaled in a grid object created by this behavior, called "Biomass Map". Each cell in biomass map contains the average biomass for all adult trees in that cell.

 Each behavior (in fact, each population and grid as well) has a name string. This string is a unique identifier. Other objects can use the name string to find and access this one. So document your name string where it can be seen. If you don't want other objects to access your code, don't skip the name string; just seal up your data. This behavior's name string will be "sample biomass behavior". Call it in the parameter file (if it were real) with the string "sample biomass behavior".

 @author Lora E. Murphy
 Copyright 2005 Charles D. Canham

 Edit history:
 -------------------
 February 16, 2005: Created (LEM)
 */

class clSampleBehavior: virtual public clBehaviorBase {
//note: you need the virtual keyword to avoid base class ambiguity.
public:

	/**
	 Constructor. (You must have a constructor, if only to initialize the base class constructors. Again, note the Doxygen commenting.)

	 This constructor sets all array pointers to NULL. Doing this avoids access violations if this object is destroyed before the arrays are allocated.

	 @param p_oSimManager Sim Manager object. (The "@param" is the Doxygen format for commenting a function argument. "@return" documents a function's return value, if applicable. "@throws" documents errors thrown by the function. Remember those, and the fact that the special comment goes right before the thing commented in the header file, and you know almost everything you need to to make your documentation Doxygen-compliant.)
	 */
	clSampleBehavior(clSimManager *p_oSimManager);

	/**
	 Destructor. (You are responsible for freeing your behavior's memory here. This function is not required if you don't need it.) The sample behavior frees the array memory here.
	 */
	~clSampleBehavior();

	/**
	 (This function is overridden from clBehaviorBase. It automatically executes once per timestep, at this behavior's turn in the behavior order.)

	 What this function does:

	 Clears the biomass map of the previous timestep's values by setting the "avg" and "count" data members to zero for the "Biomass Map" grid (grid documentation below).
	 Uses the clTreePopulation's Find() method to gain access to all adult trees.
	 For each adult tree, calculates its biomass and adds it to "avg", and increments "count" by one.
	 For each cell in the biomass map, divide the value in "avg" by "count" and assign back to "avg".

	 */
	void Action();

protected:

	/**
	 (This function is overridden from clBehaviorBase. It automatically executes at the beginning of the run. Behaviors can assume that the tree population and grid objects have already been set up from the parameter file when this function is called. All behavior setup code should go here, or in functions called from this function.

	 The p_oDoc parameter is the parsed parameter file. Each behavior is responsible for extracting what it needs from the parameter file, and throwing errors if it gets bad data. The functions defined in ParsingFunctions.h are there to help you read data from the parameter file.)

	 This function performs setup for this behavior. It reads in parameters from the parameter file and validates them. Then it sets up the "Biomass Map" grid.

	 @param p_oDoc DOM tree of the parsed input file.
	 @throws modelErr if a value for a is less than 0.
	 */
	void GetData(DOMDocument *p_oDoc);

	/**
	 Biomass map grid. This behavior will create this grid. Its name is "Biomass Map". It has two data members: the first is a float value called "avg" which will be used to first total the biomass values in the cell and then hold the average biomass. The second member is an int value called "count" which holds the number of adult trees in the cell, for calculating the average. The user can enter a new grid cell resolution if they want to.
	 */
	clGrid *mp_oBiomassGrid;

	/**Slope of the biomass equation (a). This is an array sized number of species.*/
	float *mp_fSlopeA;

	/**Exponent of the biomass equation (d). This is an array sized number of species.*/
	float *mp_fExpD;

	/**This is the data return code for the "avg" data member of the "Biomass Map" grid.*/
	short int m_iAvgCode;

	/**This is the data return code for the "count" data member of the "Biomass Map" grid.*/
	short int m_iCountCode;
};
#endif

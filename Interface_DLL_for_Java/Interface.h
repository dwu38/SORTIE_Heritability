//---------------------------------------------------------------------------

#ifndef InterfaceH
#define InterfaceH

//This section controls dll builds.  When using this file in a dll, the dll
//project should define BUILDING_INTERFACE.  When building an app that uses the
//resulting DLL and may include this as a header file, the app should define
//GENERATING_DLL_APP.  If using this file in a non-DLL app, don't define
//anything.
#if defined(BUILDING_INTERFACE)
#define __INTERFACE_EXPORT_TYPE __declspec(dllexport)
#elif defined(GENERATING_DLL_APP)
#define __INTERFACE_EXPORT_TYPE __declspec(dllimport)
#else
#define __INTERFACE_EXPORT_TYPE
#endif

#include "Messages.h"
class clSimManager;

/****************************************************************************
*****************************************************************************
  JAVA INTERFACE

  This DLL contains everything for sending brief messages to the Java GUI.

****************************************************************************
****************************************************************************/

//Model states - keep synched with the sim manager
enum simState{No_Data, //object managers are created but have no objects -
                          //system is ready to receive data
                 Initialized, //sufficient input data has been received - a
                        //run is possible - object managers have objects and all
                        //have a status of ready
                 Paused, //this indicates that there is an unfinished run
                 Run_Complete //this indicates that a run is complete
                };

/****************************************************************************
  REQUIRED FUNCTIONS
  These functions are those which are required to be defined and implemented
  in order to make a complete SORTIE interface.
****************************************************************************/

/*//////////////////////////////////////////////////////////////////////////
 DoIntroduction()
 This will be called before anything else is done in a non-batch situation.
 It allows the interface to call any introductory screens, or ignore it if
 there's nothing to be done.
//////////////////////////////////////////////////////////////////////////*/
__INTERFACE_EXPORT_TYPE void DoIntroduction();

/*//////////////////////////////////////////////////////////////////////////
 GetMessage()
 This is used when the model is expecting a user message and cannot go on
 without one.  This function will continue to listen for input until some has
 been received.
//////////////////////////////////////////////////////////////////////////*/
__INTERFACE_EXPORT_TYPE modelMsg GetMessage();

/**
 CheckForMessage()
 This is used when the model is checking for a message which may have come
 in during processing.  A return message of "no message" is acceptable if
 the user has not requested any actions.

 This checks for a file called "messages.txt" in the application's home
 directory.  A lack of file, or an empty file, signifies no messages.  A
 character of 'p' as the first character means "pause".  A character of 'q' as
 the first character means "quit".  Any others mean no messages.  The contents
 of the file are erased afterwards.
 @param cAppPath The path to the SORTIE application, with a path separator at
 the end (this is the expected location of messages.txt).  This string should
 be no longer than MAX_FILENAME_SIZE.
*/
__INTERFACE_EXPORT_TYPE modelMsg CheckForMessage(char *cAppPath);

/*//////////////////////////////////////////////////////////////////////////
 SendMessage()
 This function accepts a message from the sim manager.
//////////////////////////////////////////////////////////////////////////*/
__INTERFACE_EXPORT_TYPE void SendMessage(modelMsg msg, bool batch = false);

/*//////////////////////////////////////////////////////////////////////////
 ExternalErrorHandler()
 This function accepts an error from the sim manager and processes it by
 displaying the appropriate error message.  If this is batch, different
 messages can be written.
//////////////////////////////////////////////////////////////////////////*/
__INTERFACE_EXPORT_TYPE void ExternalErrorHandler(
                    modelErr stcError, int iModelState, bool batch = false);

/***************************************************************************
  INTERFACE-SPECIFIC FUNCTIONS
  These functions are specific to this interface.
****************************************************************************/

/*///////////////////////////////////////////////////////////////////////////
  WriteIntroductoryScreen()
  This writes the screen that appears when SORTIE is first started from the
  command line.
///////////////////////////////////////////////////////////////////////////*/
void WriteIntroductoryScreen();

/*///////////////////////////////////////////////////////////////////////////
  WriteBadCommand()
  This writes an error message if the user types a command that the model
  doesn't recognize.
///////////////////////////////////////////////////////////////////////////*/
void WriteBadCommand();

/*///////////////////////////////////////////////////////////////////////////
  WriteGeneralHelp()
  This writes the general help screen.
///////////////////////////////////////////////////////////////////////////*/
void WriteGeneralHelp();

/*///////////////////////////////////////////////////////////////////////////
  WriteInputHelp()
  This writes the input help screen.
///////////////////////////////////////////////////////////////////////////*/
void WriteInputHelp();

/*///////////////////////////////////////////////////////////////////////////
  WriteRunHelp()
  This writes the help screen for the "run" command.
///////////////////////////////////////////////////////////////////////////*/
void WriteRunHelp();

/*///////////////////////////////////////////////////////////////////////////
  WritePauseHelp()
  This writes the help screen for the "pause" command.
///////////////////////////////////////////////////////////////////////////*/
void WritePauseHelp();

/*///////////////////////////////////////////////////////////////////////////
  WriteQuitHelp()
  This writes the help screen for the "quit" command.
///////////////////////////////////////////////////////////////////////////*/
void WriteQuitHelp();

/*///////////////////////////////////////////////////////////////////////////
  WriteNeedFile()
  If the application is expecting a file name and path and doesn't get one,
  this prompts the user for one
///////////////////////////////////////////////////////////////////////////*/
void WriteNeedFile();

/*///////////////////////////////////////////////////////////////////////////
  WriteBadFile()
  If a user inputs a bad file name, this gives an appropriate error message
///////////////////////////////////////////////////////////////////////////*/
void WriteBadFile(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteModelReady()
  This writes a message indicating that the model is ready.
///////////////////////////////////////////////////////////////////////////*/
void WriteModelReady();

/*///////////////////////////////////////////////////////////////////////////
  WriteModelNotReady()
  This writes a message indicating that the model is not ready.
///////////////////////////////////////////////////////////////////////////*/
void WriteModelNotReady();

/*///////////////////////////////////////////////////////////////////////////
  WriteBadArgument()
  This writes that there is an argument that is not understood.
///////////////////////////////////////////////////////////////////////////*/
void WriteBadArgument(modelMsg stcMsg);

/*///////////////////////////////////////////////////////////////////////////
  WriteNotReadyToRun()
  If the user has attempted to run the model when the model was not ready,
  this writes an error message.
///////////////////////////////////////////////////////////////////////////*/
void WriteNotReadyToRun();

/*///////////////////////////////////////////////////////////////////////////
  WriteRunComplete()
  This writes that the requested run is complete.
///////////////////////////////////////////////////////////////////////////*/
void WriteRunComplete(int iCurrentTimestep, int iNumTimesteps, unsigned long
  iRunTime);

/*///////////////////////////////////////////////////////////////////////////
  WriteModelPaused()
  This writes that the model is paused, and what to do to start it again.
///////////////////////////////////////////////////////////////////////////*/
void WriteModelPaused(modelMsg stcMsg);

/*///////////////////////////////////////////////////////////////////////////
  WriteBadFileType()
  This writes that the file type passed is not what was expected.
///////////////////////////////////////////////////////////////////////////*/
void WriteBadFileType(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteBadXMLFile()
  This writes that the XML file passed is malformed or invalid.
///////////////////////////////////////////////////////////////////////////*/
void WriteBadXMLFile(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteDataMissingFromFile()
  This writes that data was missing from a file passed.
///////////////////////////////////////////////////////////////////////////*/
void WriteDataMissingFromFile(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteBadData()
  This writes that a file passed had bad data.
///////////////////////////////////////////////////////////////////////////*/
void WriteBadData(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteIllegalOp()
  This writes that an object attempted an illegal operation.
///////////////////////////////////////////////////////////////////////////*/
void WriteIllegalOp(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteDataReadOnly()
  This writes that something tried to access read-only data.  Not a user-
  based error.
///////////////////////////////////////////////////////////////////////////*/
void WriteDataReadOnly(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteCantFindObject()
  This writes that a needed object wasn't found.  Probably not a user-based
  error.
///////////////////////////////////////////////////////////////////////////*/
void WriteCantFindObject(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteTreeWrongType()
  This writes that a tree was not the expected type.  Probably not a user-
  based error.
///////////////////////////////////////////////////////////////////////////*/
void WriteTreeWrongType(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteAccessViolation()
  This writes that an access violation has occurred.  Not a user-based error.
///////////////////////////////////////////////////////////////////////////*/
void WriteAccessViolation(modelErr stcErr);

/*///////////////////////////////////////////////////////////////////////////
  WriteUnknownError()
  This writes that an unknown error has occurred.
///////////////////////////////////////////////////////////////////////////*/
void WriteUnknownError(modelErr stcErr);

//---------------------------------------------------------------------------
#endif



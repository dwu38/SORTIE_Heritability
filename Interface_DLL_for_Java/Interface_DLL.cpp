#include <windows.h>
extern "C"
{
#ifdef __BORLANDC__
  #pragma argsused
#endif
  int WINAPI DllMain( HINSTANCE hInst, /* Library instance handle. */
  unsigned long reason, /* Reason this function is being called. */
  void * lpReserved ) /* Not used. */
  {
    return 1;
  }
  int WINAPI DllEntryPoint( HINSTANCE hInst, unsigned long reason, void * lpReserved )
  {
    return DllMain( hInst, reason, lpReserved );
  }
}

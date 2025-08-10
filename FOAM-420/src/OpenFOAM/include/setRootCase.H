Foam::argList args(argc, argv);
if (!args.checkRootCase())
{
    Foam::FatalError.exit();
}
#if defined(WIN64) || defined(WIN32)
//Hack by blueCAPE Lda 2009
#include "include/forceLoadLibraries.H"
#endif
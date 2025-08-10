    argList::validOptions.insert("latestTime", "");
    argList::validOptions.insert("time", "time");

    argList::addOption
    (
        "dict",
        "file",
        "specify alternative dictionary for the caseSetup"
    );

    argList::addOption
    (
        "writePause",
        "scalar",
        "time delay for file system to sync (default = 0.0 s)"
    );

    argList::addBoolOption
    (
        "initialiseFields",
        "perform only field initialisation"
    );


#include "include/setRootCase.H"

    scalar writePause = args.optionLookupOrDefault<scalar>("writePause", 0.0);
    bool initialiseFields = false;
    if (args.optionFound("initialiseFields"))
    {
        initialiseFields = true;
    }

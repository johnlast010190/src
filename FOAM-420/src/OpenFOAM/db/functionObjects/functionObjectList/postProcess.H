/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2016 OpenFOAM Foundation

Global
    postProcess

Description
    Execute application functionObjects to post-process existing results.

    If the "dict" argument is specified the functionObjectList is constructed
    from that dictionary otherwise the functionObjectList is constructed from
    the "functions" sub-dictionary of "system/controlDict"

    Multiple time-steps may be processed and the standard utility time
    controls are provided.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef CREATE_TIME
    #define CREATE_TIME include/createTime.H
#endif

#ifndef CREATE_MESH
    #define CREATE_MESH include/createMesh.H
#endif

#ifndef CREATE_FIELDS
    #define CREATE_FIELDS createFields.H
#endif

#ifndef CREATE_CONTROL
    #define CREATE_CONTROL cfdTools/general/solutionControl/createControl.H
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define INCLUDE_FILE(X) INCLUDE_FILE2(X)
#define INCLUDE_FILE2(X) #X

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::argList::addBoolOption
(
    argList::postProcessOptionName,
    "Execute functionObjects only"
);

if (argList::postProcess(argc, argv))
{
    Foam::timeSelector::addOptions();
    #include "include/addRegionOption.H"
    #include "include/addFunctionObjectOptions.H"

    // Set functionObject post-processing mode
    functionObject::postProcess = true;

    #include "include/setRootCase.H"

    if (args.optionFound("list"))
    {
        functionObjectList::list();
        return 0;
    }

    const bool lowMem  = args.optionFound("lowMemory");

    #include INCLUDE_FILE(CREATE_TIME)
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include INCLUDE_FILE(CREATE_MESH)

    #ifndef NO_CONTROL
    #include INCLUDE_FILE(CREATE_CONTROL)
    #endif

    // Externally stored dictionary for functionObjectList
    // if not constructed from runTime
    dictionary functionsDict;

    HashSet<word> selectedFields;

    // Construct functionObjectList
    autoPtr<functionObjectList> functionsPtr;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        FatalIOError.throwExceptions();

        try
        {
            #include INCLUDE_FILE(CREATE_FIELDS)

            #ifdef CREATE_FIELDS_2
            #include INCLUDE_FILE(CREATE_FIELDS_2)
            #endif

            #ifdef CREATE_FIELDS_3
            #include INCLUDE_FILE(CREATE_FIELDS_3)
            #endif

            if
            (
                functionsPtr.empty() || lowMem
             || mesh.readUpdate() != polyMesh::UNCHANGED
            )
            {
                // Update functionObjects if mesh changes
                functionsPtr = functionObjectList::New
                (
                    args,
                    runTime,
                    functionsDict,
                    selectedFields
                );
            }

            functionsPtr->execute(lowMem);

            // Execute the functionObject 'end()' function for the last time
            if (!lowMem && timei == timeDirs.size()-1)
            {
                functionsPtr->end();
            }
        }
        catch (IOerror& err)
        {
            Warning<< err << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef INCLUDE_FILE
#undef INCLUDE_FILE2

#undef CREATE_MESH
#undef CREATE_FIELDS
#undef CREATE_CONTROL

// ************************************************************************* //

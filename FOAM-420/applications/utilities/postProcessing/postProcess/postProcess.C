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

Application
    postProcess

Description
    Execute the set of functionObjects specified in the selected dictionary
    (which defaults to system/controlDict) or on the command-line for the
    selected set of times on the selected set of fields.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/timeSelector.H"
#include "fields/ReadFields/ReadFields.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define ReadFields(GeoFieldType)                                               \
    readFields<GeoFieldType>(mesh, objects, selectedFields, storedObjects);

#define ReadPointFields(GeoFieldType)                                          \
    readFields<GeoFieldType>(pMesh, objects, selectedFields, storedObjects);

#define ReadUniformFields(FieldType)                                           \
    readUniformFields<FieldType>                                               \
    (constantObjects, selectedFields, storedObjects);

void executeFunctionObjects
(
    const argList& args,
    const Time& runTime,
    fvMesh& mesh,
    const HashSet<word>& selectedFields,
    functionObjectList& functions,
    bool lastTime,
    bool lowMem
)
{
    Info<< nl << "Reading fields:" << endl;

    // Maintain a stack of the stored objects to clear after executing
    // the functionObjects
    LIFOStack<regIOobject*> storedObjects;

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read volFields
    ReadFields(volScalarField);
    ReadFields(volVectorField);
    ReadFields(volSphericalTensorField);
    ReadFields(volSymmTensorField);
    ReadFields(volTensorField);

    // Read internal fields
    ReadFields(volScalarField::Internal);
    ReadFields(volVectorField::Internal);
    ReadFields(volSphericalTensorField::Internal);
    ReadFields(volSymmTensorField::Internal);
    ReadFields(volTensorField::Internal);

    // Read surface fields
    ReadFields(surfaceScalarField);
    ReadFields(surfaceVectorField);
    ReadFields(surfaceSphericalTensorField);
    ReadFields(surfaceSymmTensorField);
    ReadFields(surfaceTensorField);

    // Read point fields.
    const pointMesh& pMesh = pointMesh::New(mesh);

    ReadPointFields(pointScalarField)
    ReadPointFields(pointVectorField);
    ReadPointFields(pointSphericalTensorField);
    ReadPointFields(pointSymmTensorField);
    ReadPointFields(pointTensorField);

    // Read uniform dimensioned fields
    IOobjectList constantObjects(mesh, runTime.constant());

    ReadUniformFields(uniformDimensionedScalarField);
    ReadUniformFields(uniformDimensionedVectorField);
    ReadUniformFields(uniformDimensionedSphericalTensorField);
    ReadUniformFields(uniformDimensionedSymmTensorField);
    ReadUniformFields(uniformDimensionedTensorField);

    Info<< nl << "Executing functionObjects" << endl;

    // Execute the functionObjects in post-processing mode
    functions.execute(lowMem);

    // Execute the functionObject 'end()' function for the last time
    if (!lowMem && lastTime)
    {
        functions.end();
    }

    while (!storedObjects.empty())
    {
        storedObjects.pop()->checkOut();
    }
}


int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    #include "include/addRegionOption.H"
    #include "include/addFunctionObjectOptions.H"

    argList::addBoolOption
    (
        "lowMemory",
        "Delete function object after executed to reduce maximum memory"
    );

    // Set functionObject post-processing mode
    functionObject::postProcess = true;

    #include "include/setRootCase.H"

    if (args.optionFound("list"))
    {
        functionObjectList::list();
        return 0;
    }

    #include "include/createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include "include/createNamedMesh.H"

    // Initialize the set of selected fields from the command-line options
    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }
    if (args.optionFound("field"))
    {
        selectedFields.insert(args.optionLookup("field")());
    }

    const bool lowMem  = args.optionFound("lowMemory");

    // Externally stored dictionary for functionObjectList
    // if not constructed from runTime
    dictionary functionsDict;

    // Construct functionObjectList
    autoPtr<functionObjectList> functionsPtr;

    if (!lowMem)
    {
        functionsPtr = functionObjectList::New
        (
            args,
            runTime,
            functionsDict,
            selectedFields
        );
    }

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        if (mesh.readUpdate() != polyMesh::UNCHANGED || lowMem)
        {
            // Update functionObjectList if mesh changes
            functionsPtr = functionObjectList::New
            (
                args,
                runTime,
                functionsDict,
                selectedFields
            );
        }

        FatalIOError.throwExceptions();

        try
        {
            executeFunctionObjects
            (
                args,
                runTime,
                mesh,
                selectedFields,
                functionsPtr(),
                timei == timeDirs.size()-1,
                lowMem
            );
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


// ************************************************************************* //

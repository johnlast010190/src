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
    (c) 2010-2016, 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "stateIndex/stateIndex.H"
#include "global/etcFiles/etcFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void createMesh(const Time& runTime, autoPtr<fvMesh>& mesh)
{
    if (!mesh.valid())
    {
        mesh.reset
        (
            new fvMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            )
        );
    }
}

void updateFvOptionEntry
(
    dictionary& sourceDict,
    dictionary& targetDict,
    const word& entryName
)
{
    if (sourceDict.found(entryName))
    {
        targetDict.add
        (
            entryName,
            sourceDict.lookup(entryName)
        );

        sourceDict.remove(entryName);
    }
}

int main(int argc, char *argv[])
{
#include "addOptions.H"
#include "include/createDictDefaults.H"
#include "include/createTime.H"
#include "checkTimeOptions.H"
    runTime.functionObjects().off();

    bool collated = false;
    if (Foam::fileHandler().type() != "uncollated")
    {
        collated = true;
    }

#include "setupDict.H"


    // Construct stateIndex and compile entries
    stateIndex caseState(runTime, dict, defaults, distributed, collated);

    if (initialiseFields)
    {
        //When initialise fields read the BCs from CBC dict.
        caseState.resetBoundaries(false);

        caseState.createRegionObjects();

        caseState.createFields();
    }
    else
    {
        caseState.writeDictionaries();

        caseState.createRegionObjects();

        caseState.modifyMesh(writePause);

        caseState.createFields();
    }

    Info<< nl << nl << "end" << nl << endl;

    return(0);
}

// ************************************************************************* //


// ************************************************************************* //

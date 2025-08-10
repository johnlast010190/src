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
    (c) 2019-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/coordinateFrame.H"
#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::coordinateFrame& Foam::coordinateFrame::New
(
    const fvMesh& mesh,
    const word& frameName
)
{
    if (!mesh.foundObject<coordinateFrame>(frameName))
    {
        const word dictName = "meshObjects";

        autoPtr<IOobject> meshObjectIO
        (
            new IOobject
            (
                dictName,
                mesh.time().caseSystem(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        if (!meshObjectIO->typeHeaderOk<localIOdictionary>())
        {
            // Fallback: look for definition at the top level
            autoPtr<IOobject> defaultRegionIO
            (
                new IOobject
                (
                    dictName,
                    mesh.time().caseSystem(),
                    mesh.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            if (defaultRegionIO->typeHeaderOk<localIOdictionary>())
            {
                meshObjectIO = defaultRegionIO;
            }
        }

        // Use objectPath rather than filePath so a useful error is given if
        // file is not found
        IFstream is(meshObjectIO->objectPath());
        dictionary meshObjDict(is);
        dictionary frameDict = meshObjDict.subDict(frameName);

        const auto ctor =
            ctorTableLookup
            (
                "coordinate system type",
                dictionaryConstructorTable_(),
                frameDict.lookup<word>("type")
            );
        regIOobject::store(ctor(mesh, frameDict, frameName).ptr());
    }

    return mesh.lookupObjectRef<coordinateFrame>(frameName);
}


Foam::coordinateFrame* Foam::coordinateFrame::New
(
    const fvMesh& mesh,
    const word& frameName,
    const dictionary& dict
)
{
    if (!mesh.foundObject<coordinateFrame>(frameName))
    {
        const auto ctor =
            ctorTableLookup
            (
                "coordinate system type",
                dictionaryConstructorTable_(),
                dict.lookup<word>("type")
            );
        regIOobject::store(ctor(mesh, dict, frameName).ptr());
    }

    return &mesh.lookupObjectRef<coordinateFrame>(frameName);
}


Foam::coordinateFrame* Foam::coordinateFrame::lookupNew
(
    const fvMesh& mesh,
    const dictionary& dictWithFrameName
)
{
    return
        &coordinateFrame::New
        (
            mesh,
            dictWithFrameName.lookup<word>("referenceFrame")
        );
}


Foam::coordinateFrame* Foam::coordinateFrame::globalFrame(const fvMesh& mesh)
{
    dictionary frameDict;
    frameDict.add("type", "coordinateFrame");
    dictionary coorSys;
    coorSys.add("type", "cartesian");
    coorSys.add("origin", "(0 0 0)");
    coorSys.add("e1", "(1 0 0)");
    coorSys.add("e2", "(0 1 0)");
    frameDict.set("coordinateSystem", coorSys);
    return coordinateFrame::New(mesh, "globalFrame", frameDict);
}


// ************************************************************************* //

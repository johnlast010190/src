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

\*---------------------------------------------------------------------------*/

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fvMeshFunctionObject, 0);
}
}

Foam::fileName Foam::functionObjects::fvMeshFunctionObject::outputFileDir() const
{
    fileName baseDir = mesh_.time().path();

    if (Pstream::parRun())
    {
        // Put in undecomposed case (Note: gives problems for
        // distributed data running)
        baseDir = baseDir/".."/outputPrefix;
    }
    else
    {
        baseDir = baseDir/outputPrefix;
    }

    // Append mesh name if not default
    if (isA<polyMesh>(mesh_))
    {
        const polyMesh& mesh = dynamic_cast<const polyMesh&>(mesh_);
        if (mesh.name() != polyMesh::defaultRegion)
        {
            baseDir = baseDir/mesh.name();
        }
    }

    // Append region name if different
    if (obr_.name() != word::null && obr_.name() != mesh_.name())
    {
        baseDir = baseDir/obr_.name();
    }

    // Remove any ".."
    baseDir.clean();

    return baseDir;
}


const Foam::fvMesh&
Foam::functionObjects::fvMeshFunctionObject::selectMesh
(
    const Time& runTime,
    const objectRegistry& obr
) const
{
    word meshName;

    // The dictionary operatingPointDict only exists when obr behaves as an
    // operating point of the multi-point of FOAM-Adjoint.
    if (obr.foundObject<IOdictionary>("operatingPointDict"))
    {
        meshName =
            obr.lookupObject<IOdictionary>("operatingPointDict").
                lookup<word>
                (
                    "parentMeshName"
                );
    }
    else
    {
        meshName = meshObr_.name();
    }

    return runTime.lookupObject<fvMesh>(meshName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fvMeshFunctionObject::fvMeshFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    mesh_(selectMesh(runTime, obr_))
{}


Foam::functionObjects::fvMeshFunctionObject::fvMeshFunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    regionFunctionObject(name, obr, dict),
    mesh_(refCast<const fvMesh>(obr_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fvMeshFunctionObject::~fvMeshFunctionObject()
{}


// ************************************************************************* //

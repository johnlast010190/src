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
    (c) 2015 OpenCFD Ltd.
    (c) 2014-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "decompositionModel.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionModel, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionModel::decompositionModel
(
    const polyMesh& mesh,
    const fileName& decompDictFile
)
:
    MeshObject
    <
        polyMesh,
        Foam::UpdateableMeshObject,
        decompositionModel
    >(mesh),
    IOdictionary
    (
        selectIO
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh.local(),
#if defined(WIN64) || defined(WIN32)
                mesh,
                #else
                mesh.db(),
                #endif
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false,  //io.registerObject(),
                true    //io.globalObject()
            ),
            decompDictFile
        )
    )
{}


Foam::decompositionModel::decompositionModel
(
    const polyMesh& mesh,
    const dictionary& dict,
    const fileName& decompDictFile
)
:
    MeshObject
    <
        polyMesh,
        Foam::UpdateableMeshObject,
        decompositionModel
    >(mesh),
    IOdictionary
    (
        selectIO
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh.local(),
#if defined(WIN64) || defined(WIN32)
                mesh,
                #else
                mesh.db(),
                #endif
                (dict.size() ? IOobject::NO_READ : IOobject::MUST_READ),
                IOobject::NO_WRITE,
                false,  //io.registerObject(),
                true    //io.globalObject()
            ),
            decompDictFile
        ),
        dict
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

const Foam::decompositionModel& Foam::decompositionModel::New
(
    const polyMesh& mesh,
    const fileName& decompDictFile
)
{
    return
        MeshObject
        <
            polyMesh,
            Foam::UpdateableMeshObject,
            decompositionModel
        >::New(mesh, decompDictFile);
}


const Foam::decompositionModel& Foam::decompositionModel::New
(
    const polyMesh& mesh,
    const dictionary& dict,
    const fileName& decompDictFile
)
{
    return
        MeshObject
        <
            polyMesh,
            Foam::UpdateableMeshObject,
            decompositionModel
        >::New(mesh, dict, decompDictFile);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::IOobject Foam::decompositionModel::selectIO
(
    const IOobject& io,
    const fileName& f
)
{
    return
    (
        f.size()
      ? IOobject        // construct from filePath instead
        (
            fileName(f).toAbsolute(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            io.registerObject(),
            io.globalObject()
        )
     :  io
    );
}


// ************************************************************************* //

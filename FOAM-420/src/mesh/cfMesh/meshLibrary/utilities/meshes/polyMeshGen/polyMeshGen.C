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
    (c) Creative Fields, Ltd.
    (c) 2017 Esi Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description

\*---------------------------------------------------------------------------*/

#include "utilities/meshes/polyMeshGen/polyMeshGen.H"
#include "include/demandDrivenData.H"
#include "db/IOstreams/Fstreams/OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGen::polyMeshGen(const Time& t)
:
    polyMeshGenCells(t),
    metaDict_
    (
        IOobject
        (
            "meshMetaDict",
            runTime_.constant(),
            "polyMesh",
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{}

//- Construct from components without the boundary
polyMeshGen::polyMeshGen
(
    const Time& t,
    const pointField& points,
    const faceList& faces,
    const cellList& cells
)
:
    polyMeshGenCells(t, points, faces, cells),
    metaDict_
    (
        IOobject
        (
            "meshMetaDict",
            runTime_.constant(),
            "polyMesh",
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{}

//- Construct from components with the boundary
polyMeshGen::polyMeshGen
(
    const Time& t,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const polyBoundaryMesh& patches
)
:
    polyMeshGenCells
    (
        t,
        points,
        faces,
        cells,
        patches
    ),
    metaDict_
    (
        IOobject
        (
            "meshMetaDict",
            runTime_.constant(),
            "polyMesh",
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor
polyMeshGen::~polyMeshGen()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGen::read()
{
    polyMeshGenCells::read();
}

void polyMeshGen::write() const
{
    //- remove old mesh before writting
    const fileName meshDir = runTime_.path()/runTime_.constant()/"polyMesh";

    rm(meshDir/"points");
    rm(meshDir/"faces");
    rm(meshDir/"owner");
    rm(meshDir/"neighbour");
    rm(meshDir/"cells");
    rm(meshDir/"boundary");
    rm(meshDir/"pointZones");
    rm(meshDir/"faceZones");
    rm(meshDir/"cellZones");
    rm(meshDir/"meshModifiers");
    rm(meshDir/"parallelData");
    rm(meshDir/"meshMetaDict");

    // remove sets if they exist
    if (isDir(meshDir/"sets"))
    {
        rmDir(meshDir/"sets");
    }

    //- write the mesh
    polyMeshGenCells::write();

    //- write meta data
    OFstream fName(meshDir/"meshMetaDict");

    metaDict_.writeHeader(fName);
    metaDict_.writeData(fName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

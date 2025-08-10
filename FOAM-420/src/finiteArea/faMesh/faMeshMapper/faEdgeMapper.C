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
    (c) 1991-2005 OpenCFD Ltd.

Description
    FV edge mapper.

\*---------------------------------------------------------------------------*/

#include "faMesh/faMeshMapper/faEdgeMapper.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faEdgeMapper::calcAddressing() const
{
    if (directAddrPtr_)
    {
        FatalErrorInFunction
            << "Addressing already calculated"
            << abort(FatalError);
    }

    // Dummy mapping: take value from edge 0
    directAddrPtr_ = new labelList(size(), 0);
}


void Foam::faEdgeMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faEdgeMapper::faEdgeMapper
(
    const faMesh& mesh,
    const mapPolyMesh& mpm
)
:
    mesh_(mesh),
    mpm_(mpm),
    sizeBeforeMapping_(mesh.nInternalEdges()),
    directAddrPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faEdgeMapper::~faEdgeMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::faEdgeMapper::directAddressing() const
{
    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::faEdgeMapper::addressing() const
{
    FatalErrorInFunction
        << "Requested interpolative addressing for a direct mapper."
        << abort(FatalError);
    ::abort();
}


const Foam::scalarListList& Foam::faEdgeMapper::weights() const
{
    FatalErrorInFunction
        << "Requested interpolative weights for a direct mapper."
        << abort(FatalError);
    ::abort();
}


// ************************************************************************* //

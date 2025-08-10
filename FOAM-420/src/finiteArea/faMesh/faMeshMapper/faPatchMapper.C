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

\*---------------------------------------------------------------------------*/

#include "faMesh/faMeshMapper/faPatchMapper.H"
#include "faMesh/faPatches/faPatch/faPatch.H"
#include "faMesh/faBoundaryMesh/faBoundaryMesh.H"
#include "faMesh/faMesh.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "meshes/polyMesh/mapPolyMesh/faceMapper/faceMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faPatchMapper::calcAddressing() const
{
    if (directAddrPtr_)
    {
        FatalErrorInFunction
            << "Addressing already calculated"
            << abort(FatalError);
    }

    directAddrPtr_ = new labelList(patch_.size(), 0);
    labelList& addr = *directAddrPtr_;

    hasUnmapped_ = false;

    // Make a map of old edgeFaces, giving edge index in patch given the new
    // face label next to the patch

    // Create edge index lookup
    Map<label> edgeIndexLookup;

    const labelList& reverseFaceMap = mpm_.reverseFaceMap();

    forAll(oldEdgeFaces_, oefI)
    {
        if (reverseFaceMap[oldEdgeFaces_[oefI]] > -1)
        {
            // Face has survived.  Insert its label under new face index
            edgeIndexLookup.insert(reverseFaceMap[oldEdgeFaces_[oefI]], oefI);
        }
    }

    // Go through new edgeFaces and for each edge try to locate old index
    const labelList& ef = patch_.edgeFaces();

    forAll(ef, efI)
    {
        if (edgeIndexLookup.found(ef[efI]))
        {
            addr[efI] = edgeIndexLookup[ef[efI]];
        }
        else
        {
            // Not found: map from zero
            addr[efI] = 0;
        }
    }
}


void Foam::faPatchMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    hasUnmapped_ = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faPatchMapper::faPatchMapper
(
    const faPatch& patch,
    const mapPolyMesh& mpm
)
:
    patch_(patch),
    mpm_(mpm),
    sizeBeforeMapping_(patch.size()),
    oldEdgeFaces_(patch.edgeFaces()),
    hasUnmapped_(false),
    directAddrPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faPatchMapper::~faPatchMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::faPatchMapper::directAddressing() const
{
    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::faPatchMapper::addressing() const
{
    FatalErrorInFunction
        << "Requested interpolative addressing for a direct mapper."
        << abort(FatalError);
    ::abort();
}


const Foam::scalarListList& Foam::faPatchMapper::weights() const
{
    FatalErrorInFunction
        << "Requested interpolative weights for a direct mapper."
        << abort(FatalError);
    ::abort();
}


// ************************************************************************* //

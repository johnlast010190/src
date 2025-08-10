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
#include "faMesh/faMeshMapper/faMeshMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshMapper::faMeshMapper
(
    const faMesh& mesh,
    const mapPolyMesh& mpm
)
:
    mesh_(mesh),
    nOldPoints_(mesh.nPoints()),
    nOldEdges_(mesh.nEdges()),
    nOldInternalEdges_(mesh.nInternalEdges()),
    nOldFaces_(mesh.nFaces()),
    oldPatchSizes_(mesh.boundary().size(), 0),
    oldPatchStarts_(mesh.boundary().size(), -1),
    oldPatchEdgeFaces_(mesh.boundary().size()),
    areaMap_(mesh, mpm),
    edgeMap_(mesh, mpm),
    boundaryMap_(mesh, mpm)
{
    // Capture old patch information
    const faBoundaryMesh& patches = mesh.boundary();

    forAll(patches, patchI)
    {
        oldPatchSizes_[patchI] = patches[patchI].size();
        oldPatchStarts_[patchI] = patches[patchI].start();

        oldPatchEdgeFaces_[patchI] = patches[patchI].edgeFaces();
    }
}


// ************************************************************************* //

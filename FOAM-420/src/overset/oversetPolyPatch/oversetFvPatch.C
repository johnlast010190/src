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
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "oversetPolyPatch/oversetFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvMesh.H"
#include "primitives/transform/transform.H"
#include "cellCellStencil/cellCellStencil/cellCellStencilObject.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, oversetFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::oversetFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::oversetFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& restrictMap
) const
{
    // Store the restrictMap. This routine gets used for:
    // - GAMGAgglomeration: this is the use we want to intercept.
    // - GAMGProcAgglomeration: to find out the cell number on the other side.
    if (master())
    {
        restrictMap_ = restrictMap;
    }

    return patchInternalField(restrictMap);
}


const Foam::labelListList& Foam::oversetFvPatch::stencil() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());

    return overlap.cellStencil();
}


const Foam::mapDistribute& Foam::oversetFvPatch::cellInterpolationMap() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());
    return overlap.cellInterpolationMap();
}


const Foam::List<Foam::scalarList>&
Foam::oversetFvPatch::cellInterpolationWeights() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());
    return overlap.cellInterpolationWeights();
}


const Foam::scalarField& Foam::oversetFvPatch::normalisation() const
{
    return boundaryMesh().mesh().V().field();
}


const Foam::labelList& Foam::oversetFvPatch::interpolationCells() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());
    return overlap.interpolationCells();
}


const Foam::scalarList& Foam::oversetFvPatch::cellInterpolationWeight() const
{
    const cellCellStencilObject& overlap = Stencil::New(boundaryMesh().mesh());
    return overlap.cellInterpolationWeight();
}


// ************************************************************************* //

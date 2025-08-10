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
    (c) 2011-2016 OpenFOAM Foundation

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/primitiveMesh.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellCentresAndVols() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Calculating cell centres and cell volumes"
            << endl;
    }

    // It is an error to attempt to recalculate cellCentres
    // if the pointer is already set
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorInFunction
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(nCells());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(nCells());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Finished calculating cell centres and cell volumes"
            << endl;
    }
}


void Foam::primitiveMesh::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    // Clear the fields for accumulation
    cellCtrs = Zero;
    cellVols = 0.0;

    const labelList& owner = faceOwner();
    const labelList& neighbour = faceNeighbour();

    // first estimate the approximate cell centre as the average of
    // face centres

    vectorField cEst(nCells(), Zero);
    scalarField cellVolsMag(nCells(), Zero);
    labelField nCellFaces(nCells(), 0);

    forAll(owner, facei)
    {
        label own = owner[facei];
        cEst[own] += fCtrs[facei];
        nCellFaces[own] += 1;
    }

    forAll(neighbour, facei)
    {
        label nei = neighbour[facei];
        cEst[nei] += fCtrs[facei];
        nCellFaces[nei] += 1;
    }

    forAll(cEst, celli)
    {
        cEst[celli] /= nCellFaces[celli];
    }

    Switch version22
    (
        debug::debugSwitches().lookupOrDefault<Switch>
        ("version22", false)
    );

    forAll(owner, facei)
    {
        label own = owner[facei];
        // Calculate 3*face-pyramid volume
        scalar  pyr3Vol = 0;
        if (!version22)
        {
            pyr3Vol = fAreas[facei] & (fCtrs[facei] - cEst[own]);
        }
        else
        {
            pyr3Vol =
                max(fAreas[facei] & (fCtrs[facei] - cEst[own]), VSMALL);
        }

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst[own];

        scalar pyr3VolMag = mag(pyr3Vol);

        // Accumulate volume-weighted face-pyramid centre
        cellCtrs[own] += pyr3VolMag*pc;

        // Accumulate face-pyramid volume
        cellVols[own] += pyr3Vol;
        cellVolsMag[own] += pyr3VolMag;
    }

    forAll(neighbour, facei)
    {
        label nei = neighbour[facei];
        // Calculate 3*face-pyramid volume
        scalar  pyr3Vol = 0;
        if (!version22)
        {
            pyr3Vol =
                fAreas[facei] & (cEst[nei] - fCtrs[facei]);
        }
        else
        {
            pyr3Vol =
                max(fAreas[facei] & (cEst[nei] - fCtrs[facei]), VSMALL);
        }

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst[nei];

        scalar pyr3VolMag = mag(pyr3Vol);

        // Accumulate volume-weighted face-pyramid centre
        cellCtrs[nei] += pyr3VolMag*pc;

        // Accumulate face-pyramid volume
        cellVols[nei] += pyr3Vol;
        cellVolsMag[nei] += pyr3VolMag;
    }

    if (!version22)
    {
        forAll(cellCtrs, celli)
        {
            if (mag(cellVolsMag[celli]) > VSMALL)
            {
                cellCtrs[celli] /= cellVolsMag[celli];
            }
            else
            {
                cellCtrs[celli] = cEst[celli];
            }
        }
    }
    else
    {
        cellCtrs /= cellVolsMag;
    }

    cellVols *= (1.0/3.0);
}

// ************************************************************************* //

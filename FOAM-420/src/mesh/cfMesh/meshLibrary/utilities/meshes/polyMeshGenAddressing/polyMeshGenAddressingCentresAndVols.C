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

Class
    polyMeshGenAddressing

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressing.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcCellCentresAndVols() const
{
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorInFunction
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    const cellListPMG& cells = mesh_.cells();

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(cells.size());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(cells.size());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);
}

void polyMeshGenAddressing::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    const labelList& own = mesh_.owner();
    const cellListPMG& cells = mesh_.cells();
    const label nCells = cells.size();

    # ifdef USE_OMP
    # pragma omp parallel for if (nCells > 1000)
    # endif
    for (label cellI=0;cellI<nCells;++cellI)
    {
        const cell& c = cells[cellI];

        //- approximate the centre first
        vector cEst(vector::zero);
        forAll(c, fI)
            cEst += fCtrs[c[fI]];

        cEst /= c.size();

        //- start evaluating the volume and the cell centre
        vector cellCentre(vector::zero);
        scalar cellVol(0.0);
        scalar cellVolMag(0.0);

        forAll(c, fI)
        {
            // Calculate 3*face-pyramid volume
            scalar pyr3Vol = (fAreas[c[fI]] & (fCtrs[c[fI]] - cEst));

            if (own[c[fI]] != cellI)
                pyr3Vol *= -1.0;

            scalar pyr3VolMag = mag(pyr3Vol);

            // Calculate face-pyramid centre
            const vector pc = (3.0/4.0)*fCtrs[c[fI]] + (1.0/4.0)*cEst;

            // Accumulate volume-weighted face-pyramid centre
            cellCentre += pyr3VolMag*pc;

            // Accumulate face-pyramid volume
            cellVol += pyr3Vol;
            cellVolMag += pyr3VolMag;
        }

        if (mag(cellVolMag) > VSMALL)
        {
            cellCtrs[cellI] = cellCentre / cellVolMag;
        }
        else
        {
            cellCtrs[cellI] = cEst;
        }

        cellVols[cellI] = cellVol / 3.0;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& polyMeshGenAddressing::cellCentres() const
{
    if (!cellCentresPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcCellCentresAndVols();
    }

    return *cellCentresPtr_;
}

const scalarField& polyMeshGenAddressing::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcCellCentresAndVols();
    }

    return *cellVolumesPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

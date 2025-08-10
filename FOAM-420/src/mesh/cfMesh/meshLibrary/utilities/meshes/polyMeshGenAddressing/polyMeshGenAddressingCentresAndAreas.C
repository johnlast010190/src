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

Description
    Calulate the face centres and areas.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

\*---------------------------------------------------------------------------*/

#include "utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcFaceCentresAndAreas() const
{
    if (faceCentresPtr_ || faceAreasPtr_)
    {
        FatalErrorInFunction
            << "Face centres or face areas already calculated"
            << abort(FatalError);
    }

    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    faceCentresPtr_ = new vectorField(faces.size());
    vectorField& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorField(faces.size());
    vectorField& fAreas = *faceAreasPtr_;

    makeFaceCentresAndAreas(points, fCtrs, fAreas);
}

void polyMeshGenAddressing::makeFaceCentresAndAreas
(
    const pointFieldPMG& p,
    vectorField& fCtrs,
    vectorField& fAreas
) const
{
    const faceListPMG& fs = mesh_.faces();
    const label nFaces = fs.size();

    # ifdef USE_OMP
    # pragma omp parallel for if (nFaces > 1000)
    # endif
    for (label facei=0;facei<nFaces;++facei)
    {
        face f = fs[facei];
        label nPoints = f.size();

        bool flip = false;

        //Need face centres to be calculated in consistent order at processor
        //faces to avoid rounding error differences (particularly in SP).
        //Flip the face based on a geometric test for calculation.
        if (facei >= mesh_.nInternalFaces())
        {
            point pNext = p[f[1]];
            point pPrev = p[f[nPoints-1]];
            vector diff = pNext-pPrev;

            forAll(diff, cmpt)
            {
                if (diff[cmpt] > SMALL)
                {
                    flip = true;
                    f.flip();
                    break;
                }
                else if (diff[cmpt] < -SMALL)
                {
                    break;
                }
            }
        }

        // If the face is a triangle, do a direct calculation for efficiency
        // and to avoid round-off error-related problems
        if (nPoints == 3)
        {
            fCtrs[facei] = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
            if (flip)
            {
                fAreas[facei] = -0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
            }
            else
            {
                fAreas[facei] = 0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
            }
        }
        else
        {
            vector sumN = Zero;
            scalar sumA = 0.0;
            vector sumAc = Zero;

            point fCentre = Zero;
            forAll(f, pI)
            {
                fCentre += p[f[pI]];
            }

            fCentre /= nPoints;

            forAll(f, pI)
            {
                const point& pt = p[f[pI]];
                const point& pNext = p[f[f.fcIndex(pI)]];

                const vector n = (pNext - pt)^(fCentre - pt);
                sumN += n;
            }
            const vector sumNHat = normalised(sumN);

            forAll(f, pI)
            {
                const point& pt = p[f[pI]];
                const point& pNext = p[f[f.fcIndex(pI)]];

                const vector n = (pNext - pt)^(fCentre - pt);
                const vector c = pt + pNext + fCentre;
                scalar a = n & sumNHat;

                sumA += a;
                sumAc += a*c;
            }

            // This is to deal with zero-area faces. Mark very small faces
            // to be detected in e.g., processorPolyPatch.
            if (sumA > VSMALL)
            {
                fCtrs[facei] = (1.0/3.0)*sumAc/sumA;
            }
            else
            {
                fCtrs[facei] = fCentre;
            }

            if (flip)
            {
                fAreas[facei] = -0.5*sumN;
            }
            else
            {
                fAreas[facei] = 0.5*sumN;
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& polyMeshGenAddressing::faceCentres() const
{
    if (!faceCentresPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcFaceCentresAndAreas();
    }

    return *faceCentresPtr_;
}

const vectorField& polyMeshGenAddressing::faceAreas() const
{
    if (!faceAreasPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcFaceCentresAndAreas();
    }

    return *faceAreasPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

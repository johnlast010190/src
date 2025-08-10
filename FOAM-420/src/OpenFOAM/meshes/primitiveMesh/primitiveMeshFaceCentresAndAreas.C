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
    (c) 2017-2022 Esi Ltd.

Description
    Calulate the face centres and areas.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/primitiveMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcFaceCentresAndAreas() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcFaceCentresAndAreas() : "
            << "Calculating face centres and face areas"
            << endl;
    }

    // It is an error to attempt to recalculate faceCentres
    // if the pointer is already set
    if (faceCentresPtr_ || faceAreasPtr_ || magFaceAreasPtr_)
    {
        FatalErrorInFunction
            << "Face centres or face areas already calculated"
            << abort(FatalError);
    }

    faceCentresPtr_ = new vectorField(nFaces());
    vectorField& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorField(nFaces());
    vectorField& fAreas = *faceAreasPtr_;

    magFaceAreasPtr_ = new scalarField(nFaces());
    scalarField& magfAreas = *magFaceAreasPtr_;

    makeFaceCentresAndAreas(points(), fCtrs, fAreas, magfAreas);

    if (debug)
    {
        Pout<< "primitiveMesh::calcFaceCentresAndAreas() : "
            << "Finished calculating face centres and face areas"
            << endl;
    }
}


void Foam::primitiveMesh::makeFaceCentresAndAreas
(
    const pointField& p,
    vectorField& fCtrs,
    vectorField& fAreas,
    scalarField& magfAreas
) const
{
    const faceList& fs = faces();

    forAll(fs, facei)
    {
        const face& f = fs[facei];
        label nPoints = f.size();

        bool flip = false;

        //Need face centres to be calculated in consistent order at processor
        //faces to avoid rounding error differences (particularly in SP).
        //Flip the face based on a geometric test for calculation.

        if (facei >= nInternalFaces())
        {
            point pNext = p[f[1]];
            point pPrev = p[f[nPoints-1]];
            vector diff = pNext-pPrev;

            forAll(diff, cmpt)
            {
                if (diff[cmpt] > SMALL)
                {
                    flip = true;
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
            if (flip)
            {
                fCtrs[facei] = (1.0/3.0)*(p[f[0]] + p[f[2]] + p[f[1]]);
                fAreas[facei] = -0.5*((p[f[2]] - p[f[0]])^(p[f[1]] - p[f[0]]));
            }
            else
            {
                fCtrs[facei] = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
                fAreas[facei] = 0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
            }
        }
        else
        {
            vector sumN = Zero;
            scalar sumA = 0.0;
            vector sumAc = Zero;

            point fCentre = Zero;

            if (flip)
            {
                label nextFp(0);
                forAll(f, pI)
                {
                    fCentre += p[f[nextFp]];
                    nextFp = f.rcIndex(nextFp);
                }

                fCentre /= nPoints;

                nextFp = 0;
                forAll(f, pI)
                {
                    const point& pt = p[f[nextFp]];
                    nextFp = f.rcIndex(nextFp);
                    const point& pNext = p[f[nextFp]];

                    const vector n = (pNext - pt)^(fCentre - pt);
                    sumN += n;
                }
                const vector sumNHat = normalised(sumN);
                nextFp = 0;
                forAll(f, pI)
                {
                    const point& pt = p[f[nextFp]];
                    nextFp = f.rcIndex(nextFp);
                    const point& pNext = p[f[nextFp]];

                    const vector n = (pNext - pt)^(fCentre - pt);
                    const vector c = pt + pNext + fCentre;
                    scalar a = n & sumNHat;

                    sumA += a;
                    sumAc += a*c;
                }
            }
            else
            {
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

        magfAreas[facei] = max(mag(fAreas[facei]), VSMALL);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// TODO: Inline me
const Foam::scalarField& Foam::primitiveMesh::magFaceAreas() const
{
    if (!magFaceAreasPtr_)
    {
        calcFaceCentresAndAreas();
    }

    return *magFaceAreasPtr_;
}

// ************************************************************************* //

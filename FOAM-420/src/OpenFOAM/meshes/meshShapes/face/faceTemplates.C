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
    (c) 2011-2016,2020 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/meshShapes/face/face.H"
#include "containers/Lists/DynamicList/DynamicList.H"

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

template<class PointField>
Foam::point Foam::face::faceCentre(const PointField& ps)
{
    // If the face is a triangle, do a direct calculation
    if (ps.size() == 3)
    {
        return (1.0/3.0)*(ps[0] + ps[1] + ps[2]);
    }

    // For more complex faces, decompose into triangles ...

    // Compute an estimate of the centre as the average of the points
    point pAvg = Zero;
    forAll(ps, pi)
    {
        pAvg += ps[pi];
    }
    pAvg /= ps.size();

    // Compute the face area normal and unit normal by summing up the
    // normals of the triangles formed by connecting each edge to the
    // point average.
    vector sumA = Zero;
    forAll(ps, pi)
    {
        const point& p = ps[pi];
        const point& pNext = ps[ps.fcIndex(pi)];

        const vector a = (pNext - p)^(pAvg - p);

        sumA += a;
    }
    const vector sumAHat = normalised(sumA);

    // Compute the area-weighted sum of the triangle centres. Note use
    // the triangle area projected in the direction of the face normal
    // as the weight, *not* the triangle area magnitude. Only the
    // former makes the calculation independent of the initial estimate.
    scalar sumAn = 0;
    vector sumAnc = Zero;
    forAll(ps, pi)
    {
        const point& p = ps[pi];
        const point& pNext = ps[ps.fcIndex(pi)];

        const vector a = (pNext - p)^(pAvg - p);
        const vector c = p + pNext + pAvg;

        const scalar an = a & sumAHat;

        sumAn += an;
        sumAnc += an*c;
    }

    // Complete calculating centres and areas. If the face is too small
    // for the sums to be reliably divided then just set the centre to
    // the initial estimate.
    if (sumAn > VSMALL)
    {
        return (1.0/3.0)*sumAnc/sumAn;
    }
    else
    {
        return pAvg;
    }
}


template<class PointField>
Foam::vector Foam::face::faceArea(const PointField& ps)
{
    // If the face is a triangle, do a direct calculation
    if (ps.size() == 3)
    {
        return 0.5*((ps[1] - ps[0])^(ps[2] - ps[0]));
    }

    // For more complex faces, decompose into triangles ...

    // Compute an estimate of the centre as the average of the points
    point pAvg = Zero;
    forAll(ps, pi)
    {
        pAvg += ps[pi];
    }
    pAvg /= ps.size();

    // Compute the face area normal and unit normal by summing up the
    // normals of the triangles formed by connecting each edge to the
    // point average.
    vector sumA = Zero;
    forAll(ps, pi)
    {
        const point& p = ps[pi];
        const point& pNext = ps[ps.fcIndex(pi)];

        const vector a = (pNext - p)^(pAvg - p);

        sumA += a;
    }

    return 0.5*sumA;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::face::average
(
    const UList<point>& meshPoints,
    const Field<Type>& fld
) const
{
    // Calculate the average by breaking the face into triangles and
    // area-weighted averaging their averages

    // If the face is a triangle, do a direct calculation
    if (size() == 3)
    {
        return
            (1.0/3.0)
           *(
               fld[operator[](0)]
             + fld[operator[](1)]
             + fld[operator[](2)]
            );
    }

    // For more complex faces, decompose into triangles ...
    label nPoints = size();

    // Compute an estimate of the centre and the field average as the average
    // of the point values
    point centrePoint = Zero;
    Type cf = Zero;
    forAll(*this, pI)
    {
        centrePoint += meshPoints[operator[](pI)];
        cf += fld[operator[](pI)];
    }

    centrePoint /= nPoints;
    cf /= nPoints;

    // Compute the face area normal and unit normal by summing up the
    // normals of the triangles formed by connecting each edge to the
    // point average.
    vector sumN = Zero;
    forAll(*this, pI)
    {
        const point& p = meshPoints[operator[](pI)];
        const point& pNext = meshPoints[operator[](fcIndex(pI))];

        const vector n = (pNext - p)^(centrePoint - p);

        sumN += n;
    }
    const vector sumNHat = normalised(sumN);

    // Compute the area-weighted sum of the average field values on each
    // triangle
    scalar sumA = 0;
    Type sumAf = Zero;

    forAll(*this, pI)
    {
        const point& p = meshPoints[operator[](pI)];
        const point& pNext = meshPoints[operator[](fcIndex(pI))];

        const vector n = (pNext - p)^(centrePoint - p);

        // Calculate 3*triangle centre field value
        Type ttcf  =
        (
            fld[operator[](pI)]
          + fld[operator[](fcIndex(pI))]
          + cf
        );

        // Calculate 2*triangle area
        scalar ta = n & sumNHat;

        sumA += ta;
        sumAf += ta*ttcf;
    }

    if (sumA > VSMALL)
    {
        return (1.0/3.0)*sumAf/sumA;
    }
    else
    {
        return cf;
    }
}


// ************************************************************************* //

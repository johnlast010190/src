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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class PointField>
void Foam::triIntersect::writePolygon
(
    const word& name,
    const PointField& ps
)
{
    pointField pf(ps.size());
    forAll(ps, pointi)
    {
        pf[pointi] = ps[pointi];
    }

    Info<< indent << "Writing face to " << name + ".vtk" << endl;

    simpleVTKWriter polyVTK(labelList(identity(ps.size())), pf);
    polyVTK.write(name + ".vtk");
}


template<class Type>
Type Foam::triIntersect::srcTriInterpolate
(
    const barycentric2D& y,
    const FixedList<Type, 3>& srcPsis
)
{
    const label yMini = findMin(y);

    // Inside the triangle
    if (y[yMini] >= 0)
    {
        return y[0]*srcPsis[0] + y[1]*srcPsis[1] + y[2]*srcPsis[2];
    }

    barycentric2D z(y);
    z[yMini] = VGREAT;
    const label yMidi = findMin(z);

    // Outside an edge
    if (y[yMidi] >= 0)
    {
        const label yi0 = (yMini + 1) % 3, yi1 = (yi0 + 1) % 3;
        return (y[yi0]*srcPsis[yi0] + y[yi1]*srcPsis[yi1])/(y[yi0] + y[yi1]);
    }

    const label yMaxi = findMax(y);

    // Outside a corner
    return srcPsis[yMaxi];
}


template<class Type>
Type Foam::triIntersect::tgtTriInterpolate
(
    const barycentric2D& y,
    const FixedList<Type, 3>& tgtPsis
)
{
    return y.a()*tgtPsis[0] + y.b()*tgtPsis[1] + y.c()*tgtPsis[2];
}


// ************************************************************************* //

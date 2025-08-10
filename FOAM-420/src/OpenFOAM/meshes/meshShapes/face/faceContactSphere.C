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
    (c) 2011 OpenFOAM Foundation

Description
    Return location of contact sphere on the face

\*---------------------------------------------------------------------------*/

#include "meshes/meshShapes/face/face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::face::contactSphereDiameter
(
    const point& p,
    const vector& n,
    const UList<point>& meshPoints
) const
{
    scalar magN = Foam::mag(n);

    vector n1 = n/(magN + SMALL);
    vector n2 = areaNormal(meshPoints);

    n2 /= Foam::mag(n2) + SMALL;

    return 2*((centre(meshPoints) - p) & n2)/((n1 & n2) - 1.0);
}


// ************************************************************************* //

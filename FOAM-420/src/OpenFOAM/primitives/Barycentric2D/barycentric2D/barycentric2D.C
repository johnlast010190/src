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
    (c) 2017-2018 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/Barycentric2D/barycentric2D/barycentric2D.H"
#include "primitives/random/Random/Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::barycentric2D Foam::barycentric2D01(Random& rndGen)
{
    scalar s = rndGen.sample01<scalar>();
    scalar t = rndGen.sample01<scalar>();

    // Transform the random point in the unit square to a random point in the
    // unit tri by reflecting across the diagonal if needed
    if (s + t > 1)
    {
        s = 1 - s;
        t = 1 - t;
    }

    return Foam::barycentric2D(1 - s - t, s, t);
}


// ************************************************************************* //

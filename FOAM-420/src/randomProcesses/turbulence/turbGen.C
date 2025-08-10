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

\*---------------------------------------------------------------------------*/

#include "fft/fft.H"
#include "turbulence/turbGen.H"
#include "Kmesh/Kmesh.H"
#include "fields/Fields/primitiveFields.H"
#include "turbulence/Ek.H"
#include "global/constants/mathematical/mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbGen::turbGen(const Kmesh& k, const scalar EA, const scalar K0)
:
    K(k),
    Ea(EA),
    k0(K0),
    RanGen(label(0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vectorField Foam::turbGen::U()
{
    vectorField s(K.size());
    scalarField rndPhases(K.size());

    forAll(K, i)
    {
        s[i] = RanGen.sample01<vector>();
        rndPhases[i] = RanGen.sample01<scalar>();
    }

    s = K ^ s;
    s = s/(mag(s) + 1.0e-20);

    s = Ek(Ea, k0, mag(K))*s;

    complexVectorField up
    (
        fft::reverseTransform
        (
            ComplexField(cos(constant::mathematical::twoPi*rndPhases)*s,
            sin(constant::mathematical::twoPi*rndPhases)*s),
            K.nn()
        )
    );

    return ReImSum(up);
}


// ************************************************************************* //

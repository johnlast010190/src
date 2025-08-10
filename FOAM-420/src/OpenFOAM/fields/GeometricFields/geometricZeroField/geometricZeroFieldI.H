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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/GeometricFields/geometricZeroField/geometricZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::dimensionSet& Foam::geometricZeroField::dimensions() const
{
    return dimless;
}

inline Foam::scalar Foam::geometricZeroField::operator[](const label) const
{
    return scalar(0);
}

inline Foam::zeroField Foam::geometricZeroField::field() const
{
    return zeroField();
}

inline Foam::zeroField Foam::geometricZeroField::operator()() const
{
    return zeroField();
}

inline Foam::zeroField Foam::geometricZeroField::oldTime() const
{
    return zeroField();
}

inline Foam::zeroFieldField Foam::geometricZeroField::boundaryField() const
{
    return zeroFieldField();
}


// ************************************************************************* //

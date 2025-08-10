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

\*---------------------------------------------------------------------------*/

#include "fields/GeometricFields/geometricOneField/geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::dimensionSet& Foam::geometricOneField::dimensions() const
{
    return dimless;
}


inline Foam::one Foam::geometricOneField::operator[](const label) const
{
    return one();
}


inline Foam::oneField Foam::geometricOneField::field() const
{
    return oneField();
}


inline Foam::oneField Foam::geometricOneField::oldTime() const
{
    return oneField();
}


inline Foam::geometricOneField::Internal
Foam::geometricOneField::operator()() const
{
    return Internal();
}


inline Foam::geometricOneField::Internal
Foam::geometricOneField::v() const
{
    return Internal();
}


inline typename Foam::geometricOneField::Internal::FieldType
Foam::geometricOneField::primitiveField() const
{
    return typename Internal::FieldType();
}


inline Foam::geometricOneField::Boundary
Foam::geometricOneField::boundaryField() const
{
    return Boundary();
}


inline const Foam::geometricOneField& Foam::operator*
(
    const geometricOneField& gof,
    const geometricOneField&
)
{
    return gof;
}


inline const Foam::geometricOneField::Internal& Foam::operator*
(
    const geometricOneField::Internal& of,
    const geometricOneField&
)
{
    return of;
}


inline const Foam::geometricOneField::Internal& Foam::operator*
(
    const geometricOneField&,
    const geometricOneField::Internal& of
)
{
    return of;
}


inline const Foam::geometricOneField& Foam::operator/
(
    const geometricOneField& gof,
    const geometricOneField&
)
{
    return gof;
}


inline const Foam::geometricOneField::Internal& Foam::operator/
(
    const geometricOneField::Internal& of,
    const geometricOneField&
)
{
    return of;
}


inline const Foam::geometricOneField::Internal& Foam::operator/
(
    const geometricOneField&,
    const geometricOneField::Internal& of
)
{
    return of;
}


// ************************************************************************* //

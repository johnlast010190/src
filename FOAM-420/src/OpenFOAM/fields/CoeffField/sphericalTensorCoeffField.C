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
    (c) 2004-6 H. Jasak All rights reserved

\*---------------------------------------------------------------------------*/

#include "fields/CoeffField/sphericalTensorCoeffField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CoeffField<Foam::sphericalTensor>::CoeffField(const label size)
:
    DecoupledCoeffField<sphericalTensor>(size)
{}


Foam::CoeffField<Foam::sphericalTensor>::CoeffField
(
    const CoeffField<sphericalTensor>& f
)
:
    DecoupledCoeffField<sphericalTensor>(f)
{}


Foam::CoeffField<Foam::sphericalTensor>::CoeffField
(
    const DecoupledCoeffField<sphericalTensor>& f
)
:
    DecoupledCoeffField<sphericalTensor>(f)
{}


Foam::CoeffField<Foam::sphericalTensor>::CoeffField
(
    const tmp<DecoupledCoeffField<sphericalTensor>>& tf
)
:
    DecoupledCoeffField<sphericalTensor>(tf())
{}


Foam::CoeffField<Foam::sphericalTensor>::CoeffField(Istream& is)
:
    DecoupledCoeffField<sphericalTensor>(is)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const CoeffField<sphericalTensor>& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const tmp<CoeffField<sphericalTensor>>& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const CoeffField<sphericalTensor>::scalarTypeField& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const tmp<CoeffField<sphericalTensor>::scalarTypeField>& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const CoeffField<sphericalTensor>::linearTypeField& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const tmp<CoeffField<sphericalTensor>::linearTypeField>& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

Foam::tmp<Foam::CoeffField<Foam::sphericalTensor>> Foam::inv
(
    const CoeffField<sphericalTensor>& f
)
{
    const DecoupledCoeffField<sphericalTensor>& df = f;

    return tmp<CoeffField<sphericalTensor>>
    (
        new CoeffField<sphericalTensor>(inv(df)())
    );
}


// ************************************************************************* //

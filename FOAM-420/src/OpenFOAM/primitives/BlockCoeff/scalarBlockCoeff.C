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

#include "primitives/BlockCoeff/scalarBlockCoeff.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BlockCoeff<Foam::scalar>::BlockCoeff()
:
    scalarCoeff_(pTraits<scalar>::zero)
{}


Foam::BlockCoeff<Foam::scalar>::BlockCoeff(const BlockCoeff<scalar>& f)
:
    scalarCoeff_(f.scalarCoeff_)
{}


Foam::BlockCoeff<Foam::scalar>::BlockCoeff(Istream& is)
:
    scalarCoeff_(readScalar(is))
{}


Foam::BlockCoeff<Foam::scalar> Foam::BlockCoeff<Foam::scalar>::clone() const
{
    return BlockCoeff<scalar>(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BlockCoeff<Foam::scalar>::~BlockCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::BlockCoeff<Foam::scalar>::operator=(const BlockCoeff<scalar>& f)
{
    if (this == &f)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    scalarCoeff_ = f.scalarCoeff_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const BlockCoeff<scalar>& f)
{
    os << f.scalarCoeff_;

    return os;
}


// ************************************************************************* //

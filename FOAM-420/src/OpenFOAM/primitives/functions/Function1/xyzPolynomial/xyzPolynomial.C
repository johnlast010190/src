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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/xyzPolynomial/xyzPolynomial.H"
#include "primitives/functions/Function1/PolynomialEntry/PolynomialEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::xyzPolynomial<Type>::xyzPolynomial
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName),
    inputComp_()
{
    read(dict);
}


template<class Type>
Foam::Function1Types::xyzPolynomial<Type>::xyzPolynomial
(
    const word& entryName,
    const List<Tuple2<Type, vector>>& coeffs
)
:
    Function1<Type>(entryName),
    inputComp_(coeffs)
{}


template<class Type>
Foam::Function1Types::xyzPolynomial<Type>::xyzPolynomial
(
    const xyzPolynomial& poly
)
:
    Function1<Type>(poly),
    inputComp_(poly.inputComp_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::xyzPolynomial<Type>::~xyzPolynomial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::xyzPolynomial<Type>::convertTimeBase(const Time& t)
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::xyzPolynomial<Type>::value(const vector& xyz) const
{
    Type y(Zero);
    forAll(inputComp_, i)
    {
        scalar x = 1.0;
        const vector& xyzExp = inputComp_[i].second();
        forAll(xyzExp, compI)
        {
            x *= pow(xyz[compI], xyzExp[compI]);
        }
        y += cmptMultiply(inputComp_[i].first(), pTraits<Type>::one*x);
    }
    return y;
}


template<class Type>
Type Foam::Function1Types::xyzPolynomial<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::xyzPolynomial<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
void Foam::Function1Types::xyzPolynomial<Type>::read(const dictionary& dict)
{
    Istream& is(dict.lookup(this->name_));
    word entryType(is);
    is  >> inputComp_;
    if (!inputComp_.size())
    {
        FatalErrorInFunction
            << "Polynomial coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


template<class Type>
void Foam::Function1Types::xyzPolynomial<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << nl << indent << inputComp_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //

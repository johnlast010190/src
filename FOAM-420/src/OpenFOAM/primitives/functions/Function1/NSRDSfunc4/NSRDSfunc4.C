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
    (c) 2019-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/NSRDSfunc4/NSRDSfunc4.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc4<Type>::NSRDSfunc4
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName)
{
    Istream& is(dict.lookup(entryName));
    // Gobble type of function1
    word entryType(is);
    is >> coeffs_;
}


template<class Type>
Foam::Function1Types::NSRDSfunc4<Type>::NSRDSfunc4
(
    const word& entryName,
    const FixedList<Type, 5>& coeffs
)
:
    Function1<Type>(entryName),
    coeffs_(coeffs)
{}


template<class Type>
Foam::Function1Types::NSRDSfunc4<Type>::NSRDSfunc4(const NSRDSfunc4& poly)
:
    Function1<Type>(poly),
    coeffs_(poly.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc4<Type>::~NSRDSfunc4()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::NSRDSfunc4<Type>::indefiniteIntegral
(
    const scalar x
) const
{
    // Return the value of the vectorised indefinite integral of the function
    // a_ + b_/x + c_/pow(x, 3) + d_/pow(x, 8) + e_/pow(x, 9)
    // given by
    // a_*x + b_*ln(x) - c_/(2*pow(x, 2)) - d_/(7*pow(x, 7)) - e_/(8*pow(x, 8)

    return
        coeffs_[0]*x + coeffs_[1]*log(x) - coeffs_[2]/(2*pow(x, 2))
      - coeffs_[3]/(7*pow(x, 7)) - coeffs_[4]/(8*pow(x, 8));
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc4<Type>::indefiniteIntegralYoverX
(
    const scalar x
) const
{
    // Return the value of the vectorised indefinite integral of the function
    // a_ + b_/x + c_/pow(x, 3) + d_/pow(x, 8) + e_/pow(x, 9) divided by x, i.e.
    // a_/x + b_/pow(x,2) + c_/pow(x, 4) + d_/pow(x, 9) + e_/pow(x, 10) .
    //
    // The integral is given by
    // a_*ln(x) - b_/x - c_/(3*pow(x, 3)) - d_/(8*pow(x, 8)) - e_/(9*pow(x, 9) .

    return
        coeffs_[0]*log(x) - coeffs_[1]/x - coeffs_[2]/(3*pow(x, 3))
      - coeffs_[3]/(8*pow(x, 8)) - coeffs_[4]/(9*pow(x, 9));
}


template<class Type>
void Foam::Function1Types::NSRDSfunc4<Type>::convertTimeBase(const Time& t)
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc4<Type>::value(const scalar x) const
{
    // Evaluate the function
    // a_ + b_/x + c_/pow(x, 3) + d_/pow(x, 8) + e_/pow(x, 9)
    // and return the result

    return
        coeffs_[0] + coeffs_[1]/x + coeffs_[2]/pow(x, 3) + coeffs_[3]/pow(x, 8)
      + coeffs_[4]/pow(x, 9);
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc4<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    // Evaluate the vectorised integral of the function
    // a_ + b_/x + c_/pow(x, 3) + d_/pow(x, 8) + e_/pow(x, 9)
    // given by
    // a_*x + b_*ln(x) - c_/(2*pow(x, 2)) - d_/(7*pow(x, 7)) - e_/(8*pow(x, 8)
    // and return the result of the definite integral between x1 and x2

    return indefiniteIntegral(x2) - indefiniteIntegral(x1);
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc4<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    // Evaluate the vectorised integral of the function
    // a_ + b_/x + c_/pow(x, 3) + d_/pow(x, 8) + e_/pow(x, 9) divided by x, i.e.
    // a_/x + b_/pow(x,2) + c_/pow(x, 4) + d_/pow(x, 9) + e_/pow(x, 10) .
    //
    // The integral is given by
    // a_*ln(x) - b_/x - c_/(3*pow(x, 3)) - d_/(8*pow(x, 8)) - e_/(9*pow(x, 9) .
    //
    // The result of the definite integral between x1 and x2 is returned.

    return indefiniteIntegralYoverX(x2) - indefiniteIntegralYoverX(x1);
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc4<Type>::derivative(const scalar x) const
{
    // Evaluate the vectorised derivative of the function
    // a_ + b_/x + c_/pow(x, 3) + d_/pow(x, 8) + e_/pow(x, 9)
    // given by -b_/sqr(x) - 3*c_/pow(x, 4) - 8*d_/pow(x, 9) - 9*e_/pow(x, 10)
    // and return the result

    return
      - coeffs_[1]/sqr(x) - 3*coeffs_[2]/pow(x, 4) - 8*coeffs_[3]/pow(x, 9)
      - 9*coeffs_[4]/pow(x, 10);
}


template<class Type>
void Foam::Function1Types::NSRDSfunc4<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    // Switch to ASCII for writing
    const bool isBinary(os.format() == IOstream::BINARY);
    isBinary ? os.format(IOstream::ASCII) : os.format(IOstream::ASCII);
    os  << nl << indent << coeffs_
        << token::END_STATEMENT << nl;
    isBinary ? os.format(IOstream::BINARY) : os.format(IOstream::ASCII);
}


// ************************************************************************* //

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

#include "primitives/functions/Function1/NSRDSfunc1/NSRDSfunc1.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc1<Type>::NSRDSfunc1
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
Foam::Function1Types::NSRDSfunc1<Type>::NSRDSfunc1
(
    const word& entryName,
    const FixedList<Type, 5>& coeffs
)
:
    Function1<Type>(entryName),
    coeffs_(coeffs)
{}


template<class Type>
Foam::Function1Types::NSRDSfunc1<Type>::NSRDSfunc1(const NSRDSfunc1& poly)
:
    Function1<Type>(poly),
    coeffs_(poly.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc1<Type>::~NSRDSfunc1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::NSRDSfunc1<Type>::convertTimeBase(const Time& t)
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc1<Type>::value(const scalar x) const
{

    // Evaluate the function
    // exp(a_ + b_/x + c_*log(x) + d_*pow(x, e_))
    // and return the result

    return
        cmptExp
        (
            coeffs_[0] + coeffs_[1]/x + coeffs_[2]*log(x)
          + cmptMultiply
            (
                coeffs_[3], cmptPow(pTraits<Type>::one*x, coeffs_[4])
            )
        );
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc1<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc1<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc1<Type>::derivative(const scalar x) const
{
    // Evaluate the derivative of the function
    // exp(a_ + b_/x + c_*log(x) + d_*pow(x, e_))
    //
    // given by
    // exp(a_ + b_/x + c_*log(x) + d_*pow(x, e_))
    //      *(-b_/(x*x) + c_/x + d_*e_*pow(x, e_ - 1))
    // and return the result

    return
        cmptMultiply
        (
            cmptExp
            (
                coeffs_[0] + coeffs_[1]/x + coeffs_[2]*log(x)
              + cmptMultiply
                (
                    coeffs_[3], cmptPow(pTraits<Type>::one*x, coeffs_[4])
                )
            ),
          - coeffs_[1]/(x*x) + coeffs_[2]/x
          + cmptMultiply
            (
                coeffs_[3],
                cmptMultiply
                (
                    coeffs_[4],
                    cmptPow
                    (
                        pTraits<Type>::one*x,
                        coeffs_[4]-pTraits<Type>::one
                    )
                )
            )
        );
}


template<class Type>
void Foam::Function1Types::NSRDSfunc1<Type>::writeData(Ostream& os) const
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

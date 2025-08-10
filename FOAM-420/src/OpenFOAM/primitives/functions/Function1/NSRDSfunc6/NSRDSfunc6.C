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

#include "primitives/functions/Function1/NSRDSfunc6/NSRDSfunc6.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "db/IOstreams/Tstreams/ITstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc6<Type>::NSRDSfunc6
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName)
{
    ITstream& stream = dict.lookup(entryName);
    // Gobble type of function1
    word entryType(stream);
    stream >> Tc_;
    stream >> coeffs_;
}


template<class Type>
Foam::Function1Types::NSRDSfunc6<Type>::NSRDSfunc6
(
    const word& entryName,
    const scalar& Tc,
    const FixedList<Type, 5>& coeffs
)
:
    Function1<Type>(entryName),
    Tc_(Tc),
    coeffs_(coeffs)
{}


template<class Type>
Foam::Function1Types::NSRDSfunc6<Type>::NSRDSfunc6(const NSRDSfunc6& poly)
:
    Function1<Type>(poly),
    Tc_(poly.Tc_),
    coeffs_(poly.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc6<Type>::~NSRDSfunc6()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::NSRDSfunc6<Type>::exponent
(
    const scalar x
) const
{
    // Return the value of the exponent of the exponential function
    // a_*pow(1 - Tr, ((e_*Tr + d_)*Tr + c_)*Tr + b_) where
    // scalar Tr = x/Tc_
    // i.e. return the value of ((e_*Tr + d_)*Tr + c_)*Tr + b_)

    scalar Tr = x/mag(Tc_);

    return ((coeffs_[4]*Tr + coeffs_[3])*Tr + coeffs_[2])*Tr + coeffs_[1];
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc6<Type>::exponentDerivative
(
    const scalar x
) const
{
    // Return the value of the derivative of the exponent of the exponential
    // function
    // a_*pow(1 - Tr, ((e_*Tr + d_)*Tr + c_)*Tr + b_) where
    // scalar Tr = x/Tc_
    // i.e. return the value of (3*e_*x/(Tc_^3) + 2*d_/(Tc_^2))*x + c_/Tc_

    return
        (
            3*coeffs_[4]*x/pow(mag(Tc_),3) + 2*coeffs_[3]/pow(mag(Tc_),2)
        )*x + coeffs_[2]/mag(Tc_);
}


template<class Type>
void Foam::Function1Types::NSRDSfunc6<Type>::convertTimeBase(const Time& t)
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc6<Type>::value(const scalar x) const
{
    // Evaluate the function
    // a_*pow(1 - Tr, ((e_*Tr + d_)*Tr + c_)*Tr + b_) where
    // scalar Tr = x/Tc_
    // and return the result

    scalar Tr = x/mag(Tc_);

    return
        cmptMultiply
        (
            coeffs_[0],
            cmptPow
            (
                pTraits<Type>::one*(1 - Tr),
                ((coeffs_[4]*Tr + coeffs_[3])*Tr + coeffs_[2])*Tr + coeffs_[1]
            )
        );
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc6<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc6<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc6<Type>::derivative(const scalar x) const
{
    // Evaluate the vectorised derivative of the function
    // a_*pow(f(x), g(x)) with
    // f(x) = 1 - Tr , and
    // g(x) = exponent(x) = ((e_*Tr + d_)*Tr + c_)*Tr + b_), and where
    // scalar Tr = x/Tc_ .
    //
    // The derivative of pow(f(x), g(x)) is given by logarithmic differentiation
    // h(x) = pow(f(x), g(x))
    // log(h(x)) = g(x)log(f(x))
    // h'(x)/h(x) = g'(x)*log(f(x)) + g(x)*f'(x)/f(x)
    // h'(x) = h(x)(g'(x)*log(f(x)) + g(x)*f'(x)/f(x)
    // h'(x) = value(x)(exponentDerivative(x)*log(f(x)) + exponent(x)*f'(x)/f(x)
    //
    // Return the value of h'(x).

    scalar Tr = x/mag(Tc_);

    return
        cmptMultiply
        (
            coeffs_[0],
            cmptMultiply
            (
                value(x),
                exponentDerivative(x)*log(1-Tr)
              + exponent(x)/(-mag(Tc_)*log(1-Tr))
            )
        );
}


template<class Type>
void Foam::Function1Types::NSRDSfunc6<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    // Switch to ASCII for writing
    const bool isBinary(os.format() == IOstream::BINARY);
    isBinary ? os.format(IOstream::ASCII) : os.format(IOstream::ASCII);
    os  << nl << indent << Tc_ << token::SPACE << coeffs_
        << token::END_STATEMENT << nl;
    isBinary ? os.format(IOstream::BINARY) : os.format(IOstream::ASCII);
}


// ************************************************************************* //

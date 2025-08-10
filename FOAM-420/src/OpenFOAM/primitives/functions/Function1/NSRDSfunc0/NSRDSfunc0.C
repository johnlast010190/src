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

#include "primitives/functions/Function1/NSRDSfunc0/NSRDSfunc0.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc0<Type>::NSRDSfunc0
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
Foam::Function1Types::NSRDSfunc0<Type>::NSRDSfunc0
(
    const word& entryName,
    const FixedList<Type, 6>& coeffs
)
:
    Function1<Type>(entryName),
    coeffs_(coeffs)
{}


template<class Type>
Foam::Function1Types::NSRDSfunc0<Type>::NSRDSfunc0(const NSRDSfunc0& poly)
:
    Function1<Type>(poly),
    coeffs_(poly.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc0<Type>::~NSRDSfunc0()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::NSRDSfunc0<Type>::convertTimeBase(const Time& t)
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc0<Type>::value(const scalar x) const
{
    // Evaluate the function
    // ((((f_*x + e_)*x + d_)*x + c_)*x + b_)*x + a_
    // = f*x^5 + e*x^4 + d*x^3 + c*x^2 + b*x + a
    // and return the result

    return
        (
            (
                (
                    (
                        coeffs_[5]*x + coeffs_[4]
                    )*x + coeffs_[3]
                )*x + coeffs_[2]
            )*x + coeffs_[1]
        )*x + coeffs_[0];
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc0<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    // Evaluate the integral of the function
    // ((((f_*x + e_)*x + d_)*x + c_)*x + b_)*x + a_
    // = f*x^5 + e*x^4 + d*x^3 + c*x^2 + b*x + a
    //
    // given by
    // (((((f_*x/6 + e_/5)*x + d_/4)*x + c_/3)*x + b_/2)*x + a_)*x
    // = (f*x^6)/6 + (e*x^5)/5 + (d*x^4)/4 + (c*x^3)/3 + (b*x^2)/2 + ax
    //
    // Return the value of the definite integral between x1 and x2

    Type integralOfX2 =
        (
            (
                (
                    (
                        (
                            coeffs_[5]/6*x2 + coeffs_[4]/5
                        )*x2 + coeffs_[3]/4
                    )*x2 + coeffs_[2]/3
                )*x2 + coeffs_[1]/2
            )*x2 + coeffs_[0]
        )*x2;

    Type integralOfX1 =
        (
            (
                (
                    (
                        (
                            coeffs_[5]/6*x1 + coeffs_[4]/5
                        )*x1 + coeffs_[3]/4
                    )*x1 + coeffs_[2]/3
                )*x1 + coeffs_[1]/2
            )*x1 + coeffs_[0]
        )*x1;

    return integralOfX2 - integralOfX1;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc0<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    // Evaluate the integral of the function
    // ((((f_*x + e_)*x + d_)*x + c_)*x + b_) + a_/x
    // = f*x^4 + e*x^3 + d*x^2 + c*x + b + a/x
    //
    // given by
    // ((((f_*x/5 + e_/4)*x + d_/3)*x + c_/2)*x + b_)*x + a_*log(x)
    // = (f*x^5)/5 + (e*x^4)/4 + (d*x^3)/3 + (c*x^2)/2 + b*x + a*log(x)
    //
    // Return the value of the definite integral between x1 and x2

    Type integralOfX2 =
        (
            (
                (
                    (
                        coeffs_[5]/5*x2 + coeffs_[4]/4
                    )*x2 + coeffs_[3]/3
                )*x2 + coeffs_[2]/2
            )*x2 + coeffs_[1]
        )*x2 + (coeffs_[0]*log(x2));

    Type integralOfX1 =
        (
            (
                (
                    (
                        coeffs_[5]/5*x1 + coeffs_[4]/4
                    )*x1 + coeffs_[3]/3
                )*x1 + coeffs_[2]/2
            )*x1 + coeffs_[1]
        )*x1 + (coeffs_[0]*log(x1));

    return integralOfX2 - integralOfX1;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc0<Type>::derivative(const scalar x) const
{
    // Evaluate the derivative of the function
    // ((((f_*x + e_)*x + d_)*x + c_)*x + b_)*x + a_
    // = f*x^5 + e*x^4 + d*x^3 + c*x^2 + b*x + a
    //
    // given by
    // ((((5*f_*x + 4*e_)*x + 3*d_)*x + 2*c_)*x + b_)
    // = 5*f*x^4 + 4*e*x^3 + 3*d*x^2 + 2*c*x + b
    // and return the result

    return
        (
            (
                (
                    5*coeffs_[5]*x + 4*coeffs_[4]
                )*x + 3*coeffs_[3]
            )*x + 2*coeffs_[2]
        )*x + coeffs_[1];
}


template<class Type>
void Foam::Function1Types::NSRDSfunc0<Type>::writeData(Ostream& os) const
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

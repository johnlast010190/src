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

#include "primitives/functions/Function1/NSRDSfunc7/NSRDSfunc7.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc7<Type>::NSRDSfunc7
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
    stream >> coeffs_;
}


template<class Type>
Foam::Function1Types::NSRDSfunc7<Type>::NSRDSfunc7
(
    const word& entryName,
    const FixedList<Type, 5>& coeffs
)
:
    Function1<Type>(entryName),
    coeffs_(coeffs)
{}


template<class Type>
Foam::Function1Types::NSRDSfunc7<Type>::NSRDSfunc7(const NSRDSfunc7& poly)
:
    Function1<Type>(poly),
    coeffs_(poly.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc7<Type>::~NSRDSfunc7()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::NSRDSfunc7<Type>::sinhPart
(
    const scalar x
) const
{
    // Return the value of (c_/x)/sinh(c_/x)

    return cmptDivide(coeffs_[2]/x, cmptSinh(coeffs_[2]/x));
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc7<Type>::coshPart
(
    const scalar x
) const
{
    // Return the value of (e_/x)/cosh(e_/x)

    return cmptDivide(coeffs_[4]/x, cmptCosh(coeffs_[4]/x));
}


template<class Type>
void Foam::Function1Types::NSRDSfunc7<Type>::convertTimeBase(const Time& t)
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc7<Type>::value(const scalar x) const
{
    // Evaluate the function
    // a_ + b_*sqr((c_/x)/sinh(c_/x)) + d_*sqr((e_/x)/cosh(e_/x))
    // and return the result

    return
        coeffs_[0]
      + cmptMultiply(coeffs_[1], cmptSqr(sinhPart(x)))
      + cmptMultiply(coeffs_[3], cmptSqr(coshPart(x)));
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc7<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc7<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc7<Type>::derivative(const scalar x) const
{
    // Evaluate the vectorised derivative of the function
    // a_ + b_*sqr((c_/x)/sinh(c_/x)) + d_*sqr((e_/x)/cosh(e_/x))
    // and return the result

    Type u = coeffs_[2]/x;
    Type v = coeffs_[4]/x;

    return
        cmptMultiply
        (
            2*coeffs_[1],
            cmptMultiply
            (
                sinhPart(x),
                cmptMultiply
                (
                    u,
                    cmptDivide
                    (
                        cmptMultiply(u, cmptCosh(u)) - cmptSinh(u),
                        x*cmptSqr(cmptSinh(u))
                    )
                )
            )
        )
      + cmptMultiply
        (
            2*coeffs_[3],
            cmptMultiply
            (
                coshPart(x),
                cmptMultiply
                (
                    v,
                    cmptDivide
                    (
                        cmptMultiply(v, cmptSinh(v)) - cmptCosh(v),
                        x*cmptSqr(cmptCosh(v))
                    )
                )
            )
        );
}


template<class Type>
void Foam::Function1Types::NSRDSfunc7<Type>::writeData(Ostream& os) const
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

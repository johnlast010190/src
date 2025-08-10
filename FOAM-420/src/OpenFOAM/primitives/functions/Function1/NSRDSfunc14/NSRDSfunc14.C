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

#include "primitives/functions/Function1/NSRDSfunc14/NSRDSfunc14.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "db/IOstreams/Tstreams/ITstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc14<Type>::NSRDSfunc14
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
Foam::Function1Types::NSRDSfunc14<Type>::NSRDSfunc14
(
    const word& entryName,
    const scalar& Tc,
    const FixedList<Type, 4>& coeffs
)
:
    Function1<Type>(entryName),
    Tc_(Tc),
    coeffs_(coeffs)
{}


template<class Type>
Foam::Function1Types::NSRDSfunc14<Type>::NSRDSfunc14(const NSRDSfunc14& poly)
:
    Function1<Type>(poly),
    Tc_(poly.Tc_),
    coeffs_(poly.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::NSRDSfunc14<Type>::~NSRDSfunc14()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::NSRDSfunc14<Type>::convertTimeBase(const Time& t)
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc14<Type>::value(const scalar x) const
{
    // Evaluate the function
    //
    //    a_*a_/(t + ROOTVSMALL) + b_ - t
    //  *(
    //      2.0*a_*c_
    //      + t*(a_*d_ + t*(c_*c_/3.0 + t*(0.5*c_*d_ + 0.2*d_*d_*t)))
    //  )
    //
    // where
    // scalar t = 1.0 - Tdash/Tc_ , and
    // Tdash = min(T, Tc_ - ROOTVSMALL)
    //
    // and return the result

    scalar Tdash = min(x, Tc_ - ROOTVSMALL);

    scalar t = 1.0 - Tdash/Tc_;

    return
        cmptMultiply(coeffs_[0], coeffs_[0])/(t + ROOTVSMALL) + coeffs_[1]
      - (
            2.0*cmptMultiply(coeffs_[0], coeffs_[2])
          + (
                cmptMultiply(coeffs_[0], coeffs_[3])
              + (
                    cmptMultiply(coeffs_[2], coeffs_[2])/3.0
                  + (
                        0.5*cmptMultiply(coeffs_[2], coeffs_[3])
                      + 0.2*cmptMultiply(coeffs_[3], coeffs_[3])*t
                    )*t
                )*t
            )*t
        )*t;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc14<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc14<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1Types::NSRDSfunc14<Type>::derivative(const scalar x) const
{
    // Evaluate the vectorised derivative of the function
    //
    //    a_*a_/(t + ROOTVSMALL) + b_ - t
    //  *(
    //      2.0*a_*c_
    //      + t*(a_*d_ + t*(c_*c_/3.0 + t*(0.5*c_*d_ + 0.2*d_*d_*t)))
    //  )
    //
    // where
    // scalar t = 1.0 - Tdash/Tc_ , and
    // Tdash = min(T, Tc_ - ROOTVSMALL)
    //
    // given by
    //
    //   (a_*a_/((t + ROOTVSMALL)*(t + ROOTVSMALL)) + 2.0*a_*c_ + t
    //  *(
    //      2.0*a_*d_
    //      + t*(c_*c_ + t*(2.0*c_*d_ + t*(d_*d_)))
    //  ))/Tc_
    //
    // and return the result

    scalar Tdash = min(x, Tc_ - ROOTVSMALL);

    scalar t = 1.0 - Tdash/Tc_;

    if (Tdash == x)
    {
        return
            cmptMultiply(coeffs_[0], coeffs_[0])
           /((t + ROOTVSMALL)*(t + ROOTVSMALL))
          + 2*cmptMultiply(coeffs_[0], coeffs_[2])
          + (
                2.0*cmptMultiply(coeffs_[0], coeffs_[3])
              + (
                    cmptMultiply(coeffs_[2], coeffs_[2])
                  + (
                        cmptMultiply(coeffs_[2], coeffs_[3])*2
                      + (
                            cmptMultiply(coeffs_[3], coeffs_[3])
                        )*t
                    )*t
                )*t
            )*t/Tc_;
    }
    else
    {
        return pTraits<Type>::zero;
    }
}


template<class Type>
void Foam::Function1Types::NSRDSfunc14<Type>::writeData(Ostream& os) const
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

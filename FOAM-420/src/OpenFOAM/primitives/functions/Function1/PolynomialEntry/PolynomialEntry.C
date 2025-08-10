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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/PolynomialEntry/PolynomialEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Polynomial<Type>::Polynomial
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName),
    coeffs_(),
    limits_(-GREAT, GREAT)
{
    Istream& is(dict.lookup(entryName));
    word entryType(is);

    is  >> coeffs_;

    if (is.good())
    {
        token t;
        is.read(t);
        is.putBack(t);
        if (t == token::BEGIN_LIST || t.isLabel())
        {
            // Do not read the Pair directly from stream as it does not accept
            // an ascii input in binary mode
            List<scalar> limits(is);
            if (limits.size() != 2)
            {
                FatalIOErrorInFunction(is)
                    << "Expected a pair of values, but got list of length "
                    << limits.size() << exit(FatalIOError);
            }
            limits_.first() = limits[0];
            limits_.second() = limits[1];
            if (limits_.first() > limits_.second())
            {
                FatalErrorInFunction
                    << "Limits for polynomial entry " << this->name_
                    << " are inverted" << nl << exit(FatalError);
            }
        }
    }

    if (!coeffs_.size())
    {
        FatalErrorInFunction
            << "Polynomial coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


template<class Type>
Foam::Function1Types::Polynomial<Type>::Polynomial
(
    const word& entryName,
    const List<Tuple2<Type, Type>>& coeffs,
    const Pair<scalar>& limits
)
:
    Function1<Type>(entryName),
    coeffs_(coeffs),
    limits_(limits)
{
    if (!coeffs_.size())
    {
        FatalErrorInFunction
            << "Polynomial coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
    if (limits_.first() > limits_.second())
    {
        FatalErrorInFunction
            << "Limits for polynomial entry " << this->name_
            << " are inverted" << nl << exit(FatalError);
    }
}


template<class Type>
Foam::Function1Types::Polynomial<Type>::Polynomial(const Polynomial& poly)
:
    Function1<Type>(poly),
    coeffs_(poly.coeffs_),
    limits_(poly.limits_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Polynomial<Type>::~Polynomial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Polynomial<Type>::convertTimeBase(const Time& t)
{
    forAll(coeffs_, i)
    {
        Type value = coeffs_[i].first();
        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            setComponent(coeffs_[i].first(), cmpt) =
                t.userTimeToTime(component(value, cmpt));
        }
    }
}


template<class Type>
Type Foam::Function1Types::Polynomial<Type>::value(const scalar x) const
{
    const scalar xlim = max(limits_.first(), min(limits_.second(), x));

    Type y(Zero);
    forAll(coeffs_, i)
    {
        y += cmptMultiply
        (
            coeffs_[i].first(),
            cmptPow(pTraits<Type>::one*xlim, coeffs_[i].second())
        );
    }

    return y;
}


template<class Type>
Type Foam::Function1Types::Polynomial<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    Type intx(Zero);

    const scalar x1Lim = max(limits_.first(), min(limits_.second(), x1));
    const scalar x2Lim = max(limits_.first(), min(limits_.second(), x2));

    forAll(coeffs_, i)
    {
        // First integrate the portion between the limits
        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            const scalar coeff = component(coeffs_[i].first(), cmpt);
            const scalar exponent = component(coeffs_[i].second(), cmpt);
            if (mag(exponent + scalar(1)) < ROOTVSMALL)
            {
                setComponent(intx, cmpt) += coeff*(log(x2Lim) - log(x1Lim));
            }
            else
            {
                setComponent(intx, cmpt) +=
                    coeff / (exponent + scalar(1))
                   *(
                        pow(x2Lim, exponent + scalar(1))
                      - pow(x1Lim, exponent + scalar(1))
                    );
            }
        }
    }

    // Add integral of any portion outside the limits
    if (x1 < x2)
    {
        intx += (x1Lim-x1)*value(limits_.first());
        intx += (x2-x2Lim)*value(limits_.second());
    }
    else
    {
        intx += (x1Lim-x1)*value(limits_.second());
        intx += (x2-x2Lim)*value(limits_.first());
    }

    return intx;
}


template<class Type>
Type Foam::Function1Types::Polynomial<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    Type intx(Zero);

    const scalar x1Lim = max(limits_.first(), min(limits_.second(), x1));
    const scalar x2Lim = max(limits_.first(), min(limits_.second(), x2));

    forAll(coeffs_, i)
    {
        // First integrate the portion between the limits
        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            const scalar coeff = component(coeffs_[i].first(), cmpt);
            const scalar exponent = component(coeffs_[i].second(), cmpt);
            if (mag(exponent) < ROOTVSMALL)
            {
                setComponent(intx, cmpt) += coeff*(log(x2Lim) - log(x1Lim));
            }
            else
            {
                setComponent(intx, cmpt) +=
                    coeff / exponent
                   *(
                        pow(x2Lim, exponent)
                      - pow(x1Lim, exponent)
                    );
            }
        }
    }

    // Add integral of any portion outside the limits
    if (x1 < x2)
    {
        intx += (log(x1Lim)-log(x1))*value(limits_.first());
        intx += (log(x2)-log(x2Lim))*value(limits_.second());
    }
    else
    {
        intx += (log(x1Lim)-log(x1))*value(limits_.second());
        intx += (log(x2)-log(x2Lim))*value(limits_.first());
    }

    return intx;
}


template<class Type>
Type Foam::Function1Types::Polynomial<Type>::derivative(const scalar x) const
{
    Type deriv(Zero);

    // Outside of the limites the function is constant hence derivitive is zero
    if (limits_.first() <= x  || limits_.second() >= x)
    {
        return deriv;
    }

    forAll(coeffs_, i)
    {
        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            const scalar coeff = component(coeffs_[i].first(), cmpt);
            const scalar exponent = component(coeffs_[i].second(), cmpt);
            if (mag(exponent) > ROOTVSMALL)
            {
                if (mag(exponent - 1.0) > ROOTVSMALL)
                {
                    setComponent(deriv, cmpt) +=
                        coeff*exponent*pow(x, exponent - 1.0);
                }
                else
                {
                    setComponent(deriv, cmpt) += coeff*exponent;
                }
            }
        }
    }

    return deriv;
}


template<class Type>
void Foam::Function1Types::Polynomial<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    // Switch to ASCII for writing the polynomial and limits
    const bool isBinary(os.format() == IOstream::BINARY);
    isBinary ? os.format(IOstream::ASCII) : os.format(IOstream::ASCII);
    os  << nl << indent << coeffs_
        << nl << indent << limits_
        << token::END_STATEMENT << nl;
    isBinary ? os.format(IOstream::BINARY) : os.format(IOstream::ASCII);
}


// ************************************************************************* //

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
    (c) 2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "Function2.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2<Type>::Function2(const word& name)
:
    name_(name)
{}


template<class Type>
Foam::Function2<Type>::Function2(const Function2<Type>& de)
:
    tmp<Function2<Type>>::refCount(),
    name_(de.name_)
{}


template<class Type, class Function2Type>
Foam::FieldFunction2<Type, Function2Type>::FieldFunction2
(
    const word& name
)
:
    Function2<Type>(name)
{}


template<class Type, class Function2Type>
Foam::tmp<Foam::Function2<Type>>
Foam::FieldFunction2<Type, Function2Type>::clone() const
{
    return tmp<Function2<Type>>
    (
        new Function2Type(refCast<const Function2Type>(*this))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2<Type>::~Function2()
{}


template<class Type, class Function2Type>
Foam::FieldFunction2<Type, Function2Type>::~FieldFunction2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::Function2<Type>::name() const
{
    return name_;
}


template<class Type, class Function2Type>
Foam::tmp<Foam::Field<Type>> Foam::FieldFunction2<Type, Function2Type>::value
(
    const scalarField& x,
    const scalarField& y
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x, i)
    {
        fld[i] = refCast<const Function2Type>(*this).value(x[i], y[i]);
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function2<Type>::operator=(const Function2<Type>& f)
{
    if (this == &f)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::Function2<Type>::writeData(Ostream& os) const
{
    os.writeKeyword(name_) << type();
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Function2<Type>& f2
)
{
    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const Function2<Type>&)"
    );

    f2.write(os);

    return os;
}


// ************************************************************************* //

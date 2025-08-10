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

#include "primitives/functions/Function1/Function1/Function1.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1<Type>::Function1(const word& entryName)
:
    name_(entryName)
{}


template<class Type>
Foam::Function1<Type>::Function1(const Function1<Type>& de)
:
    tmp<Function1<Type>>::refCount(),
    name_(de.name_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1<Type>::~Function1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::Function1<Type>::name() const
{
    return name_;
}


template<class Type>
void Foam::Function1<Type>::convertTimeBase(const Time&)
{}


template<class Type>
void Foam::Function1<Type>::shiftTimeBase(const Time&)
{}


template<class Type>
Type Foam::Function1<Type>::value(const scalar x) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1<Type>::value(const vector& xyz) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1<Type>::derivative(const scalar x) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1<Type>::integrate(const scalar x1, const scalar x2) const
{
    NotImplemented;
}


template<class Type>
Type Foam::Function1<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1<Type>::value
(
    const scalarField& x
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x, i)
    {
        fld[i] = this->value(x[i]);
    }
    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1<Type>::value
(
    const Field<vector>& x
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x, i)
    {
        fld[i] = this->value(x[i]);
    }
    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1<Type>::integrate
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x1.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x1, i)
    {
        fld[i] = this->integrate(x1[i], x2[i]);
    }
    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1<Type>::integrateYoverX
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x1.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x1, i)
    {
        fld[i] = this->integrateYoverX(x1[i], x2[i]);
    }
    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1<Type>::derivative
(
    const scalarField& x1
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x1.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x1, i)
    {
        fld[i] = this->derivative(x1[i]);
    }
    return tfld;
}


template<class Type>
void Foam::Function1<Type>::writeData(Ostream& os) const
{
    os.writeKeyword(name_) << type();
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Function1<Type>& f1
)
{
    os.check(FUNCTION_NAME);

    f1.writeData(os);

    return os;
}


// ************************************************************************* //

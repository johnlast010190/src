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
    (c) 2016 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/IOobjects/GlobalIOField/GlobalIOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::GlobalIOField<Type>::GlobalIOField(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOField<Type>>();

    readHeaderOk(IOstream::BINARY, typeName);
}


template<class Type>
Foam::GlobalIOField<Type>::GlobalIOField(const IOobject& io, const label size)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOField<Type>>();

    if (!readHeaderOk(IOstream::BINARY, typeName))
    {
        Field<Type>::setSize(size);
    }
}


template<class Type>
Foam::GlobalIOField<Type>::GlobalIOField
(
    const IOobject& io,
    const Field<Type>& f
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOField<Type>>();

    if (!readHeaderOk(IOstream::BINARY, typeName))
    {
        Field<Type>::operator=(f);
    }
}


template<class Type>
Foam::GlobalIOField<Type>::GlobalIOField
(
    const IOobject& io,
    const Xfer<Field<Type>>& f
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOField<Type>>();

    Field<Type>::transfer(f());

    readHeaderOk(IOstream::BINARY, typeName);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::GlobalIOField<Type>::~GlobalIOField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::GlobalIOField<Type>::readData(Istream& is)
{
    is >> *this;
    return is.good();
}


template<class Type>
bool Foam::GlobalIOField<Type>::writeData(Ostream& os) const
{
    return (os << static_cast<const Field<Type>&>(*this)).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::GlobalIOField<Type>::operator=(const GlobalIOField<Type>& rhs)
{
    Field<Type>::operator=(rhs);
}


template<class Type>
void Foam::GlobalIOField<Type>::operator=(const Field<Type>& rhs)
{
    Field<Type>::operator=(rhs);
}


// ************************************************************************* //

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
    (c) Creative Fields, Ltd.
    (c) 2018-2021 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

\*---------------------------------------------------------------------------*/

#include "utilities/containers/DynList/DynList.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::DynList<T, SizeMin>::DynList(Istream& is)
:
    UList<T>(),
    heapList_(),
    capacity_(0)
{
    is >> *this;
}


template<class T, int SizeMin>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Foam::DynList<T, SizeMin>& list
)
{
    os << static_cast<const UList<T>&>(list);
    return os;
}


template<class T, int SizeMin>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    Foam::DynList<T, SizeMin>& list
)
{
    list.clearStorage();

    List<T> input(is);

    const label newLen = input.size();

    list.heapList_ = std::move(input);

    list.setSize(newLen);

    return is;
}


// ************************************************************************* //

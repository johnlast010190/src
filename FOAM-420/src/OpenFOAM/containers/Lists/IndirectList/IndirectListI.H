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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


inline Foam::IndirectListAddressing::IndirectListAddressing
(
    const labelUList& addr
)
:
    addressing_(addr)
{}


inline Foam::IndirectListAddressing::IndirectListAddressing
(
    const Xfer<List<label>>& addr
)
:
    addressing_(addr)
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& completeList,
    const labelUList& addr
)
:
    IndirectListAddressing(addr),
    UIndirectList<T>
    (
        completeList,
        IndirectListAddressing::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& completeList,
    const Xfer<List<label>>& addr
)
:
    IndirectListAddressing(addr),
    UIndirectList<T>
    (
        completeList,
        IndirectListAddressing::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const IndirectList<T>& lst
)
:
    IndirectListAddressing(lst.addressing()),
    UIndirectList<T>
    (
        lst.completeList(),
        IndirectListAddressing::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UIndirectList<T>& lst
)
:
    IndirectListAddressing(lst.addressing()),
    UIndirectList<T>
    (
        lst.completeList(),
        IndirectListAddressing::addressing()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::List<Foam::label>&
Foam::IndirectListAddressing::addressing() const
{
    return addressing_;
}


inline void Foam::IndirectListAddressing::resetAddressing
(
    const labelUList& addr
)
{
    addressing_ = addr;
}


inline void Foam::IndirectListAddressing::resetAddressing
(
    const Xfer<List<label>>& addr
)
{
    addressing_.transfer(addr());
}


template<class T>
inline void Foam::IndirectList<T>::shallowCopy(const IndirectList<T>& a)
{
    UIndirectList<T>::operator=(a);
}


// ************************************************************************* //

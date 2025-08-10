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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const label subSize
)
:
    UList<T>(list.v_, subSize)
{
    #ifdef FULLDEBUG
    list.checkSize(subSize);
    #endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const label subSize,
    const label startIndex
)
:
    UList<T>(&(list.v_[startIndex]), subSize)
{
    #ifdef FULLDEBUG
    // Artificially allow the start of a zero-sized subList to be
    // one past the end of the original list.
    if (subSize)
    {
        list.checkStart(startIndex);
        list.checkSize(startIndex + subSize);
    }
    else
    {
        // Start index needs to fall between 0 and size.  One position
        // behind the last element is allowed
        list.checkSize(startIndex);
    }
    #endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const labelRange& range
)
:
    SubList<T>(list.validateRange(range), list)
{}


template<class T>
inline Foam::SubList<T>::SubList
(
    const labelRange& range,
    const UList<T>& list
)
:
    SubList<T>(list, range.size(), range.start())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::SubList<T>& Foam::SubList<T>::null()
{
    return NullObjectRef<SubList<T>>();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<class T>
inline void Foam::SubList<T>::operator=(const SubList<T>& list)
{
    UList<T>::deepCopy(list);
}


template<class T>
inline void Foam::SubList<T>::operator=(const UList<T>& list)
{
    UList<T>::deepCopy(list);
}


template<class T>
inline void Foam::SubList<T>::operator=(const T& t)
{
    UList<T>::operator=(t);
}


// ************************************************************************* //

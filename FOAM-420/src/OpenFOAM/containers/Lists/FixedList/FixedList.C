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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "containers/Lists/FixedList/FixedList.H"
#include "containers/Lists/List/ListLoopM.H"

// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T, unsigned Size>
void Foam::FixedList<T, Size>::swap(FixedList<T, Size>& a)
{
    List_ACCESS(T, (*this), vp);
    List_ACCESS(T, a, ap);
    T tmp;
    List_FOR_ALL((*this), i)
        tmp = List_CELEM((*this), vp, i);
        List_ELEM((*this), vp, i) = List_CELEM(a, ap, i);
        List_ELEM(a, ap, i) = tmp;
    List_END_FOR_ALL
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned N>
Foam::label Foam::FixedList<T, N>::find(const T& val, label pos) const
{
    if (pos >= 0)
    {
        List_CONST_ACCESS(T, *this, list);

        while (pos < label(N))
        {
            if (list[pos] == val)
            {
                return pos;
            }

            ++pos;
        }
    }

    return -1;
}


template<class T, unsigned N>
Foam::label Foam::FixedList<T, N>::rfind(const T& val, label pos) const
{
    // pos == -1 has same meaning as std::string::npos - search from end
    if (pos < 0 || pos >= label(N))
    {
        pos = label(N)-1;
    }

    List_CONST_ACCESS(T, *this, list);

    while (pos >= 0)
    {
        if (list[pos] == val)
        {
            return pos;
        }

        --pos;
    }

    return -1;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator==(const FixedList<T, Size>& a) const
{
    bool equal = true;

    List_CONST_ACCESS(T, (*this), vp);
    List_CONST_ACCESS(T, (a), ap);

    List_FOR_ALL((*this), i)
        equal = (List_ELEM((*this), vp, i) == List_ELEM((a), ap, i));
        if (!equal) break;
    List_END_FOR_ALL

    return equal;
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator!=(const FixedList<T, Size>& a) const
{
    return !operator==(a);
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator<(const FixedList<T, Size>& a) const
{
    const T* const __restrict__ ptr1 = this->begin();
    const T* const __restrict__ ptr2 = a.begin();

    for (unsigned i=0; i<Size; ++i)
    {
        if (ptr1[i] < ptr2[i])
        {
            return true;
        }
        else if (ptr1[i] > ptr2[i])
        {
            return false;
        }
    }

    // Contents look to be identical.
    // The sizes are identical by definition (template parameter)
    return false;
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator>(const FixedList<T, Size>& a) const
{
    return a.operator<(*this);
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator<=(const FixedList<T, Size>& a) const
{
    return !operator>(a);
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator>=(const FixedList<T, Size>& a) const
{
    return !operator<(a);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "containers/Lists/FixedList/FixedListIO.C"

// ************************************************************************* //

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

#include "containers/LinkedLists/accessTypes/UILList/UILList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::UILList<LListBase, T>::UILList(const UILList<LListBase, T>& lst)
{
    for (const_iterator iter = lst.begin(); iter != lst.end(); ++iter)
    {
        this->append(&iter());
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::UILList<LListBase, T>::operator=(const UILList<LListBase, T>& rhs)
{
    LListBase::clear();

    for (const_iterator iter = rhs.begin(); iter != rhs.end(); ++iter)
    {
        this->append(&iter());
    }
}


template<class LListBase, class T>
bool Foam::UILList<LListBase, T>::operator==
(
    const UILList<LListBase, T>& rhs
) const
{
    bool equal = (this->size() == rhs.size());
    if (!equal)
    {
        return false;
    }

    const_iterator iter1 = this->begin();
    const_iterator iter2 = rhs.begin();

    for (; iter1 != this->end(); ++iter1, ++iter2)
    {
        equal = (iter1() == iter2());
        if (!equal) break;
    }

    return equal;
}


template<class LListBase, class T>
bool Foam::UILList<LListBase, T>::operator!=
(
    const UILList<LListBase, T>& rhs
) const
{
    return !operator==(rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "containers/LinkedLists/accessTypes/UILList/UILListIO.C"


// ************************************************************************* //

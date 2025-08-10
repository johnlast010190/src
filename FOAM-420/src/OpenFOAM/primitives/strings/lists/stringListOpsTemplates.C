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
    (c) 2017 OpenCFD Ltd.
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class UnaryMatchPredicate, class StringType>
Foam::labelList Foam::findMatchingStrings
(
    const UnaryMatchPredicate& matcher,
    const UList<StringType>& lst,
    const bool invert
)
{
    labelList indices(lst.size());

    label count = 0;
    forAll(lst, elemi)
    {
#if defined(WIN64) || defined(WIN32)
        if (matcher.match(lst[elemi]) ? !invert : invert)
#else
        if (matcher(lst[elemi]) ? !invert : invert)
#endif
        {
            indices[count++] = elemi;
        }
    }
    indices.setSize(count);

    return indices;
}


template<class UnaryMatchPredicate, class StringListType>
StringListType Foam::subsetMatchingStrings
(
    const UnaryMatchPredicate& matcher,
    const StringListType& lst,
    const bool invert
)
{
    // Create as a copy
    StringListType newLst(lst.size());

    // Ensure consistent addressable size (eg, DynamicList)
    newLst.setSize(lst.size());

    label count = 0;
    forAll(lst, elemi)
    {
        if (matcher(lst[elemi]) ? !invert : invert)
        {
            newLst[count++] = lst[elemi];
        }
    }
    newLst.setSize(count);

    return newLst;
}


template<class UnaryMatchPredicate, class StringListType>
void Foam::inplaceSubsetMatchingStrings
(
    const UnaryMatchPredicate& matcher,
    StringListType& lst,
    const bool invert
)
{
    label count = 0;
    forAll(lst, elemi)
    {
        if (matcher(lst[elemi]) ? !invert : invert)
        {
            if (count != elemi)
            {
                lst[count] = lst[elemi];
            }
            ++count;
        }
    }
    lst.setSize(count);
}


// ************************************************************************* //

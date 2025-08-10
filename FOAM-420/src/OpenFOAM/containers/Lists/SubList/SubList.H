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
    (c) 2017-2021 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::SubList

Description
    A List obtained as a section of another List.

    Since the SubList is itself unallocated, no storage is allocated or
    de-allocated during its use.  To achieve this behaviour, SubList is
    derived from UList rather than List.

SourceFiles
    SubListI.H

\*---------------------------------------------------------------------------*/

#ifndef SubList_H
#define SubList_H

#include "containers/Lists/List/List.H"
#include "primitives/ranges/labelRange/labelRange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SubList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class SubList
:
    public UList<T>
{
public:

    // Static Member Functions

        //- Return a null SubList
        inline static const SubList<T>& null();


    // Generated Methods

        //- Default construct, zero-sized and nullptr
        SubList() noexcept = default;

        //- Copy construct, shallow copy
        SubList(const SubList&) noexcept = default;


    // Constructors

        //- Construct from UList and sub-list size
        inline SubList
        (
            const UList<T>& list,
            const label subSize
        );

        //- Construct from UList, sub-list size and start index
        inline SubList
        (
            const UList<T>& list,
            const label subSize,
            const label startIndex
        );

        //- Construct from UList and a (start,size) range.
        //  The range is subsetted with the list size itself to ensure that the
        //  result always addresses a valid section of the list.
        inline SubList
        (
            const UList<T>& list,
            const labelRange& range
        );

        //- Construct from UList and a (start,size) range, but bypassing
        //  run-time range checking.
        inline SubList
        (
            const labelRange& range,
            const UList<T>& list
        );


    // Member operators

        //- Assignment of all entries to the given sub-list
        inline void operator=(const SubList<T>& list);

        //- Assignment of all entries to the given list
        inline void operator=(const UList<T>& list);

        //- Assignment of all entries to the given value
        inline void operator=(const T& t);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "containers/Lists/SubList/SubListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

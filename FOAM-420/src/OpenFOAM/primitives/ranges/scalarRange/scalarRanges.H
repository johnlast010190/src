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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::scalarRanges

Description
    A List of scalarRange.

SourceFiles
    scalarRanges.C

\*---------------------------------------------------------------------------*/

#ifndef scalarRanges_H
#define scalarRanges_H

#include "primitives/ranges/scalarRange/scalarRange.H"
#include "containers/Lists/List/List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class scalarRanges Declaration
\*---------------------------------------------------------------------------*/

class scalarRanges
:
    public List<scalarRange>
{
public:

    // Constructors

        //- Construct null
        scalarRanges();

        //- Construct from Istream.
        //  The list items are comma-delimited.
        scalarRanges(Istream& is);


    // Member Functions

        //- Return true if the given value is within the ranges
        bool selected(const scalar value) const;

        //- Return the set of selected entries in the given list
        //  that are within the ranges
        List<bool> selected(const UList<scalar>& values) const;

        //- Select a list of values that are within the ranges
        List<scalar> select(const List<scalar>& values) const;

        //- Select a list of values that are within the ranges
        void inplaceSelect(List<scalar>& values) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

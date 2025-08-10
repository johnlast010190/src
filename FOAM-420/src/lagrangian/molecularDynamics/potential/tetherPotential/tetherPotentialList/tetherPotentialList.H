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
    Foam::tetherPotentialList

Description

SourceFiles
    tetherPotentialList.C

\*---------------------------------------------------------------------------*/

#ifndef tetherPotentialList_H
#define tetherPotentialList_H

#include "containers/Lists/ListOps/ListOps.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "primitives/strings/word/word.H"
#include "tetherPotential/basic/tetherPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class tetherPotentialList Declaration
\*---------------------------------------------------------------------------*/

class tetherPotentialList
:
    public PtrList<tetherPotential>
{
    // Private data

        List<label> idMap_;


    // Private Member Functions

        inline label tetherPotentialIndex
        (
            const label a
        ) const;

        void readTetherPotentialDict
        (
            const List<word>& siteIdList,
            const dictionary& tetherPotentialDict,
            const List<word>& tetherSiteIdList
        );

        //- Disallow default bitwise assignment
        void operator=(const tetherPotentialList&);

        //- Disallow default bitwise copy construct
        tetherPotentialList(const tetherPotentialList&);


public:

    // Constructors

        tetherPotentialList();

        //- Construct from siteIdList and potental dictionaries
        tetherPotentialList
        (
            const List<word>& siteIdList,
            const dictionary& tetherPotentialDict,
            const List<word>& tetherSiteIdList
        );


    //- Destructor
    ~tetherPotentialList();


    // Member Functions

        void buildPotentials
        (
            const List<word>& siteIdList,
            const dictionary& tetherPotentialDict,
            const List<word>& tetherSiteIdList
        );

        // Access

            inline const List<word>& idMap() const;

            const tetherPotential& tetherPotentialFunction(const label a) const;

            vector force(const label a, const vector rIT) const;

            scalar energy (const label a, const vector rIT) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetherPotential/tetherPotentialList/tetherPotentialListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

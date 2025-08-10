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
    Foam::patchInteractionDataList

Description
    List container for patchInteractionData class

\*---------------------------------------------------------------------------*/

#ifndef patchInteractionDataList_H
#define patchInteractionDataList_H

#include "submodels/Kinematic/PatchInteractionModel/LocalInteraction/patchInteractionData.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/dictionary/dictionary.H"
#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class patchInteractionDataList Declaration
\*---------------------------------------------------------------------------*/

class patchInteractionDataList
:
    public List<patchInteractionData>
{
private:

    // Private data

        //- List of patch IDs for each patch group
        labelListList patchGroupIDs_;


public:

    // Constructor

        //- Construct null
        patchInteractionDataList();

        //- Construct copy
        patchInteractionDataList(const patchInteractionDataList& pidl);

        //- Construct from Istream
        patchInteractionDataList(const polyMesh& mesh, const dictionary& dict);


    // Member functions

        //- Return label of group containing patch id
        //  Returns -1 if patch id is not present
        label applyToPatch(const label id) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

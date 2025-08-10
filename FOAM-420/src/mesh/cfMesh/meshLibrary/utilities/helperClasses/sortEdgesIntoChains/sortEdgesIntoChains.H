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

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    sortEdgesIntoChains

Description
    Creates boundary faces from a chain of open edges

SourceFiles
    sortEdgesIntoChains.C

\*---------------------------------------------------------------------------*/

#ifndef sortEdgesIntoChains_H
#define sortEdgesIntoChains_H

#include "utilities/containers/DynList/DynList.H"
#include "primitives/ints/lists/labelList.H"
#include "meshes/meshShapes/edge/edge.H"
#include "containers/HashTables/Map/Map.H"
#include "primitives/bools/lists/boolList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sortEdgesIntoChains Declaration
\*---------------------------------------------------------------------------*/

class sortEdgesIntoChains
{
    // Members
        const DynList<edge>& bEdges_;

        bool openEdges_;

        Map<label> newNodeLabel_;

        List<DynList<label>> edgesAtPoint_;

        DynList<labelList> createdChains_;

    // Private member functions
        void createNodeLabels();

        bool findPointsBelongingToTheChain
        (
            const label currPos,
            boolList& chainEdges
        ) const;

        void shrinkEdges(const boolList& chainEdges);

        void createChainFromEdges(const boolList& chainEdges);

        void sortEdges();

    public:

        sortEdgesIntoChains
        (
            const DynList<edge>& bEdges
        );

    // Member functions

        //- a list of points which have not yet been resolved
        const DynList<labelList>& sortedChains() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

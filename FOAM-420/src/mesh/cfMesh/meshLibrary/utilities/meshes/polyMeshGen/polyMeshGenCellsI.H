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

Description


\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/meshes/polyMeshGen/polyMeshGenCells.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const cellListPMG& polyMeshGenCells::cells() const
{
    return cells_;
}

inline void polyMeshGenCells::addCellToSubset
(
    const label selID,
    const label cellI
)
{
    std::map<label, meshSubset>::iterator it = cellSubsets_.find(selID);
    if (it == cellSubsets_.end())
        return;

    it->second.addElement(cellI);
}

inline void polyMeshGenCells::removeCellFromSubset
(
    const label selID,
    const label cellI
)
{
    std::map<label, meshSubset>::iterator it = cellSubsets_.find(selID);
    if (it == cellSubsets_.end())
        return;

    it->second.removeElement(cellI);
}

inline void polyMeshGenCells::cellInSubsets
(
    const label cellI,
    DynList<label>& cellSubsets
) const
{
    cellSubsets.clear();

    std::map<label, meshSubset>::const_iterator it;
    for
    (
        it=cellSubsets_.begin();
        it!=cellSubsets_.end();
        ++it
    )
    {
        if (it->second.contains(cellI))
            cellSubsets.append(it->first);
    }
}

inline void polyMeshGenCells::cellSubsetIndices(DynList<label>& indices) const
{
    indices.clear();

    std::map<label, meshSubset>::const_iterator it;
    for
    (
        it=cellSubsets_.begin();
        it!=cellSubsets_.end();
        ++it
    )
        indices.append(it->first);
}

template<class ListType>
inline void polyMeshGenCells::cellsInSubset
(
    const label selID,
    ListType& cellLabels
) const
{
    cellLabels.clear();

    std::map<label, meshSubset>::const_iterator it =
        cellSubsets_.find(selID);
    if (it == cellSubsets_.end())
        return;

    it->second.containedElements(cellLabels);
}

template<class ListType>
inline void polyMeshGenCells::updateCellSubsets
(
    const ListType& newCellLabels
)
{
    for
    (
        std::map<label, meshSubset>::iterator it=cellSubsets_.begin();
        it!=cellSubsets_.end();
        ++it
    )
        it->second.updateSubset(newCellLabels);
}

inline void polyMeshGenCells::updateCellSubsets
(
    const VRWGraph& newCellsForCell
)
{
    for
    (
        std::map<label, meshSubset>::iterator it=cellSubsets_.begin();
        it!=cellSubsets_.end();
        ++it
    )
        it->second.updateSubset(newCellsForCell);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

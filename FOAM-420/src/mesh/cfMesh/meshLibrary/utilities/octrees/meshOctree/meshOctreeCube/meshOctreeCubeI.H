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

#include "utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCube.H"
#include "triSurface/triSurface.H"
#include "utilities/octrees/meshOctree/meshOctreeSlot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline meshOctreeCube::meshOctreeCube()
:
    meshOctreeCubeBasic(),
    activeSlotPtr_(nullptr),
    subCubesPtr_(nullptr),
    cubeLabel_(-1),
    containedElementsLabel_(-1),
    containedEdgesLabel_(-1)
{}

inline meshOctreeCube::meshOctreeCube(const meshOctreeCube& moc)
:
    meshOctreeCubeBasic(moc.coordinates(), moc.cubeType(), moc.procNo()),
    activeSlotPtr_(moc.activeSlotPtr_),
    subCubesPtr_(moc.subCubesPtr_),
    cubeLabel_(moc.cubeLabel_),
    containedElementsLabel_(moc.containedElementsLabel_),
    containedEdgesLabel_(moc.containedEdgesLabel_)
{}

inline const meshOctreeSlot* meshOctreeCube::slotPtr() const
{
    return activeSlotPtr_;
}

inline bool meshOctreeCube::isLeaf() const
{
    if (!subCubesPtr_)
        return true;

    return false;
}

inline meshOctreeCube* meshOctreeCube::subCube(const label scI) const
{
    return subCubesPtr_[scI];
}

inline label meshOctreeCube::cubeLabel() const
{
    return cubeLabel_;
}

inline bool meshOctreeCube::hasContainedElements() const
{
    if (containedElementsLabel_ == -1)
        return false;

    return true;
}

inline label meshOctreeCube::containedElements() const
{
    return containedElementsLabel_;
}

inline bool meshOctreeCube::hasContainedEdges() const
{
    if (containedEdgesLabel_ == -1)
        return false;;

    return true;
}

inline label meshOctreeCube::containedEdges() const
{
    return containedEdgesLabel_;
}

inline void meshOctreeCube::operator=(const meshOctreeCube& moc)
{
    static_cast<meshOctreeCubeCoordinates&>(*this) = moc.coordinates();
    activeSlotPtr_ = moc.activeSlotPtr_;
    subCubesPtr_ = moc.subCubesPtr_;
    cubeLabel_ = moc.cubeLabel_;
    containedElementsLabel_ = moc.containedElementsLabel_;
    containedEdgesLabel_ = moc.containedEdgesLabel_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

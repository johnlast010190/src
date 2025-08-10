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
    And the precice adapter solver object is based on the preCICE-
    adapter for OpenFOAM.

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
    (c) 2017-2023 Gerasimos Chourdakis
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "Displacement.H"

using namespace Foam;

preciceAdapter::FSI::Displacement::Displacement(
    const Foam::fvMesh& mesh,
    const std::string namePointDisplacement,
    const std::string nameCellDisplacement)
: pointDisplacement_(
    namePointDisplacement == "unused"
        ? nullptr
        : const_cast<pointVectorField*>(
            &mesh.lookupObject<pointVectorField>(namePointDisplacement))),
  cellDisplacement_(
      const_cast<volVectorField*>(
          &mesh.lookupObject<volVectorField>(nameCellDisplacement))),
  mesh_(mesh)
{
    dataType_ = vector;
}

// We cannot do this step in the constructor by design of the adapter since the information of the CouplingDataUser is
// defined later. Hence, we call this method after the CouplingDaaUser has been configured
void preciceAdapter::FSI::Displacement::initialize()
{
    // Initialize appropriate objects for each interface patch, namely the volField and the interpolation object
    // this is only necessary for face based FSI
    if (this->locationType_ == LocationType::faceCenters)
    {
        primitivePatches_.clear();
        interpolationObjects_.clear();
        for (unsigned int j = 0; j < patchIDs_.size(); ++j)
        {
            const unsigned int patchID = patchIDs_.at(j);

            primitivePatches_.append
            (
                polyPatchToPrimitivePatch(mesh_.boundaryMesh()[patchID])
            );
            interpolationObjects_.append
            (
                new primitivePatchInterpolation(primitivePatches_[j])
            );
        }
    }
}


void preciceAdapter::FSI::Displacement::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    /* TODO: Implement
     * We need two nested for-loops for each patch,
     * the outer for the locations and the inner for the dimensions.
     * See the preCICE writeBlockVectorData() implementation.
     */

    // Copy the displacement field from OpenFOAM to the buffer

    if (this->locationType_ == LocationType::faceCenters)
    {
        // For every boundary patch of the interface
        for (const label patchID : patchIDs_)
        {
            // Write the displacement to the preCICE buffer
            // For every cell of the patch
            forAll(cellDisplacement_->boundaryField()[patchID], i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    buffer[i * dim + d] =
                        cellDisplacement_->boundaryField()[patchID][i][d];
            }
        }
    }
    else if (this->locationType_ == LocationType::faceNodes)
    {
        DEBUG(adapterInfo(
            "Please be aware of issues with using 'locationType faceNodes' "
            "in parallel. \n"
            "See https://github.com/precice/openfoam-adapter/issues/153.",
            "warning"));

        // For every boundary patch of the interface
        for (const label patchID : patchIDs_)
        {
            // Write the displacement to the preCICE buffer
            // For every cell of the patch
            forAll(pointDisplacement_->boundaryField()[patchID], i)
            {
                const labelList& meshPoints =
                    mesh_.boundaryMesh()[patchID].meshPoints();

                for (unsigned int d = 0; d < dim; ++d)
                    buffer[i * dim + d] =
                        pointDisplacement_->internalField()[meshPoints[i]][d];
            }
        }
    }
}


// return the displacement to use later in the velocity?
void preciceAdapter::FSI::Displacement::read(double* buffer, const unsigned int dim)
{
    for (unsigned int j = 0; j < patchIDs_.size(); j++)
    {
        // Get the ID of the current patch
        const unsigned int patchID = patchIDs_.at(j);

        if (this->locationType_ == LocationType::faceCenters)
        {

            // the boundaryCellDisplacement is a vector and ordered according to the iterator j
            // and not according to the patchID
            // First, copy the buffer data into the center based vectorFields on each interface patch
            forAll(cellDisplacement_->boundaryField()[patchID], i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    cellDisplacement_->boundaryFieldRef()[patchID][i][d] = buffer[i * dim + d];
            }

            if (pointDisplacement_ != nullptr)
            {
                // Get a reference to the displacement on the point patch in order to overwrite it
                vectorField& pointDisplacementFluidPatch(
                    refCast<vectorField>(
                        pointDisplacement_->boundaryFieldRef()[patchID]));

                // Overwrite the node based patch using the interpolation objects and the cell based vector field
                // Afterwards, continue as usual
                pointDisplacementFluidPatch = interpolationObjects_[j].faceToPointInterpolate(cellDisplacement_->boundaryField()[patchID]);
            }
        }
        else if (this->locationType_ == LocationType::faceNodes)
        {

            // Get the displacement on the patch
            fixedValuePointPatchVectorField& pointDisplacementFluidPatch(
                refCast<fixedValuePointPatchVectorField>(
                    pointDisplacement_->boundaryFieldRef()[patchID]));

            // Overwrite the nodes on the interface directly
            forAll(pointDisplacement_->boundaryFieldRef()[patchID], i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    pointDisplacementFluidPatch[i][d] = buffer[i * dim + d];
            }
        }
    }
}

bool preciceAdapter::FSI::Displacement::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters || this->locationType_ == LocationType::faceNodes);
}

std::string preciceAdapter::FSI::Displacement::getDataName() const
{
    return "Displacement";
}

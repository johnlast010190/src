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

Class
    preciceAdapter::FSI::Displacement

SourceFiles
    Displacement.C

\*---------------------------------------------------------------------------*/

#ifndef FSI_DISPLACEMENT_H
#define FSI_DISPLACEMENT_H

#include "CouplingDataUser.H"

#include "cfdTools/general/include/fvCFD.H"
#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchFields.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads displacement
class Displacement : public CouplingDataUser
{

private:
    // Displacement pointVectorField
    Foam::pointVectorField* pointDisplacement_;

    // cellDisplacement field
    Foam::volVectorField* cellDisplacement_;

    const Foam::fvMesh& mesh_;
    // NOTE: this allocation could be avoided in case we directly write to the
    // Needs to be a pointer since the class disables assignment and copy constructors
    PtrList<const primitivePatchInterpolation> interpolationObjects_;

    // List of patches that are used for the interpolation
    PtrList<primitivePatch> primitivePatches_;

public:
    //- Constructor
    Displacement(
        const Foam::fvMesh& mesh,
        const std::string namePointDisplacement,
        const std::string nameCellDisplacement);

    //- Write the displacement values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read the displacement values from the buffer
    void read(double* buffer, const unsigned int dim) final;

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;

    //- We need to initialize the cell-based vector and the interpolation object
    // in case we want to use the faceCenter location for the coupling
    void initialize() final;
};

}
}

#endif

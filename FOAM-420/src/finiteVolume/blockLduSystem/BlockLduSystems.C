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
    (c) 2022 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "blockLduSystem/BlockLduSystem.H"
#include "blockLduSystem/BlockLduSystem.C"
#include "blockLduSystem/BlockLduSystemOperations.C"

// * * * * * * * * * * * * * * * Specialisations * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockLduSystem<vector, scalar>::addContinuityCoupledBC
(
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const volScalarField& p,
    const surfaceTensorField& rDf
)
{
    const typename GeometricField<vector, fvPatchField, volMesh>::
        Boundary& bFields = U.boundaryField();
    forAll(bFields, patchi)
    {
        bFields[patchi].addContinuityCoupledBC(*this, p, rDf);
    }
}

}

// * * * * * * * * * * * * * * * Instantiations  * * * * * * * * * * * * * * //

namespace Foam
{
template class BlockLduSystem<scalar,scalar>;
template class BlockLduSystem<scalar,vector>;
template class BlockLduSystem<vector,scalar>;
template class BlockLduSystem<vector,vector>;
template class BlockLduSystem<vector,tensor>;
template class BlockLduSystem<tensor,tensor>;
template class BlockLduSystem<symmTensor,symmTensor>;
template class BlockLduSystem<sphericalTensor,sphericalTensor>;
template class BlockLduSystem<vector4,vector4>;
}

// ************************************************************************* //

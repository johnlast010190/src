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
    (c) 2011-2014 OpenFOAM Foundation

Typedef
    Foam::tensor

Description
    Tensor of scalars.

SourceFiles
    tensor.C

\*---------------------------------------------------------------------------*/

#ifndef tensor_H
#define tensor_H

#include "primitives/Tensor/Tensor.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "primitives/contiguous/contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Tensor<scalar> tensor;

vector eigenValues(const tensor& T);
vector eigenVector
(
    const tensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
);
tensor eigenVectors(const tensor& T, const vector& lambdas);
tensor eigenVectors(const tensor& T);

vector eigenValues(const symmTensor& T);
vector eigenVector
(
    const symmTensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
);
tensor eigenVectors(const symmTensor& T, const vector& lambdas);
tensor eigenVectors(const symmTensor& T);

// Tensor inverse with singular value decomposition
// (Householder transformations)
// HJ, 24/Oct/2009
tensor hinv(const tensor& t);
symmTensor hinv(const symmTensor& t);


//- Data associated with tensor type are contiguous
template<>
constexpr bool contiguous<tensor>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

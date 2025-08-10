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
    (c) 2004-6 H. Jasak All rights reserved

Class
    BlockLduMatrix

Description
    Template specialisation for sphericalTensor block matrix

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    sphericalTensorBlockLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalTensorBlockLduMatrix_H
#define sphericalTensorBlockLduMatrix_H

#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockLduMatrix<sphericalTensor>::sumDiag();

template<>
void BlockLduMatrix<sphericalTensor>::negSumDiag();

template<>
void BlockLduMatrix<sphericalTensor>::check() const;

template<>
void BlockLduMatrix<sphericalTensor>::relax
(
    const sphericalTensorField& x,
    sphericalTensorField& b,
    const scalar alpha
);

template<>
void BlockLduMatrix<sphericalTensor>::operator*=(const scalarField& sf);

template<>
void BlockLduMatrix<sphericalTensor>::AmulCore
(
    sphericalTensorField& mul,
    const sphericalTensorField& x
) const;

template<>
void BlockLduMatrix<sphericalTensor>::TmulCore
(
    sphericalTensorField& mul,
    const sphericalTensorField& x
) const;

template<>
void BlockLduMatrix<sphericalTensor>::segregateB
(
    sphericalTensorField& mul,
    const sphericalTensorField& x
) const;


template<>
tmp<sphericalTensorField>
BlockLduMatrix<sphericalTensor>::H(const sphericalTensorField& x) const;

template<>
tmp<sphericalTensorField>
BlockLduMatrix<sphericalTensor>::faceH(const sphericalTensorField& x) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduMatrix/tensorBlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::BlockLduMatrix<Foam::tensor>::sumDiag()
{
    // Decoupled version
    this->decoupledSumDiag();
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::negSumDiag()
{
    // Decoupled version
    this->decoupledNegSumDiag();
}


template<>
void Foam::BlockLduMatrix<tensor>::check() const
{
    // Decoupled version
    this->decoupledCheck();
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::relax
(
    const tensorField& x,
    tensorField& b,
    const scalar alpha
)
{
    // Decoupled version
    this->decoupledRelax(x, b, alpha);
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::operator*=(const scalarField& sf)
{
    // Decoupled version
    this->decoupledMultEqOp(sf);
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::AmulCore
(
    tensorField& Ax,
    const tensorField& x
) const
{
    decoupledAmulCore(Ax, x);
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::TmulCore
(
    tensorField& Tx,
    const tensorField& x
) const
{
    // Decoupled version
    decoupledTmulCore(Tx, x);
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::segregateB
(
    tensorField&,
    const tensorField&
) const
{
    FatalErrorInFunction
        << "Requested decoupling of tensor matrix - never coupled"
        << abort(FatalError);
}


template<>
Foam::tmp<Foam::tensorField>
Foam::BlockLduMatrix<Foam::tensor>::H(const tensorField& x) const
{
    // Decoupled version
    return decoupledH(x);
}


template<>
Foam::tmp<Foam::tensorField>
Foam::BlockLduMatrix<Foam::tensor>::faceH(const tensorField& x) const
{
    // Decoupled version
    return decoupledFaceH(x);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

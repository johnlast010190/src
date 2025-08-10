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

#include "matrices/blockLduMatrix/BlockLduMatrix/scalarBlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockLduMatrix<scalar>::sumDiag()
{
    scalarField& activeDiag = this->diag();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    if (symmetric())
    {
        const scalarField& activeUpper = this->upper();

        for (label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiag[l[coeffI]] += activeUpper[coeffI];
            activeDiag[u[coeffI]] += activeUpper[coeffI];
        }
    }
    else if (asymmetric())
    {
        const scalarField& activeLower = this->lower();
        const scalarField& activeUpper = this->upper();

        for (label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiag[l[coeffI]] += activeLower[coeffI];
            activeDiag[u[coeffI]] += activeUpper[coeffI];
        }
    }
    else
    {
        FatalErrorInFunction
            << "No off-diagonal available"
            << abort(FatalError);
    }
}


template<>
void BlockLduMatrix<scalar>::negSumDiag()
{
    scalarField& activeDiag = this->diag();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    if (symmetric())
    {
        const scalarField& activeUpper = this->upper();

        for (label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiag[l[coeffI]] -= activeUpper[coeffI];
            activeDiag[u[coeffI]] -= activeUpper[coeffI];
        }
    }
    else if (asymmetric())
    {
        const scalarField& activeLower = this->lower();
        const scalarField& activeUpper = this->upper();

        for (label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiag[l[coeffI]] -= activeLower[coeffI];
            activeDiag[u[coeffI]] -= activeUpper[coeffI];
        }
    }
    else
    {
        FatalErrorInFunction
            << "No off-diagonal available"
            << abort(FatalError);
    }
}


template<>
void BlockLduMatrix<scalar>::check() const
{
    // Copy the diagonal
    scalarField activeDiagCopy = this->diag();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    if (symmetric())
    {
        const scalarField& activeUpper = this->upper();

        for (label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiagCopy[l[coeffI]] -= activeUpper[coeffI];
            activeDiagCopy[u[coeffI]] -= activeUpper[coeffI];
        }

        Info<< "void BlockLduMatrix<scalar>::check() const : "
            << "Symmetric matrix: raw matrix difference: "
            << sum(mag(activeDiagCopy))
            << " scaled: "
            << sum(mag(activeDiagCopy))/sum(mag(this->diag()))
            << endl;
    }
    else if (asymmetric())
    {
        const scalarField& activeLower = this->lower();
        const scalarField& activeUpper = this->upper();

        for (label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiagCopy[l[coeffI]] -= activeLower[coeffI];
            activeDiagCopy[u[coeffI]] -= activeUpper[coeffI];
        }

        Info<< "void BlockLduMatrix<scalar>::check() const : "
            << "Asymmetric matrix: raw matrix difference: "
            << sum(mag(activeDiagCopy))
            << " scaled: "
            << sum(mag(activeDiagCopy))/sum(mag(this->diag()))
            << endl;
    }
    else
    {
        Info<< "void BlockLduMatrix<scalar>::check() const : "
            << "Diagonal matrix" << endl;
    }
}


template<>
void BlockLduMatrix<scalar>::relax
(
    const scalarField& x,
    scalarField& b,
    const scalar alpha
)
{
    scalarField& activeDiag = this->diag();

    scalarField activeDiagOld = this->diag();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    scalarField sumOff(activeDiag.size(), 0.0);

    if (symmetric())
    {
        const scalarField& activeUpper = this->upper();

        for (label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            sumOff[u[coeffI]] += mag(activeUpper[coeffI]);
            sumOff[l[coeffI]] += mag(activeUpper[coeffI]);
        }
    }
    else if (asymmetric())
    {
        const scalarField& activeLower = this->lower();
        const scalarField& activeUpper = this->upper();

        for (label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            sumOff[u[coeffI]] += mag(activeLower[coeffI]);
            sumOff[l[coeffI]] += mag(activeUpper[coeffI]);
        }
    }

    activeDiag = max(activeDiag, sumOff);
    activeDiag *= 1.0/alpha;

    // Add the relaxation contribution to b
    b += (activeDiag - activeDiagOld)*x;
}


template<>
void BlockLduMatrix<scalar>::operator*=(const scalarField& sf)
{
    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    if (upperPtr_)
    {
        scalarField& activeUpper = *upperPtr_;

        const labelUList& l = lduAddr().lowerAddr();

        for (label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeUpper[coeffI] *= sf[l[coeffI]];
        }
    }

    if (lowerPtr_)
    {
        scalarField& activeLower = *lowerPtr_;

        const labelUList& u = lduAddr().upperAddr();

        for (label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            activeLower[coeffI] *= sf[u[coeffI]];
        }
    }
}


template<>
void BlockLduMatrix<scalar>::AmulCore
(
    scalarField& Ax,
    const scalarField& x
) const
{
    const scalarField& Diag = diag();

    for (label rowI = 0; rowI < x.size(); rowI++)
    {
        Ax[rowI] = Diag[rowI]*x[rowI];
    }

    // Note: pointer looping
    const label* const __restrict__ U = lduAddr().upperAddr().begin();
    const label* const __restrict__ L = lduAddr().lowerAddr().begin();

    if (symmetric())
    {
        const scalar* const __restrict__ Upper = upper().begin();
        const scalar* const __restrict__ X = x.begin();

        scalar* __restrict__ AX = Ax.begin();

        for (label coeffI = 0; coeffI < upper().size(); coeffI++)
        {
            AX[U[coeffI]] += Upper[coeffI]*X[L[coeffI]];
            AX[L[coeffI]] += Upper[coeffI]*X[U[coeffI]];
        }
    }
    else if (asymmetric())
    {
        const scalar* const __restrict__ Lower = lower().begin();
        const scalar* const __restrict__  Upper = upper().begin();
        const scalar* const __restrict__ X = x.begin();

        scalar* __restrict__ AX = Ax.begin();

        for (label coeffI = 0; coeffI < upper().size(); coeffI++)
        {
            AX[U[coeffI]] += Lower[coeffI]*X[L[coeffI]];
            AX[L[coeffI]] += Upper[coeffI]*X[U[coeffI]];
        }
    }
}


template<>
void BlockLduMatrix<scalar>::TmulCore
(
    scalarField& Tx,
    const scalarField& x
) const
{
    const scalarField& Diag = diag();

    for (label rowI = 0; rowI < x.size(); rowI++)
    {
        Tx[rowI] = Diag[rowI]*x[rowI];
    }

    // Note: pointer looping
    const label* const __restrict__ U = lduAddr().upperAddr().begin();
    const label* const __restrict__ L = lduAddr().lowerAddr().begin();

    if (symmetric())
    {
        const scalar* const __restrict__ Upper = upper().begin();
        const scalar* const __restrict__ X = x.begin();

        scalar* __restrict__ TX = Tx.begin();

        for (label coeffI = 0; coeffI < upper().size(); coeffI++)
        {
            TX[U[coeffI]] += Upper[coeffI]*X[L[coeffI]];
            TX[L[coeffI]] += Upper[coeffI]*X[U[coeffI]];
        }
    }
    else if (asymmetric())
    {
        const scalar* const __restrict__ Lower = lower().begin();
        const scalar* const __restrict__ Upper = upper().begin();
        const scalar* const __restrict__ X = x.begin();

        scalar* __restrict__ TX = Tx.begin();

        for (label coeffI = 0; coeffI < upper().size(); coeffI++)
        {
            TX[U[coeffI]] += Upper[coeffI]*X[L[coeffI]];
            TX[L[coeffI]] += Lower[coeffI]*X[U[coeffI]];
        }
    }
}


template<>
void BlockLduMatrix<scalar>::segregateB
(
    scalarField&,
    const scalarField&
) const
{
    FatalErrorInFunction
        << "Requested decoupling of scalar matrix - never coupled"
        << abort(FatalError);
}


template<>
tmp<scalarField> BlockLduMatrix<scalar>::H(const scalarField& x) const
{
    tmp<scalarField> tresult(new scalarField(lduAddr().size(), 0));
    scalarField& result = tresult.ref();

    if (thereIsUpper() || thereIsLower())
    {
        // Note: pointer looping
        const label* const __restrict__ U = lduAddr().upperAddr().begin();
        const label* const __restrict__ L = lduAddr().lowerAddr().begin();

        const scalar* const __restrict__ Lower = lower().begin();
        const scalar* const __restrict__ Upper = upper().begin();
        const scalar* const __restrict__ X = x.begin();

        scalar* __restrict__ R = result.begin();

        for (label coeffI = 0; coeffI < upper().size(); coeffI++)
        {
            R[U[coeffI]] -= Upper[coeffI]*X[U[coeffI]];
        }

        for (label coeffI = 0; coeffI < upper().size(); coeffI++)
        {
            R[L[coeffI]] -= Lower[coeffI]*X[L[coeffI]];
        }
    }

    return tresult;
}


template<>
tmp<scalarField> BlockLduMatrix<scalar>::faceH(const scalarField& x) const
{
    tmp<scalarField> tresult(new scalarField(lduAddr().upperAddr().size(), 0));
    scalarField& result = tresult.ref();

    if (thereIsUpper() || thereIsLower())
    {
        // Note: pointer looping
        const label* const __restrict__ U = lduAddr().upperAddr().begin();
        const label* const __restrict__ L = lduAddr().lowerAddr().begin();

        const scalar* const __restrict__ Lower = lower().begin();
        const scalar* const __restrict__ Upper = upper().begin();
        const scalar* const __restrict__ X = x.begin();

        scalar* __restrict__ R = result.begin();

        for (label coeffI = 0; coeffI < upper().size(); coeffI++)
        {
            R[coeffI] = Upper[coeffI]*X[U[coeffI]] - Lower[coeffI]*X[L[coeffI]];
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

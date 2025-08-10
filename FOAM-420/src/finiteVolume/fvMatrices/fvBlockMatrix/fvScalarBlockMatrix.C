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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMatrices/fvBlockMatrix/fvScalarBlockMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
template<class Type2>
void Foam::fvBlockMatrix<Foam::scalar>::insertEquation
(
    const direction& dir,
    fvMatrix<Type2>& matrix
)
{
    FatalErrorInFunction
        << "Trying to assemble a scalar matrix using a scalarFvatrix."
        << "Only scalar matrix is allowed."
        << abort(FatalError);
}


//- specializaton of creating a scalar block matrix from a sergegated matrix
//  Probably non needed. Segregated is sufficient for scalar
template<>
template<>
void Foam::fvBlockMatrix<Foam::scalar>::insertEquation
(
    const direction& dir,
    fvMatrix<Foam::scalar>& matrix
)
{
    const volScalarField& psi = matrix.psi();

    insertSolutionVector(dir, psi.primitiveField());

    insertDiagSource(dir, matrix);

    //- upper and lower
    if (!matrix.diagonal())
    {
        if (matrix.hasUpper())
        {
            const scalarField& mUpper = matrix.upper();
            scalarField& upper = this->upper().asLinear();
            upper = mUpper;
        }
        if (matrix.symmetric() && this->symmetric())
        {
            bool resetField = false;
            if (!this->thereIsLower())
            {
                resetField = true;
            }
            const scalarField& mLower = matrix.lower();
            scalarField& lower = this->lower().asLinear();

            if (resetField)
            {
                lower = 0;
            }
            lower = mLower;
        }
    }


    //- coupled BCs
    forAll(psi.boundaryField(), patchI)
    {
        forAll(psi.boundaryField(), patchI)
        {
            const fvPatchField<scalar>& pf = psi.boundaryField()[patchI];
            const fvPatch& patch = pf.patch();

            if (patch.coupled())
            {
                scalarField& pcoupleUpper =
                    this->coupleUpper()[patchI].asLinear();

                scalarField& pcoupleLower =
                    this->coupleLower()[patchI].asLinear();

                const scalarField& icp = matrix.internalCoeffs()[patchI];
                const scalarField& bcp = matrix.boundaryCoeffs()[patchI];

                pcoupleLower += icp;
                pcoupleUpper += bcp;
            }
        }
    }

    //- Store BCs for fluxes
    if (psi_.mesh().schemes().fluxRequired(matrix.psi().name()))
    {
        const FieldField<Field, scalar>& iC = matrix.internalCoeffs();
        const FieldField<Field, scalar>& bC = matrix.boundaryCoeffs();

        forAll(this->internalCoeffs_, pI)
        {
            if (!this->internalCoeffs_.set(pI))
            {
                this->internalCoeffs_.set
                (
                    pI,
                    new CoeffField<scalar>
                    (
                        iC[pI].size()
                    )
                );
            }
            if (!this->boundaryCoeffs_.set(pI))
            {
                this->boundaryCoeffs_.set
                (
                    pI,
                    new Field<scalar>
                    (
                        bC[pI].size(),
                        0
                    )
                );
            }

            Field<scalar>& blockB = this->boundaryCoeffs_[pI];

            scalarField& blockI = this->internalCoeffs_[pI].asLinear();

            blockI += iC[pI];
            blockB += iC[pI];
        }
        if (matrix.faceFluxCorrectionPtr())
        {
            setFaceFluxCorrection(*(matrix.faceFluxCorrectionPtr()), dir);
        }
    }

    checkAndSetDimensions(dir, 1, matrix.dimensions());
}



// ************************************************************************* //

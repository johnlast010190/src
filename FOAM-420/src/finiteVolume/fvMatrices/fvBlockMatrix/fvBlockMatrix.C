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
    (c) 2011-2016 OpenFOAM Foundation
    (c) Vuko Vukcevic, FMENA Zagreb.
    (c) Update by Hrvoje Jasak
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"
#include "db/IOstreams/IOstreams.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "fields/fvPatchFields/derived/fixedJumpAMI/fixedJumpAMIFvPatchField.H"
#include "fvMesh/fvPatches/derived/inactive/inactiveFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class Type>
template<class fieldType>
void Foam::fvBlockMatrix<Type>::insertSolutionVector
(
    const direction dir,
    const Field<fieldType>& xSingle
)
{
    // Get number of field components and local copy of direction, for
    // consistency with member functions where directions need to be reset.
    const direction nCmpts = pTraits<fieldType>::nComponents;
    direction localDir = dir;

    Field<Type>& psiIn = psi_.primitiveFieldRef();

    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        psiIn.replace(localDir, xSingle.component(cmptI));

        localDir++;
    }
}


template<class Type>
template<class fieldType>
void Foam::fvBlockMatrix<Type>::modifyBCCompSpecificMembers
(
    const direction dir,
    const GeometricField<fieldType, fvPatchField, volMesh>& matrixPsi
)
{
    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& bFields = psi_.boundaryFieldRef();

    const typename GeometricField<fieldType, fvPatchField, volMesh>::
        Boundary& bFieldsMatrix = matrixPsi.boundaryField();

    forAll(bFieldsMatrix, patchi)
    {
        if ((isA<fixedJumpAMIFvPatchField<fieldType>>(bFieldsMatrix[patchi])))
        {
            const fixedJumpAMIFvPatchField<fieldType>& matrixpft =
                refCast<const fixedJumpAMIFvPatchField<fieldType>>(bFieldsMatrix[patchi]);

            fixedJumpAMIFvPatchField<Type>& pft =
                refCast<fixedJumpAMIFvPatchField<Type>>(bFields[patchi]);

            if (pft.cyclicAMIPatch().owner())
            {
                tmp<Field<fieldType>> tmpJump(matrixpft.jump());
                Field<Type>& jump = pft.jumpRef();

                const direction nCmpts = pTraits<fieldType>::nComponents;
                direction localDir = dir;
                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    jump.replace(localDir, tmpJump().component(cmptI));

                    localDir++;
                }
            }
        }
    }
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::insertDiagSource
(
    const direction dir,
    fvMatrix<Type2>& matrix
)
{
    applyBoundaryManipulate(matrix);

    // Save a copy for different components
    scalarField diag = matrix.diag();
    const scalarField& saveDiag = matrix.diag();

    // Add source boundary contribution
    Field<Type2> source = matrix.source();
    matrix.addBoundarySource(source, matrix.psi(), false);

    const direction nCmpts = pTraits<Type>::nComponents;
    const direction nCmpts2 = pTraits<Type2>::nComponents;
    direction localDir = dir;

    // Get reference to this source field of block system
    Field<Type>& b = this->source();

    if
    (
        this->diag().activeType() != blockCoeffBase::SQUARE
    )
    {
        typename CoeffField<Type>::linearTypeField& blockDiag =
            this->diag().asLinear();

        for (direction cmptI = 0; cmptI < nCmpts2; cmptI++)
        {

            matrix.addBoundaryDiag(diag, cmptI);
            diag += blockDiag.component(localDir);

            scalarField sourceCmpt
            (
                source.component(cmptI) + b.component(localDir)
            );

            blockDiag.replace(localDir, diag);
            b.replace(localDir, sourceCmpt);

            localDir++;

            // Reset diagonal
            diag = saveDiag;
        }
    }
    else if (this->diag().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<Type>::squareTypeField& blockDiag =
            this->diag().asSquare();

        for (direction cmptI = 0; cmptI < nCmpts2; cmptI++)
        {
            const direction tensorDir = nCmpts*localDir + localDir;

            matrix.addBoundaryDiag(diag, cmptI);
            diag += blockDiag.component(tensorDir);

            scalarField sourceCmpt
            (
                source.component(cmptI) + b.component(localDir)
            );

            blockDiag.replace(tensorDir, diag);
            b.replace(localDir, sourceCmpt);

            localDir++;

            // Reset diagonal
            diag = saveDiag;
        }
    }
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::applyBoundaryManipulate
(
    fvMatrix<Type2>& matrix
)
{
    bool gibExists = false;
    const fvMesh& fvMesh = psi_.mesh();
    forAll(fvMesh.boundary(), pI)
    {
        const fvPatch& fvp =  psi_.mesh().boundary()[pI];
        if (isA<inactiveFvPatch>(fvp))
        {
            gibExists = true;
            break;
        }
    }
    if (gibExists)
    {
        //-To check
        GeometricField<Type2, fvPatchField, volMesh>& psiMatrix =
           const_cast
           <
               GeometricField<Type2, fvPatchField, volMesh>&
           >(matrix.psi());
        matrix.boundaryManipulate(psiMatrix.boundaryFieldRef());
    }
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::insertDiagSource
(
    const direction dir,
    const fvBlockMatrix<Type2>& matrix
)
{
    const direction nCmpts = pTraits<Type2>::nComponents;

    this->diag().expandCoeffs(matrix.diag());

    if (matrix.diag().activeType() != blockCoeffBase::UNALLOCATED)
    {
        if (this->diag().activeType() != blockCoeffBase::SQUARE)
        {
            typename CoeffField<Type>::linearTypeField& blockDiag =
                this->diag().asLinear();

            const typename CoeffField<Type2>::linearTypeField& matrixDiag =
                matrix.diag().asLinear();


            forAll(blockDiag, cI)
            {
                label localDir1 = dir;
                for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    blockDiag[cI](localDir1) += matrixDiag[cI](cmptI);
                    this->source()[cI](localDir1) += matrix.source()[cI](cmptI);
                    localDir1++;
                }
            }
        }
        else if (this->diag().activeType() == blockCoeffBase::SQUARE)
        {
            typename CoeffField<Type>::squareTypeField& blockDiag =
                this->diag().asSquare();

            const typename CoeffField<Type2>::squareTypeField& matrixDiag =
                matrix.diag().asSquare();

            Field<Type>& b = this->source();
            const Field<Type2>& matrixSource =  matrix.source();

            forAll(blockDiag, cI)
            {
                label localDir1 = dir;
                for (label cmpt1I = 0; cmpt1I < nCmpts; cmpt1I++)
                {
                    label localDir2 = dir;

                    for (label cmpt2I = 0; cmpt2I < nCmpts; cmpt2I++)
                    {
                        blockDiag[cI](localDir1, localDir2) +=
                            matrixDiag[cI](cmpt1I, cmpt2I);
                        localDir2++;
                    }

                    b[cI](localDir1) += matrixSource[cI](cmpt1I);

                    localDir1++;
                }
            }
        }
    }
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::insertUpperLower
(
    const direction dir,
    const fvMatrix<Type2>& matrix
)
{
    if (matrix.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    const direction nCmpts = pTraits<Type2>::nComponents;

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();

        /*
        if (this->upper().activeType() == blockCoeffBase::UNALLOCATED)
        {
            this->upper().asScalar() = upper;
        }
        else
        */
        if
        (
            this->upper().activeType() == blockCoeffBase::UNALLOCATED
         || this->upper().activeType() == blockCoeffBase::SCALAR
         || this->upper().activeType() == blockCoeffBase::LINEAR
        )
        {
            typename CoeffField<Type>::linearTypeField& blockUpper =
                this->upper().asLinear();

            forAll(upper, fI)
            {
                label localDir = dir;
                for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    blockUpper[fI](localDir) += upper[fI];
                    localDir++;
                }
            }
        }
        else if (this->upper().activeType() == blockCoeffBase::SQUARE)
        {
            typename CoeffField<Type>::squareTypeField& blockUpper =
                this->upper().asSquare();

            forAll(upper, fI)
            {
                label localDir = dir;
                for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    blockUpper[fI](localDir, localDir) += upper[fI];
                    localDir++;
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Error in matrix insertion: problem with block structure."
            << abort(FatalError);
    }

    if (matrix.symmetric() && this->symmetric())
    {
        Info<< "Both matrices are symmetric: inserting only upper triangle"
            << endl;
    }
    else
    {
        bool resetField = false;
        if (!this->thereIsLower())
        {
            resetField = true;
        }

        // Either scalar or block matrix is asymmetric: insert lower triangle
        const scalarField& lower = matrix.lower();
/*
        if (this->lower().activeType() == blockCoeffBase::UNALLOCATED)
        {
            this->lower().asScalar() = lower;
        }
        else
        */
        if
        (
            this->lower().activeType() == blockCoeffBase::UNALLOCATED
         || this->lower().activeType() == blockCoeffBase::SCALAR
         || this->lower().activeType() == blockCoeffBase::LINEAR
        )
        {
            typename CoeffField<Type>::linearTypeField& blockLower =
                this->lower().asLinear();

            if (resetField)
            {
                blockLower =
                    pTraits<typename CoeffField<Type>::linearType>::zero;
            }

            forAll(lower, fI)
            {
                label localDir = dir;
                for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    blockLower[fI](localDir) += lower[fI];
                    localDir++;
                }
            }
        }
        else if (this->lower().activeType() == blockCoeffBase::SQUARE)
        {
            typename CoeffField<Type>::squareTypeField& blockLower =
                this->lower().asSquare();

            if (resetField)
            {
                blockLower =
                    pTraits<typename CoeffField<Type>::squareType>::zero;
            }

            forAll(lower, fI)
            {
                label localDir = dir;
                for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    blockLower[fI](localDir, localDir) += lower[fI];
                    localDir++;
                }
            }
        }
    }
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::insertUpperLower
(
    const direction dir,
    fvBlockMatrix<Type2>& matrix
)
{

    if (matrix.diagonal())
    {
        return;
    }

    this->upper().expandCoeffs(matrix.upper());

    const direction nCmpts = pTraits<Type2>::nComponents;

    if
    (
        this->upper().activeType() == blockCoeffBase::SCALAR
     || this->upper().activeType() == blockCoeffBase::LINEAR
    )
    {

        typename CoeffField<Type>::linearTypeField& thisUpper =
            this->upper().asLinear();

        typename CoeffField<Type2>::linearTypeField& matrixUpper =
            matrix.upper().asLinear();

        forAll(thisUpper, fI)
        {
            label localDir = dir;
            for (label cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                thisUpper[fI](localDir) += matrixUpper[fI](cmptI);
                localDir++;
            }
        }
    }
    else if (this->upper().activeType() == blockCoeffBase::SQUARE)
    {

        typename CoeffField<Type>::squareTypeField& thisUpper =
            this->upper().asSquare();

        typename CoeffField<Type2>::squareTypeField& matrixUpper =
            matrix.upper().asSquare();

        forAll(thisUpper, fI)
        {
            label localDir1 = dir;
            for (label cmpt1I = 0; cmpt1I < nCmpts; cmpt1I++)
            {
                label localDir2 = dir;
                for (label cmpt2I = 0; cmpt2I < nCmpts; cmpt2I++)
                {
                    thisUpper[fI](localDir1, localDir2) +=
                        matrixUpper[fI](cmpt1I, cmpt2I);
                    localDir2++;
                }
                localDir1++;
            }
        }
    }

    if (matrix.symmetric() && this->symmetric())
    {
        Info<< "Both matrices are symmetric: inserting only upper triangle"
            << endl;
    }
    else
    {
        bool resetField = false;
        if (!this->thereIsLower())
        {
            resetField = true;
        }

        this->lower().expandCoeffs(matrix.lower());

        if
        (
            this->lower().activeType() == blockCoeffBase::SCALAR
         || this->lower().activeType() == blockCoeffBase::LINEAR
        )
        {

            typename CoeffField<Type>::linearTypeField& thisLower =
                this->lower().asLinear();

            typename CoeffField<Type2>::linearTypeField& matrixLower =
                matrix.lower().asLinear();

            if (resetField)
            {
                thisLower =
                    pTraits<typename CoeffField<Type>::linearType>::zero;
            }


            forAll(thisLower, fI)
            {
                label localDir = dir;

                for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    thisLower[fI](localDir) += matrixLower[fI](cmptI);
                    localDir++;
                }
            }
        }
        else if (this->lower().activeType() == blockCoeffBase::SQUARE)
        {

            typename CoeffField<Type>::squareTypeField& thisLower =
                this->lower().asSquare();

            typename CoeffField<Type2>::squareTypeField& matrixLower =
                matrix.lower().asSquare();

            if (resetField)
            {
                thisLower =
                    pTraits<typename CoeffField<Type>::squareType>::zero;
            }


            forAll(thisLower, fI)
            {
                label localDir1 = dir;

                for (label cmpt1I = 0; cmpt1I < nCmpts; cmpt1I++)
                {
                    label localDir2 = dir;

                    for (label cmpt2I = 0; cmpt2I < nCmpts; cmpt2I++)
                    {
                        thisLower[fI](localDir1, localDir2) +=
                            matrixLower[fI](cmpt1I, cmpt2I);

                        localDir2++;
                    }

                    localDir1++;
                }
            }
        }
    }
}



template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::updateCouplingCoeffs
(
    const direction dir,
    const fvMatrix<Type2>& matrix
)
{

    if (matrix.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    const direction nCmpts = pTraits<Type2>::nComponents;
    direction localDir = dir;

    const GeometricField<Type2, fvPatchField, volMesh>& psi =
        matrix.psi();

    forAll(psi.boundaryField(), patchI)
    {
        const fvPatchField<Type2>& pf = psi.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();

        if (patch.coupled())
        {
            const Field<Type2>& icp = matrix.internalCoeffs()[patchI];
            const Field<Type2>& bcp = matrix.boundaryCoeffs()[patchI];

            if
            (
                this->coupleUpper()[patchI].activeType()
             != blockCoeffBase::SQUARE
            )
            {
                typename CoeffField<Type>::linearTypeField& pcoupleUpper =
                    this->coupleUpper()[patchI].asLinear();

                typename CoeffField<Type>::linearTypeField& pcoupleLower =
                    this->coupleLower()[patchI].asLinear();

                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    scalarField icpCmpt ( icp.component(cmptI) );
                    scalarField bcpCmpt ( bcp.component(cmptI) );

                    forAll(pf, fI)
                    {
                        pcoupleUpper[fI](localDir) += bcpCmpt[fI];
                        pcoupleLower[fI](localDir) += icpCmpt[fI];
                    }

                    localDir++;
                }

                localDir = dir;
            }
            else if
            (
                this->coupleUpper()[patchI].activeType()
              == blockCoeffBase::SQUARE
            )
            {
                typename CoeffField<Type>::squareTypeField& pcoupleUpper =
                    this->coupleUpper()[patchI].asSquare();

                typename CoeffField<Type>::squareTypeField& pcoupleLower =
                    this->coupleLower()[patchI].asSquare();

                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    scalarField icpCmpt ( icp.component(cmptI) );
                    scalarField bcpCmpt ( bcp.component(cmptI) );

                    forAll(pf, fI)
                    {
                        pcoupleUpper[fI](localDir, localDir) += bcpCmpt[fI];
                        pcoupleLower[fI](localDir, localDir) += icpCmpt[fI];
                    }

                    localDir++;
                }

                localDir = dir;
            }
        }
    }
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::updateCouplingCoeffs
(
    const direction dir,
    const fvBlockMatrix<Type2>& matrix
)
{
    if (matrix.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    const direction nCmpts = pTraits<Type2>::nComponents;

    const GeometricField<Type2, fvPatchField, volMesh>& psi =
        matrix.psi();

    forAll(psi.boundaryField(), pI)
    {
        const fvPatchField<Type2>& pf = psi.boundaryField()[pI];
        const fvPatch& patch = pf.patch();

        if (patch.coupled())
        {

            this->coupleUpper()[pI].expandCoeffs(matrix.coupleUpper()[pI]);

            if
            (
                this->coupleUpper()[pI].activeType()
                  != blockCoeffBase::SQUARE
            )
            {
                typename CoeffField<Type>::linearTypeField& thispcUpper =
                    this->coupleUpper()[pI].asLinear();

                typename CoeffField<Type>::linearTypeField& thispcLower =
                    this->coupleLower()[pI].asLinear();

                const typename CoeffField<Type2>::linearTypeField&
                    matrixpcUpper = matrix.coupleUpper()[pI].asLinear();

                const typename CoeffField<Type2>::linearTypeField&
                    matrixpcLower = matrix.coupleLower()[pI].asLinear();

                forAll(pf, fI)
                {
                    label localDir = dir;
                    for (label cmptI = 0; cmptI < nCmpts; cmptI++)
                    {
                        thispcUpper[fI](localDir) += matrixpcUpper[fI](cmptI);
                        thispcLower[fI](localDir) += matrixpcLower[fI](cmptI);

                        localDir++;
                    }
                }
            }
            else if
            (
                this->coupleUpper()[pI].activeType()
              == blockCoeffBase::SQUARE
            )
            {
                typename CoeffField<Type>::squareTypeField& thispcUpper =
                    this->coupleUpper()[pI].asSquare();

                typename CoeffField<Type>::squareTypeField& thispcLower =
                    this->coupleLower()[pI].asSquare();

                const typename CoeffField<Type2>::squareTypeField&
                    matrixpcUpper = matrix.coupleUpper()[pI].asSquare();

                const typename CoeffField<Type2>::squareTypeField&
                    matrixpcLower = matrix.coupleLower()[pI].asSquare();


                forAll(pf, fI)
                {
                    for (label cmpt1I = 0; cmpt1I < nCmpts; cmpt1I++)
                    {
                        label localDir1 = dir;

                        for (label cmpt2I = 0; cmpt2I < nCmpts; cmpt2I++)
                        {
                            label localDir2 = dir;

                            thispcUpper[fI](localDir1, localDir2) +=
                                matrixpcUpper[fI](cmpt1I, cmpt2I);
                            thispcLower[fI](localDir1, localDir2) +=
                                matrixpcLower[fI](cmpt1I, cmpt2I);

                            localDir2++;
                        }

                        localDir1++;
                    }
                }
            }
        }
    }
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::insertBoundaryCoeffsIfRequired
(
    const direction dir,
    const fvMatrix<Type2>& matrix
)
{
    if (psi_.mesh().schemes().fluxRequired(matrix.psi().name()))
    {
        const direction nCmpts = pTraits<Type2>::nComponents;

        const FieldField<Field, Type2>& iC = matrix.internalCoeffs();
        const FieldField<Field, Type2>& bC = matrix.boundaryCoeffs();

        forAll(this->internalCoeffs_, pI)
        {
            if (!this->internalCoeffs_.set(pI))
            {
                this->internalCoeffs_.set
                (
                    pI,
                    new CoeffField<Type>
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
                    new Field<Type>
                    (
                        bC[pI].size(),
                        pTraits<Type>::zero
                    )
                );
            }

            Field<Type>& blockB = this->boundaryCoeffs_[pI];


            if (this->diag().activeType() != blockCoeffBase::SQUARE)
            {
                typename CoeffField<Type>::linearTypeField& blockI =
                    this->internalCoeffs_[pI].asLinear();

                direction localDir = dir;

                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    scalarField iCpIcmntI ( iC[pI].component(cmptI) );
                    scalarField bCpIcmntI ( bC[pI].component(cmptI) );

                    forAll(blockI, fI)
                    {
                        blockI[fI](localDir) += iCpIcmntI[fI];
                        blockB[fI](localDir) += bCpIcmntI[fI];
                    }

                    localDir++;
                }
            }
            else if (this->diag().activeType() == blockCoeffBase::SQUARE)
            {
                typename CoeffField<Type>::squareTypeField& blockI =
                    this->internalCoeffs_[pI].asSquare();

                direction localDir = dir;

                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    scalarField iCpIcmntI ( iC[pI].component(cmptI) );
                    scalarField bCpIcmntI ( bC[pI].component(cmptI) );

                    forAll(blockI, fI)
                    {
                        blockI[fI](localDir, localDir) += iCpIcmntI[fI];
                        blockB[fI](localDir) += bCpIcmntI[fI];
                    }

                    localDir++;
                }
            }
        }
        if (matrix.faceFluxCorrectionPtr())
        {
            setFaceFluxCorrection(*(matrix.faceFluxCorrectionPtr()), dir);
        }
    }
}



template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::insertBoundaryCoeffsIfRequired
(
    const direction dir,
    const fvBlockMatrix<Type2>& bm
)
{
    const direction nCmpts = pTraits<Type2>::nComponents;

    const PtrList<CoeffField<Type2>>& iC = bm.internalCoeffs();
    const PtrList<Field<Type2>>& bC = bm.boundaryCoeffs();

    forAll(this->internalCoeffs_, pI)
    {
        if (iC.set(pI) && bC.set(pI))
        {
            if (!this->internalCoeffs_.set(pI))
            {
                this->internalCoeffs_.set
                (
                    pI,
                    new CoeffField<Type>
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
                    new Field<Type>
                    (
                        bC[pI].size(),
                        pTraits<Type>::zero
                    )
                );
            }

            this->internalCoeffs_[pI].expandCoeffs(iC[pI]);

            if (this->diag().activeType() != blockCoeffBase::SQUARE)
            {
                typename CoeffField<Type>::linearTypeField& blockI =
                    this->internalCoeffs_[pI].asLinear();

                Field<Type>& blockB = this->boundaryCoeffs_[pI];

                const typename CoeffField<Type2>::linearTypeField& iCI =
                    iC[pI].asLinear();

                const Field<Type2>& bCI = bC[pI];

                direction localDir = dir;

                for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                {
                    const word nameFieldI(bIndex_.variableName(localDir));
                    if (psi_.mesh().schemes().fluxRequired(nameFieldI))
                    {
                        forAll(blockI, fI)
                        {
                            blockI[fI](localDir) += iCI[fI](cmptI);
                            blockB[fI](localDir) += bCI[fI](cmptI);
                        }
                    }

                    localDir++;
                }
            }
            else if (this->diag().activeType() == blockCoeffBase::SQUARE)
            {
                typename CoeffField<Type>::squareTypeField& blockI =
                    this->internalCoeffs_[pI].asSquare();

                Field<Type>& blockB = this->boundaryCoeffs_[pI];

                const typename CoeffField<Type2>::squareTypeField& iCI =
                    iC[pI].asSquare();

                const Field<Type2>& bCI = bC[pI];

                direction localDir1 = dir;

                for (direction cmpt1I = 0; cmpt1I < nCmpts; cmpt1I++)
                {
                    const word nameFieldI(bIndex_.variableName(localDir1));
                    if (psi_.mesh().schemes().fluxRequired(nameFieldI))
                    {
                        direction localDir2 = dir;

                        for (direction cmpt2I = 0; cmpt2I < nCmpts; cmpt2I++)
                        {
                            forAll(blockI, fI)
                            {
                                blockI[fI](localDir1, localDir2) +=
                                    iCI[fI](cmpt1I, cmpt2I);
                            }
                            localDir2++;
                        }

                        forAll(blockB, fI)
                        {
                            blockB[fI](localDir1) += bCI[fI](cmpt1I);
                        }

                    }
                    localDir1++;
                }
            }
        }
    }
    if (psi_.mesh().schemes().fluxRequired(bIndex_.variableName(dir)))
    {
        if (bm.faceFluxCorrectionPtr())
        {
           // I have to implement this
            setFaceFluxCorrection(*(bm.faceFluxCorrectionPtr()), dir);
        }
    }
}


template<class Type>
template<class blockType, class fieldType>
void Foam::fvBlockMatrix<Type>::insertBoundaryCoeffsIfRequired
(
    const direction dirI,
    const direction dirJ,
    const BlockLduSystem<blockType, fieldType>& bs,
    const bool incrementColumnDir
)
{
    const PtrList<CoeffField<blockType>>& iC = bs.internalCoeffs();
    const PtrList<Field<fieldType>>& bC = bs.boundaryCoeffs();
    forAll(this->internalCoeffs_, pI)
    {
        if (iC.set(pI) && bC.set(pI))
        {
            if (!this->internalCoeffs_.set(pI))
            {
                this->internalCoeffs_.set
                (
                    pI,
                    new CoeffField<Type>
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
                    new Field<Type>
                    (
                        bC[pI].size(),
                        pTraits<Type>::zero
                    )
                );
            }

            const direction nCmpts = pTraits<blockType>::nComponents;
            const direction nCmpts2 = pTraits<fieldType>::nComponents;

            direction localDirI = dirI;
            direction localDirJ = dirJ;


            typename CoeffField<Type>::squareTypeField& blockI =
                this->internalCoeffs_[pI].asSquare();

            Field<Type>& blockB = this->boundaryCoeffs_[pI];

            const typename CoeffField<blockType>::linearTypeField& block2I =
                iC[pI].asLinear();

            const Field<fieldType>& block2B = bC[pI];

            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                if (psi_.mesh().schemes().fluxRequired(bIndex_.variableName(localDirI)))
                {
                    forAll(blockI, fI)
                    {
                        blockI[fI](localDirI, localDirJ) +=
                            block2I[fI](cmptI);
                    }
                }

                if (incrementColumnDir)
                {
                    localDirI++;
                }
                else
                {
                    localDirJ++;
                }
            }
            for (direction cmptI = 0; cmptI < nCmpts2; cmptI++)
            {
                if (psi_.mesh().schemes().fluxRequired(bIndex_.variableName(localDirI)))
                {
                    forAll(blockB, fI)
                    {
                        blockB[fI](localDirI) += block2B[fI](cmptI);
                    }
                }
            }
        }
        if (bs.faceFluxCorrectionPtr())
        {
            //- missing Implementation
            //setFaceFluxCorrection(*(bm.faceFluxCorrectionPtr()), dir);
        }
    }
}

template<class Type>
template<class blockType>
void Foam::fvBlockMatrix<Type>::insertBoundaryCoeffsIfRequired
(
    const direction dirI,
    const direction dirJ,
    const BlockLduSystem<blockType, scalar>& bs,
    const bool incrementColumnDir
)
{
    const PtrList<CoeffField<blockType>>& iC = bs.internalCoeffs();
    const PtrList<Field<scalar>>& bC = bs.boundaryCoeffs();
    forAll(this->internalCoeffs_, pI)
    {
        if (iC.set(pI) && bC.set(pI))
        {
            if (!this->internalCoeffs_.set(pI))
            {
                this->internalCoeffs_.set
                (
                    pI,
                    new CoeffField<Type>
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
                    new Field<Type>
                    (
                        bC[pI].size(),
                        pTraits<Type>::zero
                    )
                );
            }

            const direction nCmpts = pTraits<blockType>::nComponents;

            direction localDirI = dirI;
            direction localDirJ = dirJ;


            typename CoeffField<Type>::squareTypeField& blockI =
                this->internalCoeffs_[pI].asSquare();

            Field<Type>& blockB = this->boundaryCoeffs_[pI];

            const typename CoeffField<blockType>::linearTypeField& block2I =
                iC[pI].asLinear();

            const scalarField& block2B = bC[pI];

            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                if (psi_.mesh().schemes().fluxRequired(bIndex_.variableName(localDirI)))
                {
                    forAll(blockI, fI)
                    {
                        blockI[fI](localDirI, localDirJ) +=
                            block2I[fI](cmptI);
                    }
                }

                if (incrementColumnDir)
                {
                    localDirI++;
                }
                else
                {
                    localDirJ++;
                }
            }

            localDirI = dirI;

            if (psi_.mesh().schemes().fluxRequired(bIndex_.variableName(localDirI)))
            {
                forAll(blockB, fI)
                {
                    blockB[fI](localDirI) += block2B[fI];
                }
            }
        }
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::checkAndSetDimensions
(
    const direction dir,
    const direction nCmpts,
    const dimensionSet& dSet
)
{
    PtrList<dimensionSet>& dimSetList = this->dimensionSets();

    const direction nDir = dir+nCmpts;

    for (direction dI = dir; dI < nDir; dI++)
    {
        if (dimSetList.set(dI))
        {
            if (dimSetList[dI] != dSet)
            {
                FatalErrorInFunction
                    << "Inconsistent dimensions"
                    << dSet << " != "
                    << dimSetList[dI]
                    << abort(FatalError);
            }
        }
        else
        {
            dimSetList.set(dI, new dimensionSet(dSet));
        }
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::checkAndSetDimensions
(
    const direction dir,
    const PtrList<dimensionSet>& dSets
)
{
    PtrList<dimensionSet>& dimSetList = this->dimensionSets();

    forAll(dSets, dI)
    {
        if (dimSetList.set(dir+dI))
        {
            if (dimSetList[dir+dI] != dSets[dI])
            {
                FatalErrorInFunction
                    << "Inconsistent dimensions"
                    << dimSetList[dir+dI] << " != "
                    << dSets[dI]
                    << abort(FatalError);
            }
        }
        else
        {
            dimSetList.set(dir+dI, new dimensionSet(dSets[dI]));
        }
    }
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::addBoundaryFlux
(
    Foam::GeometricField<Type2, Foam::fvsPatchField, Foam::surfaceMesh>& field,
    const direction sI
) const
{
    typename GeometricField<Type2, fvsPatchField, surfaceMesh>::
        Boundary& ffbf = field.boundaryFieldRef();

    const direction nCmpts2 = pTraits<Type2>::nComponents;

    typename BlockCoeff<Type>::multiply mult;

    //- This function is not used currently.
    //  Only the scalar specialization of Type2=scalar which is the next one
    //  which is used in flux
    forAll(ffbf, pI)
    {
        tmp<Field<Type>> psiPIf(psi_.boundaryField()[pI].patchInternalField());

        const Field<Type>& blockB = this->boundaryCoeffs()[pI];

        label fieldSize = ffbf[pI].size();

        for (direction cmptI = 0; cmptI < nCmpts2; cmptI++)
        {
            Field<scalar> nC(this->boundaryCoeffs_[pI].size(), 0);
            Field<scalar> iC(this->internalCoeffs_[pI].size(), 0);

            if (this->internalCoeffs_[pI].activeType() != blockCoeffBase::SQUARE)
            {
                const typename CoeffField<Type>::linearTypeField& blockI =
                    this->internalCoeffs()[pI].asLinear();
                for (label fI = 0; fI < fieldSize; fI++)
                {
                    iC[fI] =  mult(sI, blockI[fI], psiPIf()[fI]);
                }
                if (psi_.boundaryField()[pI].coupled())
                {
                    tmp<Field<Type>> psiPNf
                    (
                        psi_.boundaryField()[pI].patchNeighbourField()
                    );
                    const typename CoeffField<Type>::linearTypeField& pcUpper =
                        this->coupleUpper()[pI].asLinear();

                    for (label fI = 0; fI < fieldSize; fI++)
                    {
                        nC[fI] = mult(sI, pcUpper[fI], psiPNf()[fI]);
                    }
                }
                else
                {
                    for (label fI = 0; fI < ffbf[pI].size(); fI++)
                    {
                        nC[fI] = blockB[fI](sI+cmptI);
                    }
                }
            }
            else
            {
                const typename CoeffField<Type>::squareTypeField& blockI =
                    this->internalCoeffs()[pI].asSquare();

                for (label fI = 0; fI < fieldSize; fI++)
                {
                    iC[fI] = mult(sI+cmptI, blockI[fI], psiPIf()[fI]);
                }
                if (psi_.boundaryField()[pI].coupled())
                {
                    tmp<Field<Type>> psiPNf
                    (
                        psi_.boundaryField()[pI].patchNeighbourField()
                    );

                    const typename CoeffField<Type>::squareTypeField& pcUpper =
                        this->coupleUpper()[pI].asSquare();

                    for (label fI = 0; fI < fieldSize; fI++)
                    {
                        nC[fI] = mult(sI, pcUpper[fI], psiPNf()[fI]);
                    }
                }
                else
                {
                    for (label fI = 0; fI < ffbf[pI].size(); fI++)
                    {
                        nC[fI] = blockB[fI](sI+cmptI);
                    }
                }
            }
            ffbf[pI].replace(cmptI, iC - nC);
        }
    }
}

template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::setFaceFluxCorrection
(
    const Foam::GeometricField<Type2, Foam::fvsPatchField, Foam::surfaceMesh>&
        matrixFluxCorr,
    const direction sI
)
{
    if (!this->faceFluxCorrectionPtr_)
    {
        constructFaceFluxCorrection();
    }
    Field<Type>& fluxCorrIn =
        (*this->faceFluxCorrectionPtr_).primitiveFieldRef();

    const Field<Type2>& mfluxCorrIn = matrixFluxCorr.primitiveField();

    const direction nCmpts2 = pTraits<Type2>::nComponents;

    for (direction cmptI = 0; cmptI < nCmpts2; cmptI++)
    {
        scalarField intField
        (
            mfluxCorrIn.component(sI+cmptI)+fluxCorrIn.component(cmptI)
        );

        fluxCorrIn.replace(sI+cmptI, intField);

        forAll(matrixFluxCorr.boundaryField(), pI)
        {
            Field<Type>& pFluxCorr  =
                (*this->faceFluxCorrectionPtr_).boundaryFieldRef()[pI];

            const Field<Type2>& pmFluxCorr = matrixFluxCorr.boundaryField()[pI];

            scalarField bField
            (
                pFluxCorr.component(sI+cmptI)
              + pmFluxCorr.component(cmptI)
            );

            pFluxCorr.replace(sI+cmptI, bField);
        }
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::constructFaceFluxCorrection()
{
    this->faceFluxCorrectionPtr_ =
        new GeometricField<Type, fvsPatchField, surfaceMesh>
    (
        IOobject
        (
            "fluxCorr("+psi_.name()+')',
            psi_.instance(),
            psi_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        psi_.mesh(),
        dimensioned<Type>("zero", dimless, pTraits<Type>::zero),
        fvsPatchField<Type>::calculatedType(),
        psi_.boundaryField().patchTypes()
    );
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::addFaceFluxCorrection
(
    Foam::GeometricField<Type2, Foam::fvsPatchField, Foam::surfaceMesh>&
        matrixFluxCorr,
    const direction sI
) const
{
    if (this->faceFluxCorrectionPtr_)
    {
        const Field<Type>& fluxCorrIn =
            (*this->faceFluxCorrectionPtr_).primitiveField();

        Field<Type2>& mfluxCorrIn = matrixFluxCorr.primitiveFieldRef();

        const direction nCmpts2 = pTraits<Type2>::nComponents;

        for (direction cmptI = 0; cmptI < nCmpts2; cmptI++)
        {
            scalarField intField
            (
                mfluxCorrIn.component(cmptI)
              + fluxCorrIn.component(sI+cmptI)
            );

            mfluxCorrIn.replace(cmptI, intField);

            forAll(matrixFluxCorr.boundaryField(), pI)
            {
                const Field<Type>& pFluxCorr  =
                    (*this->faceFluxCorrectionPtr_).boundaryField()[pI];

                Field<Type2>& pmFluxCorr  =
                    matrixFluxCorr.boundaryFieldRef()[pI];

                scalarField bField
                (
                    pFluxCorr.component(sI+cmptI)
                  + pmFluxCorr.component(cmptI)
                );

                pmFluxCorr.replace(cmptI, bField);
            }
        }
    }
}


template<class Type>
template<class blockType, class fieldType>
void Foam::fvBlockMatrix<Type>::insertBlock
(
    const direction dirI,
    const direction dirJ,
    const BlockLduSystem<blockType, fieldType>& blockSystem,
    const bool incrementColumnDir
)
{
    // Sanity checks
    {
        const direction blockMatrixSize = pTraits<blockType>::nComponents;
        const direction blockMatrixASize = pTraits<Type>::nComponents;

        if (blockMatrixSize > blockMatrixASize)
        {
            FatalErrorInFunction
                << "Trying to insert a block matrix from BlockLduSystem into "
                << "smaller one from fvBlockMatrix."
                << abort(FatalError);
        }

        if (dirI == dirJ)
        {
            FatalErrorInFunction
                << "Trying to insert coupling in the position where equation "
                << "should be, since dirI = dirJ. Try using insertEquation "
                << "member function."
                << abort(FatalError);
        }
    }

    const direction nCmpts = pTraits<blockType>::nComponents;

    // Get references to ldu fields of blockMatrix always as linear
    const typename CoeffField<blockType>::linearTypeField& bmd =
        blockSystem.diag().asLinear();
    const typename CoeffField<blockType>::linearTypeField& bmu =
        blockSystem.upper().asLinear();
    const typename CoeffField<blockType>::linearTypeField& bml =
        blockSystem.lower().asLinear();

    // Get references to ldu fields of this block matrix always as square
    typename CoeffField<Type>::squareTypeField& blockDiag =
         this->diag().asSquare();
    typename CoeffField<Type>::squareTypeField& blockUpper =
         this->upper().asSquare();
    typename CoeffField<Type>::squareTypeField& blockLower =
         this->lower().asSquare();

    // Insert blockMatrix that represents coupling into larger system matrix
    // Here the loop can be a bit faster
    // cells->components instead of components->cells

    forAll(bmd, cI)
    {
        label localDirI = dirI;
        label localDirJ = dirJ;
        for (label cmptI = 0; cmptI < nCmpts; cmptI++)
        {
            blockDiag[cI](localDirI, localDirJ) += bmd[cI].component(cmptI);
            if (incrementColumnDir)
            {
                localDirI++;
            }
            else
            {
                localDirJ++;
            }
        }
    }

    forAll(bmu, fI)
    {
        label localDirI = dirI;
        label localDirJ = dirJ;
        for (label cmptI = 0; cmptI < nCmpts; cmptI++)
        {
            blockUpper[fI](localDirI, localDirJ) += bmu[fI].component(cmptI);
            blockLower[fI](localDirI, localDirJ) += bml[fI].component(cmptI);
            if (incrementColumnDir)
            {
                localDirI++;
            }
            else
            {
                localDirJ++;
            }
        }
    }
}


template<class Type>
template<class blockType, class fieldType>
void Foam::fvBlockMatrix<Type>::insertBoundaryContributions
(
    const direction dirI,
    const direction dirJ,
    const BlockLduSystem<blockType, fieldType>& blockSystem,
    const bool incrementColumnDir
)
{
    // Need to get reference to fvMesh instead of lduMesh
    const fvBoundaryMesh& bmesh =
        refCast<const fvMesh>(this->mesh()).boundary();

    const direction nCmpts = pTraits<blockType>::nComponents;
    const direction nSrcCmpts = pTraits<fieldType>::nComponents;
    direction localDirI = dirI;
    direction localDirJ = dirJ;

    const Field<fieldType>& source = blockSystem.source();

    // Get reference to this source field of block system
    Field<Type>& b = this->source();

    // Insert source from block system to this system's rhs
    for (direction cmptI = 0; cmptI < nSrcCmpts; cmptI++)
    {
        scalarField sourceCmpt(source.component(cmptI));

        forAll(b, cI)
        {
            b[cI](localDirI) += sourceCmpt[cI];
        }

        if (incrementColumnDir)
        {
            localDirI++;
        }
        else
        {
            localDirJ++;
        }
    }

    // Reset local directions for coupling contributions
    localDirI = dirI;
    localDirJ = dirJ;

    // Insert coupling contributions into block matrix
    forAll(bmesh, patchI)
    {
        if (bmesh[patchI].coupled())
        {
            typename CoeffField<Type>::squareTypeField& pcoupleUpper =
                this->coupleUpper()[patchI].asSquare();

            typename CoeffField<Type>::squareTypeField& pcoupleLower =
                this->coupleLower()[patchI].asSquare();

            const typename CoeffField<blockType>::linearTypeField& bmcu =
                blockSystem.coupleUpper()[patchI].asLinear();

            const typename CoeffField<blockType>::linearTypeField& bmcl =
                blockSystem.coupleLower()[patchI].asLinear();

            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                forAll(bmcu, fI)
                {
                    pcoupleUpper[fI](localDirI, localDirJ) +=
                          bmcu[fI].component(cmptI);

                    pcoupleLower[fI](localDirI, localDirJ) +=
                          bmcl[fI].component(cmptI);
                }

                if (incrementColumnDir)
                {
                    localDirI++;
                }
                else
                {
                    localDirJ++;
                }
            }

            // Reset local directions for other patches
            localDirI = dirI;
            localDirJ = dirJ;
        }
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::insertCouplingDiag
(
    const direction dirI,
    const direction dirJ,
    const scalarField& coeffIJ
)
{
    // Get reference to block diagonal of the block system
    typename CoeffField<Type>::squareTypeField& blockDiag =
        this->diag().asSquare();

    // Set off-diagonal coefficient
    forAll(coeffIJ, cI)
    {
        blockDiag[cI](dirI, dirJ) += coeffIJ[cI];
    }

    // Source compensation is done in function updateSourceCoupling()
    // after all coupling terms are added.  HJ, 27/Apr/2015
}


template<class Type>
void Foam::fvBlockMatrix<Type>::insertCouplingDiagSource
(
    const direction dirI,
    const direction dirJ,
    const fvScalarMatrix& matrix
)
{
    // Prepare the diagonal and source
    scalarField diag = matrix.diag();
    scalarField source = matrix.source();

    // Add boundary source contribution
    matrix.addBoundaryDiag(diag, 0);
    matrix.addBoundarySource(source, matrix.psi(), false);

    // Get reference to block diagonal of the block system
    typename CoeffField<Type>::squareTypeField& blockDiag =
        this->diag().asSquare();

    // Get reference to this source field of the block system
    Field<Type>& b = this->source();

    // Set off-diagonal coefficient
    forAll(diag, cI)
    {
        blockDiag[cI](dirI, dirJ) += diag[cI];
        b[cI](dirI) += source[cI];
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::insertCouplingUpperLower
(
    const direction dirI,
    const direction dirJ,
    const fvScalarMatrix& matrix
)
{
    if (matrix.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    if (matrix.symmetric() && this->symmetric())
    {
        Info<< "Both fvScalarMatrix and block matrix are symmetric: " << nl
            << "inserting only upper triangle"
            << endl;
    }
    else
    {
        // Either scalar or block matrix is asymmetric: insert lower triangle
        const scalarField& lower = matrix.lower();

        typename CoeffField<Type>::squareTypeField& blockLower =
            this->lower().asSquare();

        forAll(lower, cI)
        {
            blockLower[cI](dirI, dirJ) += lower[cI];
        }
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();

        typename CoeffField<Type>::squareTypeField& blockUpper =
            this->upper().asSquare();

        forAll(upper, cI)
        {
            blockUpper[cI](dirI, dirJ) += upper[cI];
        }
    }
    else
    {
        FatalErrorInFunction
            << "Error in matrix insertion: problem with block structure."
            << abort(FatalError);
    }

    // Insert block interface fields
    forAll(this->interfaces(), patchI)
    {
        if (this->interfaces().set(patchI))
        {
            // Couple upper and lower
            const scalarField& cUpper = matrix.boundaryCoeffs()[patchI];
            const scalarField& cLower = matrix.internalCoeffs()[patchI];

            typename CoeffField<Type>::squareTypeField& blockUpper =
                this->coupleUpper()[patchI].asSquare();

            typename CoeffField<Type>::squareTypeField& blockLower =
                this->coupleLower()[patchI].asSquare();

            forAll(cUpper, fI)
            {
                blockUpper[fI](dirI, dirJ) += cUpper[fI];
                blockLower[fI](dirI, dirJ) += cLower[fI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvBlockMatrix<Type>::fvBlockMatrix
(
    GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    BlockLduSystem<Type, Type>(psi.mesh()),
    psi_(psi),
    bIndex_()
{
    this->interfaces() = psi.boundaryField().blockInterfaces();
}


template<class Type>
Foam::fvBlockMatrix<Type>::fvBlockMatrix
(
    GeometricField<Type, fvPatchField, volMesh>& psi,
    const dimensionSet& dSet
)
:
    BlockLduSystem<Type, Type>(psi.mesh(), dSet),
    psi_(psi),
    bIndex_()
{
    this->interfaces() = psi.boundaryField().blockInterfaces();
}


template<class Type>
Foam::fvBlockMatrix<Type>::fvBlockMatrix
(
    const fvBlockMatrix<Type>& bxs
)
:
    BlockLduSystem<Type, Type>(bxs),
    psi_
    (
        const_cast<GeometricField<Type, fvPatchField, volMesh>&>
        (bxs.psi())
    ),
    bIndex_(bxs.bIndex_)
{}


template<class Type>
Foam::fvBlockMatrix<Type>::fvBlockMatrix
(
    fvMatrix<Type>& matrix
)
:
    BlockLduSystem<Type, Type>(matrix.psi().mesh(), matrix.dimensions()),
    psi_
    (
        const_cast<GeometricField<Type, fvPatchField, volMesh>&>
        (matrix.psi())
    ),
    bIndex_()
{
    insertEquation(0, matrix);
}


template<class Type>
Foam::fvBlockMatrix<Type>::fvBlockMatrix
(
    tmp<fvMatrix<Type>> matrix
)
:
    BlockLduSystem<Type, Type>(matrix.ref().psi().mesh(), matrix.ref().dimensions()),
    psi_
    (
        const_cast<GeometricField<Type, fvPatchField, volMesh>&>
        (matrix.ref().psi())
    ),
    bIndex_()
{
    insertEquation(0, matrix.ref());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::setValuesFromList
(
    const word name,
    const labelUList& cellLabels,
    const Field<Type2>& values
)
{
    direction sI = bIndex_.findStartIndex(name);

    for (direction cmpt=0; cmpt<pTraits<Type2>::nComponents; cmpt++)
    {
        scalarField valuesCmpt(values.component(cmpt));
        this->setValuesFromList(sI, cellLabels, valuesCmpt);
        sI += 1;
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::setValuesFromList
(
    const direction sI,
    const labelUList& cellLabels,
    const scalarField& values
)
{
    const fvMesh& mesh = psi_.mesh();

    const cellList& cells = mesh.cells();

    Field<Type>& b = this->source();
    Field<Type>& psiIn = psi_.primitiveFieldRef();

    const direction nCmpts = pTraits<Type>::nComponents;

    boolList markCells(cells.size(), false);
    forAll(cellLabels, i)
    {
        const label cI = cellLabels[i];
        markCells[cI] = true;
    }
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    forAll(cellLabels, i)
    {
        const label cI = cellLabels[i];
        const scalar& value = values[i];

        //- set diag, source and field
        psiIn[cI](sI) = value;

        if (this->diag().activeType() != blockCoeffBase::SQUARE)
        {
            typename CoeffField<Type>::linearTypeField& blockDiag =
                this->diag().asLinear();

            b[cI](sI) = blockDiag[cI](sI)*value;

        }
        else
        {
            typename CoeffField<Type>::squareTypeField& blockDiag =
                this->diag().asSquare();

            b[cI](sI) = blockDiag[cI](sI, sI)*value;
            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
            {
                if (sI!=cmptI)
                {
                    blockDiag[cI](sI, cmptI) = 0;
                }
            }
        }

        /*
        Info<< "d: " << blockDiag[cellI] << endl;
        Info<< "b: " << b[cellI] << endl;
        Info<< "psi: " << psiIn[cellI] << endl;
        */

        //- set off diagonals
        if (this->symmetric() || this->asymmetric())
        {
            const cell& c = cells[cI];

            forAll(c, j)
            {
                const label fI = c[j];

                if (mesh.isInternalFace(fI))
                {
                    if (this->symmetric())
                    {
                        if (this->diag().activeType() != blockCoeffBase::SQUARE)
                        {
                            typename CoeffField<Type>::linearTypeField&
                                blockUpper = this->upper().asLinear();

                            blockUpper[fI](sI) = 0;
                        }
                        else
                        {
                            typename CoeffField<Type>::squareTypeField&
                                blockUpper = this->upper().asSquare();

                            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                            {
                                blockUpper[fI](sI, cmptI) = 0;
                            }
                        }
                    }
                    else
                    {
                        if (this->diag().activeType() != blockCoeffBase::SQUARE)
                        {
                            typename CoeffField<Type>::linearTypeField&
                                blockUpper = this->upper().asLinear();

                            typename CoeffField<Type>::linearTypeField&
                                blockLower = this->lower().asLinear();

                            if (markCells[own[fI]])
                            {
                                blockUpper[fI](sI) = 0;
                            }
                            if (markCells[nei[fI]])
                            {
                                blockLower[fI](sI) = 0;
                            }
                        }
                        else
                        {
                            typename CoeffField<Type>::squareTypeField&
                                blockUpper = this->upper().asSquare();

                            typename CoeffField<Type>::squareTypeField&
                                blockLower = this->lower().asSquare();

                            for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                            {
                                if (markCells[own[fI]])
                                {
                                    blockUpper[fI](sI, cmptI) = 0;
                                }
                                if (markCells[nei[fI]])
                                {
                                    blockLower[fI](sI, cmptI) = 0;
                                }
                            }
                        }
                    }
                }
                else
                {
                    const label pI = mesh.boundaryMesh().whichPatch(fI);
                    const fvPatch& patch = mesh.boundary()[pI];

                    if (patch.coupled())
                    {
                        label pFI =
                            mesh.boundaryMesh()[pI].whichFace(fI);
                        if (this->symmetric())
                        {
                            if
                            (
                                this->coupleUpper()[pI].activeType()
                                  != blockCoeffBase::SQUARE
                            )
                            {
                                typename CoeffField<Type>::linearTypeField&
                                    cUpper = this->coupleUpper()[pI].asLinear();
                                cUpper[pFI](sI) = 0;
                            }
                            else
                            {
                                typename CoeffField<Type>::squareTypeField&
                                    cUpper = this->coupleUpper()[pI].asSquare();
                                for
                                (
                                    direction cmptI = 0; cmptI < nCmpts; cmptI++
                                )
                                {
                                    cUpper[pFI](sI, cmptI) = 0;
                                }
                            }
                        }
                        else
                        {
                            if
                            (
                                this->coupleUpper()[pI].activeType()
                                  != blockCoeffBase::SQUARE
                            )
                            {
                                typename CoeffField<Type>::linearTypeField&
                                    cUpper =
                                        this->coupleUpper()[pI].asLinear();

                                typename CoeffField<Type>::linearTypeField&
                                    cLower =
                                        this->coupleLower()[pI].asLinear();
                                cUpper[pFI](sI) = 0;
                                cLower[pFI](sI) = 0;
                            }
                            else
                            {
                                typename CoeffField<Type>::squareTypeField&
                                    cUpper = this->coupleUpper()[pI].asSquare();

                                typename CoeffField<Type>::squareTypeField&
                                    cLower = this->coupleLower()[pI].asSquare();
                                for
                                (
                                    direction cmptI = 0; cmptI < nCmpts; cmptI++
                                )
                                {
                                    cUpper[pFI](sI, cmptI) = 0;
                                    cLower[pFI](sI, cmptI) = 0;
                                }
                            }
                        }
                    }
                    else
                    {
                        //- Maybe something internal and boundary coeffs
                        //  But not for momentum
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::setReference
(
    const word name,
    const label cI,
    const scalar& value
)
{
    const fvMesh& mesh = psi_.mesh();

    //- Implementation only for shifting volScalarField in closed systems

    const volScalarField& p
    (
        psi_.db().template lookupObject<volScalarField>(name)
    );

    if (p.needReference())
    {
        scalar vdiff = 0;

        direction sI = bIndex_.findStartIndex(name);
        Field<Type>& psiIn = psi_.primitiveFieldRef();

        if (Pstream::parRun())
        {
            List<scalar> procValue(Pstream::nProcs(), 0.0);
            List<bool> procLoc(Pstream::nProcs(), false);
            if (cI!=-1)
            {
                procValue[Pstream::myProcNo()] = psiIn[cI](sI);
                procLoc[Pstream::myProcNo()] = true;
            }
            Pstream::allGatherList(procValue);
            Pstream::allGatherList(procLoc);

            label procI = -1;
            bool found = false;

            do
            {
                procI++;

                if (procLoc[procI]==true)
                {
                    found = true;
                }
            }while (procI < procLoc.size() && !found);
            vdiff = procValue[procI] - value;
        }
        else
        {
            if (mesh.cells().size()>0)
            {
                vdiff = psiIn[cI](sI) - value;
            }
            else
            {
                FatalErrorInFunction
                    << "Something is wrong. Run in serial mode without cells"
                    << abort(FatalError);
            }
        }

        forAll(mesh.cells(), cI)
        {
            psiIn[cI](sI) -= vdiff;
        }
        forAll(psi_.boundaryField(), pI)
        {
            Field<Type>& pf = psi_.boundaryFieldRef()[pI];
            forAll(pf, fI)
            {
                pf[fI](sI) -= vdiff;
            }
        }
        psi_.correctBoundaryConditions();
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvBlockMatrix<Type>::
flux() const
{
    if (!psi_.mesh().schemes().fluxRequired(psi_.name()))
    {
        FatalErrorInFunction
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary"
               " of fvSchemes."
            << abort(FatalError);
    }

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tfieldFlux
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "flux("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimless,
            fvsPatchField<Type>::calculatedType(),
            psi_.boundaryField().patchTypes()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& fieldFlux =
        tfieldFlux.ref();

    fieldFlux.primitiveFieldRef() = this->faceH(psi_);

    return tfieldFlux;
}


template<class Type>
template<class Type2>
Foam::tmp<Foam::GeometricField<Type2, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvBlockMatrix<Type>::
flux
(
    const Foam::GeometricField<Type2, Foam::fvPatchField, Foam::volMesh>& field
) const
{
    if (!psi_.mesh().schemes().fluxRequired(field.name()))
    {
        FatalErrorInFunction
            << "flux requested but " << field.name()
            << " not specified in the fluxRequired sub-dictionary"
               " of fvSchemes."
            << abort(FatalError);
    }
    // A bit doggy implementation because we have double information
    // psi and the field. However we do it because we need the field Type

    direction sI = bIndex_.findStartIndex(field.name());

    tmp<GeometricField<Type2, fvsPatchField, surfaceMesh>> tfieldFlux
    (
        new GeometricField<Type2, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "flux("+field.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            field.mesh(),
            dimensioned<Type2>
            (
                "zero",
                this->dimensionSets()[sI] ,
                pTraits<Type2>::zero
            ),
            fvsPatchField<Type>::calculatedType(),
            field.boundaryField().patchTypes()
        )
    );
    GeometricField<Type2, fvsPatchField, surfaceMesh>& fieldFlux =
        tfieldFlux.ref();

    //- Internal faces
    fieldFlux.primitiveFieldRef() = this->faceH(sI, field, psi_);

    //- Boundary faces based on BC contributions
    addBoundaryFlux(fieldFlux, sI);

    //- face flux correction in laplacian
    addFaceFluxCorrection(fieldFlux, sI);

    return tfieldFlux;
}


template<class Type>
template<class fieldType>
void Foam::fvBlockMatrix<Type>::retrieveSolution
(
    const direction dir,
    Field<fieldType>& xSingle
) const
{
    const direction nCmpts = pTraits<fieldType>::nComponents;
    direction localDir = dir;

    const Field<Type> psiIn = psi_.primitiveField();

    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        xSingle.replace(cmptI, psiIn.component(localDir));

        localDir++;
    }
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::insertEquation
(
    const direction& dir,
    fvMatrix<Type2>& matrix
)
{
    if (matrix.psi().name() != this->psi().name())
    {
        insertSolutionVector(dir, matrix.psi().primitiveField());
        #if !defined( WIN32 ) && !defined( WIN64 )
        modifyBCCompSpecificMembers(dir, matrix.psi());
        #endif
    }
    insertDiagSource(dir, matrix);
    insertUpperLower(dir, matrix);
    updateCouplingCoeffs(dir, matrix);
    insertBoundaryCoeffsIfRequired(dir, matrix);
    checkAndSetDimensions(dir, pTraits<Type2>::nComponents, matrix.dimensions());
}

template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::insertEquation
(
    const direction& dir,
    tmp<fvMatrix<Type2>> matrix
)
{
    this->insertEquation(dir, matrix.ref());
}


template<class Type>
template<class Type2>
void Foam::fvBlockMatrix<Type>::insertEquation
(
    const direction& dir,
    fvBlockMatrix<Type2>& matrix
)
{
    if (matrix.psi().name() != this->psi().name())
    {
        insertSolutionVector(dir, matrix.psi().primitiveField());
    }
    insertDiagSource(dir, matrix);
    insertUpperLower(dir, matrix);
    updateCouplingCoeffs(dir, matrix);
    insertBoundaryCoeffsIfRequired(dir, matrix);
    checkAndSetDimensions(dir, matrix.dimensionSets());
}


template<class Type>
template<class blockType, class fieldType>
void Foam::fvBlockMatrix<Type>::insertBlockCoupling
(
    const direction dirI,
    const direction dirJ,
    const BlockLduSystem<blockType, fieldType>& blockSystem,
    const bool incrementColumnDir
)
{
    insertBlock(dirI, dirJ, blockSystem, incrementColumnDir);
    insertBoundaryContributions(dirI, dirJ, blockSystem, incrementColumnDir);
    insertBoundaryCoeffsIfRequired(dirI, dirJ, blockSystem, incrementColumnDir);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::insertEquationCoupling
(
    const direction dirI,
    const direction dirJ,
    const scalarField& coeffIJ
)
{
    // Multiply coefficients by volume
    scalarField coeffIJVol = coeffIJ*psi_.mesh().V();

    insertCouplingDiag(dirI, dirJ, coeffIJVol);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::insertEquationCoupling
(
    const direction dirI,
    const direction dirJ,
    const fvScalarMatrix& matrix
)
{
    insertCouplingDiagSource(dirI, dirJ, matrix);
    insertCouplingUpperLower(dirI, dirJ, matrix);
}


template<class Type>
void Foam::fvBlockMatrix<Type>::blockAdd
(
    const direction dir,
    const scalarField& xSingle,
    Field<Type>& blockX
)
{
    forAll(xSingle, cI)
    {
        blockX[cI](dir) += xSingle[cI];
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::updateSourceCoupling()
{
    // Eliminate off-diagonal block coefficients from the square diagonal
    // With this change, the segregated matrix can be assembled with complete
    // source terms and linearisation can be provided independently.
    // Once the linearisation coefficients are set (off-diagonal entries
    // in the square block matrix, they are multiplied by the current value
    // of the field and subtracted from the source term

    if (this->diag().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<Type>::squareTypeField& blockDiag =
            this->diag().asSquare();

        typename CoeffField<Type>::linearTypeField lf(blockDiag.size());
        typename CoeffField<Type>::squareTypeField sf(blockDiag.size());

        // Expand and contract

        // Take out the diagonal entries from the square coefficient
        contractLinear(lf, blockDiag);

        // Expand the diagonal for full square, with zeroes in the off-diagonal
        expandLinear(sf, lf);

        // Subtract from the source the difference between the full block
        // diagonal and the diagonal terms only
        // Sign is the same as in the derivative
        this->source() += (blockDiag - sf) & psi_.internalField();
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::addMomentumGradPCoupledBC
(
    const volScalarField& p
)
{
    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& bFields = psi_.boundaryFieldRef();
    forAll(bFields, patchi)
    {
        bFields[patchi].addMomentumGradPCoupledBC(*this, p);
    }
}



template<class Type>
Foam::BlockSolverPerformance<Type> Foam::fvBlockMatrix<Type>::solve
(
    const dictionary& solverControls
)
{
    // Solver call
    BlockSolverPerformance<Type> solverPerf =
        BlockLduSolver<Type>::New
        (
            psi_.name(),
            *this,
            solverControls
        )->solve(psi_.primitiveFieldRef(), this->source());

    solverPerf.print();

    psi_.correctBoundaryConditions();

    psi_.mesh().setBlockSolverPerformance(psi_.name(), solverPerf);

    return solverPerf;
}


template<class Type>
Foam::BlockSolverPerformance<Type> Foam::fvBlockMatrix<Type>::solve()
{
    return solve(solverDict());
}


template<class Type>
void Foam::fvBlockMatrix<Type>::relaxAndCheckDominance()
{
    word psiName = psi_.select
    (
        psi_.mesh().data::template lookupOrDefault<bool>
        ("finalIteration", false)
    );
    if (psi_.mesh().solution().relaxEquation(psiName))
    {
        this->BlockLduMatrix<Type>::relaxAndCheckDominance
        (
            psi_.primitiveField(),
            this->source(),
            psi_.mesh().solution().equationRelaxationFactor(psiName)
        );
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::relax()
{
    word psiName = psi_.select
    (
        psi_.mesh().data::template lookupOrDefault<bool>
        ("finalIteration", false)
    );
    if (psi_.mesh().solution().relaxEquation(psiName))
    {
        this->BlockLduMatrix<Type>::relax
        (
            psi_.primitiveField(),
            this->source(),
            psi_.mesh().solution().equationRelaxationFactor(psiName)
        );
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::relax(const Field<Type>& relaxField)
{
    this->BlockLduMatrix<Type>::relaxVec
    (
        psi_.primitiveField(),
        this->source(),
        relaxField
    );
}


template<class Type>
void Foam::fvBlockMatrix<Type>::relax
(
    const scalar relax,
    const labelList& cells
)
{
    this->BlockLduMatrix<Type>::relax
    (
        psi_.primitiveField(),
        this->source(),
        cells,
        relax
    );
}


template<class Type>
void Foam::fvBlockMatrix<Type>::boundaryRelax()
{
    const typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& bFields = psi_.boundaryField();
    forAll(bFields, patchi)
    {
        bFields[patchi].boundaryRelaxMatrix(*this);
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::boundaryProcRelax
(
    const scalar relax
)
{

    if ((relax <= 0) || (relax >=1))
    {
        return;
    }

    CoeffField<Type>& diagCoeff = this->diag();
    label nCoeffs = pTraits<Type>::nComponents;

    tensorField& bDiag = diagCoeff.asSquare();
    const vectorField& psiInt = psi_();

    vectorField& source = this->source();

    const scalar oneMrelax = 1-relax;
    const scalar twoMrelax = 2-relax;

    forAll(psi_.mesh().boundary(), pI)
    {
        const fvPatch& fvp =  psi_.mesh().boundary()[pI];
        if (isA<processorFvPatch>(fvp))
        {

            boolList markCells(psiInt.size(), false);

            const labelUList& faceCells = fvp.faceCells();

            forAll(faceCells, fI)
            {
                const label cI = faceCells[fI];

                if (!markCells[cI])
                {
                    for (label iC=0; iC<nCoeffs; iC++)
                    {
                        source[cI][iC] +=
                            bDiag[cI][iC*nCoeffs+iC]*psiInt[cI][iC]*oneMrelax;
                    }
                    markCells[cI] = true;
                }
            }

            markCells = false;
            forAll(faceCells, fI)
            {
                const label cI = faceCells[fI];
                if (!markCells[cI])
                {
                    for (label iC=0; iC<nCoeffs; iC++)
                    {
                        bDiag[cI][iC*nCoeffs+iC] *= twoMrelax;
                    }
                    markCells[cI] = true;
                }
            }
        }
    }
}


template<class Type>
void Foam::fvBlockMatrix<Type>::boundaryRelax
(
    const scalar relax,
    const label pI
)
{
    CoeffField<Type>& diagCoeff = this->diag();
    label nCoeffs = pTraits<Type>::nComponents;

    tensorField& bDiag = diagCoeff.asSquare();
    const vectorField& psiInt = psi_();

    vectorField& source = this->source();

    const scalar oneMrelax = 1-relax;
    const scalar twoMrelax = 2-relax;

    boolList markCells(psiInt.size(), false);
    const fvPatch& fvp =  psi_.mesh().boundary()[pI];

    const labelUList& faceCells = fvp.faceCells();

    forAll(faceCells, fI)
    {
        const label cI = faceCells[fI];

        if (!markCells[cI])
        {
            for (label iC=0; iC<nCoeffs; iC++)
            {
                source[cI][iC] +=
                    bDiag[cI][iC*nCoeffs+iC]*psiInt[cI][iC]*oneMrelax;
            }
            markCells[cI] = true;
        }
    }

    markCells = false;
    forAll(faceCells, fI)
    {
        const label cI = faceCells[fI];
        if (!markCells[cI])
        {
            for (label iC=0; iC<nCoeffs; iC++)
            {
                bDiag[cI][iC*nCoeffs+iC] *= twoMrelax;
            }
            markCells[cI] = true;
        }
    }
}


template<class Type>
const Foam::dictionary& Foam::fvBlockMatrix<Type>::solverDict() const
{
    return psi_.mesh().solution().solverDict
    (
        psi_.select
        (
            psi_.mesh().data::template lookupOrDefault<bool>
            ("finalIteration", false)
        )
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const fvBlockMatrix<Type>& bxs
)
{
    os  << static_cast<const BlockLduSystem<Type, Type>&>(bxs) << nl
        << bxs.psi() << endl;

    return os;
}


// ************************************************************************* //

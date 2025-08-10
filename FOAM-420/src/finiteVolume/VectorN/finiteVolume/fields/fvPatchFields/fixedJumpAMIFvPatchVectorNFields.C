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
    (c) 2010 Ivor Clifford
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fixedJumpAMIFvPatchVectorNFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/CoeffField/CoeffField.H"
#include "fields/fvPatchFields/derived/fixedJumpAMI/fixedJumpAMIFvPatchField.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#define makeVectorNCyclicAMI(Type, args...)                            \
                                                                               \
template<>                                                                     \
void fixedJumpAMIFvPatchField<Type>::updateInterfaceMatrix                     \
(										                                       \
    const Field<Type>& psiInternal,                                            \
    Field<Type>& result,                                                       \
    const BlockLduMatrix<Type>& matrix,                                        \
    const CoeffField<Type>& coeffs,                                            \
    const Pstream::commsTypes commsType                                        \
) const                                                                        \
{                                                                              \
    const cyclicAMIFvPatch& cAMIfvPatch = this->cyclicAMIPatch();              \
    const labelUList& nbrFaceCells =                                           \
        cAMIfvPatch.cyclicAMIPatch().nbrPatch().faceCells();                   \
                                                                               \
                                                                               \
    Field<Type> pnf(psiInternal, nbrFaceCells);                                \
    Field<Type> pif(psiInternal, cAMIfvPatch.faceCells());                     \
                                                                               \
    if (cAMIfvPatch.applyLowWeightCorrection())                                \
    {                                                                          \
        tmp<vectorField> normalFt = this->patch().nf();                        \
        const vectorField& normalF = normalFt();                               \
        Field<tensor> trans(this->size(), Zero);                               \
        forAll(trans, cI)                                                      \
        {                                                                      \
            trans[cI] = I - 2.0*sqr(normalF[cI]);                              \
        }                                                                      \
                                                                               \
        Field<vector> pifv(this->size(), vector::zero);                        \
        Field<scalar> pifs(this->size(), Zero);                                \
        forAll(pifv, cI)                                                       \
        {                                                                      \
            pifv[cI][0] = pif[cI][0];                                          \
            pifv[cI][1] = pif[cI][1];                                          \
            pifv[cI][2] = pif[cI][2];                                          \
            pifs[cI] = pif[cI][3];                                             \
        }                                                                      \
        Field<vector> mirrorFieldv(Foam::transform(trans, pifv));              \
        Field<scalar> mirrorFields(Foam::transform(trans, pifs));              \
        Field<Type> zeroField(mirrorFieldv.size(), pTraits<Type>::zero);       \
        forAll(pifv, cI)                                                       \
        {                                                                      \
            zeroField[cI][0] = mirrorFieldv[cI][0];                            \
            zeroField[cI][1] = mirrorFieldv[cI][1];                            \
            zeroField[cI][2] = mirrorFieldv[cI][2];                            \
            zeroField[cI][3] = mirrorFields[cI];                               \
        }                                                                      \
                                                                               \
        pnf = cAMIfvPatch.interpolate(pnf, zeroField);                         \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        pnf = cAMIfvPatch.interpolate(pnf);                                    \
    }                                                                          \
                                                                               \
                                                                               \
                                                                               \
    if (&psiInternal == &this->primitiveField())                               \
    {                                                                          \
        Field<Type> jf(this->jump());                                          \
        if (!this->cyclicAMIPatch().owner())                                   \
        {                                                                      \
            jf *= -1.0;                                                        \
        }                                                                      \
                                                                               \
        pnf -= jf;                                                             \
    }                                                                          \
                                                                               \
    if (transforms())                                                          \
    {                                                                          \
        transform().transform(pnf, pnf);                                       \
    }                                                                          \
                                                                               \
    if (this->cyclicAMIPatch().cyclicAMIPatch().implicitNonCovered())           \
    {                                                                          \
        this->addToInternalCoupledCoveredField(result, pnf, coeffs);           \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        this->addToInternalCoupledField(result, pnf, coeffs);                  \
    }                                                                          \
                                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
void fixedJumpAMIFvPatchField<Type>::initInterfaceMatrixUpdate                 \
(                                                                              \
    const Field<Type>& psiInternal,                                            \
    Field<Type>& result,                                                       \
    const BlockLduMatrix<Type>&,                                               \
    const CoeffField<Type>& coeffs,                                            \
    const Pstream::commsTypes commsType                                        \
) const                                                                        \
{}                                                                             \
                                                                               \
template<>                                                                     \
void fixedJumpAMIFvPatchField<Type>::addToInternalCoupledField                 \
(										                                       \
    Field<Type>& result,                                                       \
    const Field<Type>& pnf,                                                    \
    const CoeffField<Type>& coeffs                                             \
) const                                                                        \
{                                                                              \
    typename BlockCoeff<Type>::multiply mult;                                  \
    const cyclicAMIFvPatch& cAMIfvPatch = this->cyclicAMIPatch();              \
    const labelUList& fc = cAMIfvPatch.faceCells();                            \
                                                                               \
    if (cAMIfvPatch.owner())                                                   \
    {                                                                          \
        forAll(fc, cfI)                                                        \
        {                                                                      \
            label cI = fc[cfI];                                                \
                                                                               \
            if (coeffs.activeType() == blockCoeffBase::SCALAR)                 \
            {                                                                  \
                result[cI] -= mult(coeffs.asScalar()[cfI], pnf[cfI]);          \
            }                                                                  \
            else if (coeffs.activeType() == blockCoeffBase::LINEAR)            \
            {                                                                  \
                result[cI] -= mult(coeffs.asLinear()[cfI], pnf[cfI]);          \
            }                                                                  \
            else if (coeffs.activeType() == blockCoeffBase::SQUARE)            \
            {                                                                  \
                result[cI] -= mult(coeffs.asSquare()[cfI], pnf[cfI]);          \
            }                                                                  \
        }                                                                      \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        forAll(fc, cfI)                                                        \
        {                                                                      \
            label cI = fc[cfI];                                                \
                                                                               \
            if (coeffs.activeType() == blockCoeffBase::SCALAR)                 \
            {                                                                  \
                result[cI] -= mult(coeffs.asScalar()[cfI], pnf[cfI]);          \
            }                                                                  \
            else if (coeffs.activeType() == blockCoeffBase::LINEAR)            \
            {                                                                  \
                result[cI] -= mult(coeffs.asLinear()[cfI], pnf[cfI]);          \
            }                                                                  \
            else if (coeffs.activeType() == blockCoeffBase::SQUARE)            \
            {                                                                  \
                result[cI] -= mult(coeffs.asSquare()[cfI], pnf[cfI]);          \
            }                                                                  \
        }                                                                      \
    }                                                                          \
}                                                                              \
                                                                               \
template<>                                                                     \
void fixedJumpAMIFvPatchField<Type>::addToInternalCoupledCoveredField         \
(										                                       \
    Field<Type>& result,                                                       \
    const Field<Type>& pnf,                                                    \
    const CoeffField<Type>& coeffs                                             \
) const                                                                        \
{                                                                              \
    typename BlockCoeff<Type>::multiply mult;                                  \
    const cyclicAMIFvPatch& cAMIfvPatch = this->cyclicAMIPatch();              \
    const labelUList& fc = cAMIfvPatch.faceCells();                            \
                                                                               \
    if (cAMIfvPatch.owner())                                                   \
    {                                                                          \
        const labelList& srcFaces =                                            \
            cAMIfvPatch.AMI().overlapSourceFaces();                            \
                                                                               \
        forAll(srcFaces, fI)                                                   \
        {                                                                      \
            label cfI = srcFaces[fI];                                          \
            label cI = fc[cfI];                                                \
                                                                               \
            if (coeffs.activeType() == blockCoeffBase::SCALAR)                 \
            {                                                                  \
                result[cI] -= mult(coeffs.asScalar()[cfI], pnf[cfI]);          \
            }                                                                  \
            else if (coeffs.activeType() == blockCoeffBase::LINEAR)            \
            {                                                                  \
                result[cI] -= mult(coeffs.asLinear()[cfI], pnf[cfI]);          \
            }                                                                  \
            else if (coeffs.activeType() == blockCoeffBase::SQUARE)            \
            {                                                                  \
                result[cI] -= mult(coeffs.asSquare()[cfI], pnf[cfI]);          \
            }                                                                  \
        }                                                                      \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        const labelList& tgtFaces =                                            \
            cAMIfvPatch.nbrPatch().AMI().overlapTargetFaces();                 \
                                                                               \
        forAll(tgtFaces, fI)                                                   \
        {                                                                      \
            label cfI = tgtFaces[fI];                                          \
            label cI = fc[cfI];                                                \
                                                                               \
            if (coeffs.activeType() == blockCoeffBase::SCALAR)                 \
            {                                                                  \
                result[cI] -= mult(coeffs.asScalar()[cfI], pnf[cfI]);          \
            }                                                                  \
            else if (coeffs.activeType() == blockCoeffBase::LINEAR)            \
            {                                                                  \
                result[cI] -= mult(coeffs.asLinear()[cfI], pnf[cfI]);          \
            }                                                                  \
            else if (coeffs.activeType() == blockCoeffBase::SQUARE)            \
            {                                                                  \
                result[cI] -= mult(coeffs.asSquare()[cfI], pnf[cfI]);          \
            }                                                                  \
        }                                                                      \
    }                                                                          \
}

#define doMakePatchTypeField(type, Type, args...)                              \
    makeVectorNCyclicAMI(type)                                         \
                                                                               \
    makeTemplatePatchTypeField(type, fixedJumpAMI);

forAllVectorNTypes(doMakePatchTypeField)

#undef doMakePatchTypeField

#undef makeVectorNCyclicAMI

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

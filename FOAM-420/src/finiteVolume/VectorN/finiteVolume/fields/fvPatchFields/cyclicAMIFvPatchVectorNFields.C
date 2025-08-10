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
    (c) 2020-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cyclicAMIFvPatchVectorNFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/CoeffField/CoeffField.H"
#include "fields/fvPatchFields/constraint/cyclicAMI/cyclicAMIFvPatchField.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#define makeVectorNCyclicAMI(Type, args...)                                    \
template<>                                                                     \
Foam::tmp<Foam::Field<Type>>                                                   \
Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField                         \
(                                                                              \
    const Pstream::commsTypes                                                  \
) const                                                                        \
{                                                                              \
    const Field<Type>& iField = this->primitiveField();                        \
    const labelUList& nbrFaceCells =                                           \
        cyclicAMIPatch_.cyclicAMIPatch().nbrPatch().faceCells();               \
                                                                               \
    Field<Type> pnf(iField, nbrFaceCells);                                     \
                                                                               \
    tmp<Field<Type>> tpnf;                                                     \
                                                                               \
    if (cyclicAMIPatch_.applyLowWeightCorrection())                            \
    {                                                                          \
        tmp<vectorField> normalFt = this->patch().nf();                        \
        const vectorField& normalF = normalFt();                               \
        Field<tensor> trans(this->size(), Zero);                               \
        forAll(trans, cI)                                                      \
        {                                                                      \
            trans[cI] = I - 2.0*sqr(normalF[cI]);                              \
        }                                                                      \
                                                                               \
        Field<Type> mirrorField( this->patchInternalField()());                \
                                                                               \
        /*  Hardcoded for v4 Up now. needs refactoring         */              \
        if (pTraits<Type>::nComponents == 4)                                   \
        {                                                                      \
            forAll(mirrorField, fI)                                            \
            {                                                                  \
                vector v1(Zero);                                               \
                v1[0] = mirrorField[fI][0];                                    \
                v1[1] = mirrorField[fI][1];                                    \
                v1[2] = mirrorField[fI][2];                                    \
                v1 = Foam::transform(trans[fI], v1);                           \
                mirrorField[fI][0] = v1[0];                                    \
                mirrorField[fI][1] = v1[1];                                    \
                mirrorField[fI][2] = v1[2];                                    \
                mirrorField[fI][3] =                                           \
                    Foam::transform(trans[fI], mirrorField[fI][3]);            \
            }                                                                  \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            forAll(mirrorField, fI)                                            \
            {                                                                  \
                for                                                            \
                (                                                              \
                    direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++\
                )                                                              \
                {                                                              \
                    mirrorField[fI][cmpt] =                                    \
                        Foam::transform(trans[fI], mirrorField[fI][cmpt]);     \
                }                                                              \
            }                                                                  \
        }                                                                      \
                                                                               \
        tpnf = cyclicAMIPatch_.interpolate(pnf, mirrorField);                  \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        tpnf = cyclicAMIPatch_.interpolate(pnf);                               \
    }                                                                          \
                                                                               \
    if (transforms())                                                          \
    {                                                                          \
        transform().transform(tpnf.ref(), tpnf());                             \
    }                                                                          \
                                                                               \
    return tpnf;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::Field<Type>>                                                  \
Foam::cyclicAMIFvPatchField<Type>::patchInternalField() const                  \
{                                                                              \
    const Field<Type>& iField = this->internalField();                         \
    const labelUList& faceCells =                                              \
        cyclicAMIPatch_.cyclicAMIPatch().faceCells();                          \
                                                                               \
    tmp<Field<Type>> tpnf(new Field<Type>(iField, faceCells));                \
    Field<Type>& pnf = tpnf.ref();                                             \
                                                                               \
    if (cyclicAMIPatch_.applyLowWeightCorrection())                            \
    {                                                                          \
        const scalarField& weights = cyclicAMIPatch_.cyclicAMIPatch().         \
            lowCorrAMIWeights();                                               \
        tmp<vectorField> normalFt = this->patch().nf();                        \
        const vectorField& normalF = normalFt();                               \
        Field<tensor> trans(this->size(), Zero);                               \
        forAll(trans, cI)                                                      \
        {                                                                      \
            trans[cI] = I - 2.0*sqr(normalF[cI]);                              \
        }                                                                      \
                                                                               \
        Field<Type> mirrorField(pnf.size(), Zero);                             \
                                                                               \
        /*  Hardcoded for v4 Up now. needs refactoring         */              \
        if (pTraits<Type>::nComponents == 4)                                   \
        {                                                                      \
            forAll(pnf, fI)                                                    \
            {                                                                  \
                vector v1(Zero);                                               \
                v1[0] = pnf[fI][0];                                            \
                v1[1] = pnf[fI][1];                                            \
                v1[2] = pnf[fI][2];                                            \
                v1 = Foam::transform(trans[fI], v1);                           \
                mirrorField[fI][0] = v1[0];                                    \
                mirrorField[fI][1] = v1[1];                                    \
                mirrorField[fI][2] = v1[2];                                    \
                mirrorField[fI][3] = pnf[fI][3];                               \
            }                                                                  \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            forAll(mirrorField, fI)                                            \
            {                                                                  \
                for                                                            \
                (                                                              \
                    direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++\
                )                                                              \
                {                                                              \
                    mirrorField[fI][cmpt] = pnf[fI][cmpt];                     \
                }                                                              \
            }                                                                  \
        }                                                                      \
                                                                               \
                                                                               \
        pnf = (1-weights)*mirrorField + weights*pnf;                           \
    }                                                                          \
                                                                               \
    return tpnf;                                                               \
}                                                                              \
                                                                               \
template<>                                                                     \
void cyclicAMIFvPatchField<Type>::updateInterfaceMatrix                        \
(                                                                              \
    const Field<Type>& psiInternal,                                            \
    Field<Type>& result,                                                       \
    const BlockLduMatrix<Type>& matrix,                                        \
    const CoeffField<Type>& coeffs,                                            \
    const Pstream::commsTypes commsType                                        \
) const                                                                        \
{                                                                              \
    const labelUList& nbrFaceCells =                                           \
        cyclicAMIPatch_.cyclicAMIPatch().nbrPatch().faceCells();               \
                                                                               \
                                                                               \
    Field<Type> pnf(psiInternal, nbrFaceCells);                                \
    Field<Type> pif(psiInternal, cyclicAMIPatch_.faceCells());                 \
                                                                               \
    if (cyclicAMIPatch_.applyLowWeightCorrection())                            \
    {                                                                          \
        tmp<vectorField> normalFt = this->patch().nf();                        \
        const vectorField& normalF = normalFt();                               \
        Field<tensor> trans(this->size(), Zero);                               \
        forAll(trans, cI)                                                      \
        {                                                                      \
            trans[cI] = I - 2.0*sqr(normalF[cI]);                              \
        }                                                                      \
        Field<Type> zeroField(this->size(), pTraits<Type>::zero);              \
                                                                               \
        /*  Hardcoded for v4 Up now. needs refactoring         */              \
        if (pTraits<Type>::nComponents == 4)                                   \
        {                                                                      \
            Field<vector> pifv(this->size(), vector::zero);                    \
            Field<scalar> pifs(this->size(), Zero);                            \
            forAll(pifv, cI)                                                   \
            {                                                                  \
                pifv[cI][0] = pif[cI][0];                                      \
                pifv[cI][1] = pif[cI][1];                                      \
                pifv[cI][2] = pif[cI][2];                                      \
                pifs[cI] = pif[cI][3];                                         \
            }                                                                  \
            Field<vector> mirrorFieldv(Foam::transform(trans, pifv));          \
            Field<scalar> mirrorFields(Foam::transform(trans, pifs));          \
                                                                               \
            forAll(pifv, cI)                                                   \
            {                                                                  \
                zeroField[cI][0] = mirrorFieldv[cI][0];                        \
                zeroField[cI][1] = mirrorFieldv[cI][1];                        \
                zeroField[cI][2] = mirrorFieldv[cI][2];                        \
                zeroField[cI][3] = mirrorFields[cI];                           \
            }                                                                  \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            forAll(zeroField, fI)                                              \
            {                                                                  \
                for                                                            \
                (                                                              \
                    direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++\
                )                                                              \
                {                                                              \
                    zeroField[fI][cmpt] =                                      \
                        Foam::transform(trans[fI], zeroField[fI][cmpt]);       \
                }                                                              \
            }                                                                  \
        }                                                                      \
                                                                               \
        pnf = cyclicAMIPatch_.interpolate(pnf, zeroField);                     \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        pnf = cyclicAMIPatch_.interpolate(pnf);                                \
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
void cyclicAMIFvPatchField<Type>::initInterfaceMatrixUpdate                    \
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
void cyclicAMIFvPatchField<Type>::addToInternalCoupledField                    \
(                                                                              \
    Field<Type>& result,                                                       \
    const Field<Type>& pnf,                                                    \
    const CoeffField<Type>& coeffs                                             \
) const                                                                        \
{                                                                              \
    typename BlockCoeff<Type>::multiply mult;                                  \
    const labelUList& fc = cyclicAMIPatch_.faceCells();                        \
                                                                               \
    if (cyclicAMIPatch_.owner())                                               \
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
void cyclicAMIFvPatchField<Type>::addToInternalCoupledCoveredField             \
(                                                                              \
    Field<Type>& result,                                                       \
    const Field<Type>& pnf,                                                    \
    const CoeffField<Type>& coeffs                                             \
) const                                                                        \
{                                                                              \
    typename BlockCoeff<Type>::multiply mult;                                  \
    const labelUList& fc = cyclicAMIPatch_.faceCells();                        \
                                                                               \
    if (cyclicAMIPatch_.owner())                                               \
    {                                                                          \
        const labelList& srcFaces =                                            \
            cyclicAMIPatch_.AMI().overlapSourceFaces();                        \
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
            cyclicAMIPatch_.nbrPatch().AMI().overlapTargetFaces();             \
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
    makeTemplatePatchTypeField(type, cyclicAMI);                               \
    makeVectorNCyclicAMI(type);

forAllVectorNTypes(doMakePatchTypeField)

#define doMakePatchFieldTypeName(type, Type, args...)                          \
    makePatchFieldTypeName(type, cyclicAMI);

forAllTensorNTypes(doMakePatchFieldTypeName)
forAllDiagTensorNTypes(doMakePatchFieldTypeName)
forAllSphericalTensorNTypes(doMakePatchFieldTypeName)

#undef doMakePatchTypeField
#undef doMakePatchFieldTypeName

#undef makeVectorNCyclicAMI

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

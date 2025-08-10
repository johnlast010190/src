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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/solutionControl/foamCoupledControl/foamCoupledControl.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::implicitCoeffs() const
{
    bool doImplicitCoeffs = false;

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    if (mesh.foundObject<foamCoupledControl>(solutionControl::typeName))
    {
        Switch implicitAMI
        (
            debug::optimisationSwitches().lookupOrDefault<Switch>
            ("implicitNonCoveredAMI", true)
        );
        const word& fieldName = this->internalField().name();

        if
        (
            implicitAMI &&
            (
                (fieldName == ("U")) || fieldName == ("p")
            )
        )
        {
            doImplicitCoeffs = true;
        }
    }
    return doImplicitCoeffs;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    calculatedMode_(true)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    calculatedMode_(true)
{
    if (!isA<cyclicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value") && this->coupled())
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    calculatedMode_(true)
{
    if (!isA<cyclicAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_),
    calculatedMode_(true)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_),
    calculatedMode_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::coupled() const
{
    return cyclicAMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    if (calculatedMode_)
    {
        // Because the underlying patch is coupled, the boundary values were
        // obtained by linear interpolation from internal and neighbour values
        // (in the default implementation of evaluate()). We use this to deduce
        // the neighbour field values. This is necessary because a derived type
        // may have changed the implemention of patchNeighbourField below,
        // but will be substituted by this constraint type in field calculations
        return
            (*this - this->patch().weights()*this->patchInternalField())
           /(1-this->patch().weights());
    }
    else
    {
        const Field<Type>& iField = this->primitiveField();
        const labelUList& nbrFaceCells =
            cyclicAMIPatch_.cyclicAMIPatch().nbrPatch().faceCells();

        Field<Type> pnf(iField, nbrFaceCells);

        tmp<Field<Type>> tpnf;

        if (cyclicAMIPatch_.applyLowWeightCorrection())
        {
            tmp<vectorField> normalFt = this->patch().nf();
            const vectorField& normalF = normalFt();
            Field<tensor> trans(this->size(), Zero);
            forAll(trans, cI)
            {
                trans[cI] = I - 2.0*sqr(normalF[cI]);
            }
            Field<Type> mirrorField
            (
                Foam::transform(trans, this->patchInternalField()())
            );

            tpnf = cyclicAMIPatch_.interpolate(pnf, mirrorField);
        }
        else
        {
            tpnf = cyclicAMIPatch_.interpolate(pnf);
        }

        transform().transform(tpnf.ref(), tpnf());

        return tpnf;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::patchInternalField() const
{
    const Field<Type>& iField = this->internalField();
    const labelUList& faceCells =
        cyclicAMIPatch_.cyclicAMIPatch().faceCells();

    tmp<Field<Type>> tpnf(new Field<Type>(iField, faceCells));
    Field<Type>& pnf = tpnf.ref();

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        const scalarField& weights = cyclicAMIPatch_.cyclicAMIPatch().
            lowCorrAMIWeights();
        tmp<vectorField> normalFt = this->patch().nf();
        const vectorField& normalF = normalFt();
        Field<tensor> trans(this->size(), Zero);
        forAll(trans, cI)
        {
            trans[cI] = I - 2.0*sqr(normalF[cI]);
        }
        Field<Type> mirrorField(Foam::transform(trans, pnf));

        pnf = (1-weights)*mirrorField + weights*pnf;
    }

    return tpnf;
}


template<class Type>
const Foam::cyclicAMIFvPatchField<Type>&
Foam::cyclicAMIFvPatchField<Type>::nbrPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const cyclicAMIFvPatchField<Type>>
    (
        fld.boundaryField()[cyclicAMIPatch_.nbrPatchID()]
    );
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    calculatedMode_ = false;
    coupledFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().nbrPatch().faceCells();

    scalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // patch Internal field
    Field<scalar> pf(psiInternal,cyclicAMIPatch_.faceCells());

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        scalarField pif(psiInternal, cyclicAMIPatch_.faceCells());
        tmp<vectorField> normalFt = this->patch().nf();
        const vectorField& normalF = normalFt();
        Field<tensor> trans(this->size(), Zero);
        forAll(trans, cI)
        {
            trans[cI] = I - 2.0*sqr(normalF[cI]);
        }
        Field<scalar> mirrorField(Foam::transform(trans, pf));

        pnf = cyclicAMIPatch_.interpolate(pnf, mirrorField);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().nbrPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        Field<Type> pif(psiInternal, cyclicAMIPatch_.faceCells());
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    tmp<Field<Type>> tvIC =
        this->coupledFvPatchField<Type>::valueInternalCoeffs(w);
    Field<Type>& vIC = tvIC.ref();

    if
    (
        cyclicAMIPatch_.applyLowWeightCorrection()
      && this->cyclicAMIPatch().cyclicAMIPatch().implicitNonCovered()
    )
    {
        if (implicitCoeffs())
        {
            //- Mimic implicit zero Neumann or Dirichlet
            // with cyclicAMI at non-overlaping faces
            // A bit hacky implementation,
            // Implementation needs revision
            // Only for coupled solver currently
            Type defValue = pTraits<Type>::zero;
            if (this->internalField().name() == "p")
            {
                defValue = pTraits<Type>::one;
            }

            if (cyclicAMIPatch_.cyclicAMIPatch().owner())
            {
                const scalarList& srcWSum =
                    cyclicAMIPatch_.AMI().srcWeightsSum();
                forAll(vIC, fI)
                {
                    scalar w2 = 1-srcWSum[fI];
                    vIC[fI] = srcWSum[fI]*vIC[fI]+w2*defValue;
                }
            }
            else
            {
                const scalarList& tgtWSum =
                    cyclicAMIPatch_.nbrPatch().AMI().tgtWeightsSum();

                forAll(vIC, fI)
                {
                    scalar w2 = 1-tgtWSum[fI];
                    vIC[fI] = tgtWSum[fI]*vIC[fI]+w2*defValue;
                }
            }
        }
    }
    return tvIC;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    tmp<Field<Type>> tvBC =
        this->coupledFvPatchField<Type>::valueBoundaryCoeffs(w);
    Field<Type>& vBC = tvBC.ref();

    if
    (
        cyclicAMIPatch_.applyLowWeightCorrection()
      && this->cyclicAMIPatch().cyclicAMIPatch().implicitNonCovered()
    )
    {
        if (implicitCoeffs())
        {
            Type defValue = pTraits<Type>::zero;
            if (cyclicAMIPatch_.cyclicAMIPatch().owner())
            {
                const labelList& nsrcFaces =
                    cyclicAMIPatch_.AMI().nonOverlapSourceFaces();

                forAll(nsrcFaces, fI)
                {
                    vBC[nsrcFaces[fI]] = defValue;
                }
            }
            else
            {
                const labelList& ntgtFaces = cyclicAMIPatch_.
                    nbrPatch().AMI().nonOverlapTargetFaces();

                forAll(ntgtFaces, fI)
                {
                    vBC[ntgtFaces[fI]] = defValue;
                }
            }
        }
    }
    return tvBC;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    tmp<Field<Type>> tvIC =
        this->coupledFvPatchField<Type>::gradientInternalCoeffs(deltaCoeffs);
    Field<Type>& vIC = tvIC.ref();

    if
    (
        cyclicAMIPatch_.applyLowWeightCorrection()
      && this->cyclicAMIPatch().cyclicAMIPatch().implicitNonCovered()
    )
    {
        if (implicitCoeffs())
        {
            Type defValue = pTraits<Type>::zero;
            if (cyclicAMIPatch_.cyclicAMIPatch().owner())
            {
                const labelList& nsrcFaces =
                    cyclicAMIPatch_.AMI().nonOverlapSourceFaces();
                forAll(nsrcFaces, fI)
                {
                    vIC[nsrcFaces[fI]] = defValue;
                }
            }
            else
            {
                const labelList& ntgtFaces =
                    cyclicAMIPatch_.nbrPatch().AMI().
                    nonOverlapTargetFaces();
                forAll(ntgtFaces, fI)
                {
                    vIC[ntgtFaces[fI]] = defValue;
                }
            }
        }
    }
    return tvIC;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    tmp<Field<Type>> tvBC =
        this->coupledFvPatchField<Type>::gradientBoundaryCoeffs(deltaCoeffs);
    Field<Type>& vBC = tvBC.ref();

    if
    (
        cyclicAMIPatch_.applyLowWeightCorrection()
      && this->cyclicAMIPatch().cyclicAMIPatch().implicitNonCovered()
    )
    {
        if (implicitCoeffs())
        {
            Type defValue = pTraits<Type>::zero;
            if (cyclicAMIPatch_.cyclicAMIPatch().owner())
            {
                const labelList& nsrcFaces =
                    cyclicAMIPatch_.AMI().nonOverlapSourceFaces();
                forAll(nsrcFaces, fI)
                {
                    vBC[nsrcFaces[fI]] = defValue;
                }
            }
            else
            {
                const labelList& ntgtFaces =
                    cyclicAMIPatch_.nbrPatch().AMI().
                    nonOverlapTargetFaces();
                forAll(ntgtFaces, fI)
                {
                    vBC[ntgtFaces[fI]] = defValue;
                }
            }
        }
    }
    return tvBC;
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //

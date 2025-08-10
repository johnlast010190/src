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
    (c) 2015 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/IOobject/IOobject.H"
#include "db/dictionary/dictionary.H"
#include "fvMesh/fvMesh.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/directFvPatchFieldMapper.H"
#include "volMesh/volMesh.H"
#include "meshes/polyMesh/zones/faceZone/faceZone.H"
#include "meshes/polyMesh/polyPatches/directPolyPatch/directPolyPatch.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"
#include "VectorN/primitives/VectorN.H"
#include "VectorN/primitives/TensorN.H"

// * * * * * * * * * * * * * Static member function  * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::fvPatchField<Type>::calculatedType()
{
    return calculatedFvPatchField<Type>::typeName;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    Field<Type>(p.size(), pTraits<Type>::zero),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(word::null),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(pTraits<Type>::zero)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Type& value
)
:
    Field<Type>(p.size(), value),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(word::null),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(pTraits<Type>::zero)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const word& patchType
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(patchType),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(pTraits<Type>::zero)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(word::null),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(pTraits<Type>::zero)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(dict.lookupOrDefault<word>("patchType", word::null)),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_
    (
        dict.lookupOrDefault<Type>
        (
            "defaultGIBValue",
            pTraits<Type>::zero
        )
    )
{
    if (valueRequired)
    {
        if (dict.found("value"))
        {
            Field<Type>::operator=
            (
                Field<Type>("value", dict, p.size())
            );
        }
        else
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Essential entry 'value' missing"
                << exit(FatalIOError);
        }
    }
    else
    {
        Field<Type>::operator=(pTraits<Type>::zero);
    }
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    Field<Type>(p.size(), Zero),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(ptf.defaultGIBValue_)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        Field<Type>& f = *this;

        if (isA<directFvPatchFieldMapper>(mapper))
        {
            const directFvPatchFieldMapper& dm =
                dynamic_cast<const directFvPatchFieldMapper&>(mapper);

            if (dm.addressing().size())
            {
                Field<Type> pif(this->patchInternalField());

                const labelList& mapAddressing = dm.addressing();

                forAll(mapAddressing, i)
                {
                    if (mapAddressing[i] < 0)
                    {
                        f[i] = pif[i];
                    }
                }
            }
        }
        else if (isA<generalFvPatchFieldMapper>(mapper) && !mapper.indirect())
        {
            const generalFvPatchFieldMapper& gm =
                dynamic_cast<const generalFvPatchFieldMapper&>(mapper);

            if (gm.direct() && gm.directAddressing().size())
            {
                // Direct addressing
                Field<Type> pif(this->patchInternalField());

                const labelList& mapAddressing = gm.directAddressing();

                forAll(mapAddressing, i)
                {
                    if (mapAddressing[i] < 0)
                    {
                        f[i] = pif[i];
                    }
                }
            }
            else if (!gm.direct() && gm.addressing().size())
            {
                // Weighted addressing
                Field<Type> pif(this->patchInternalField());

                const labelListList& mapAddressing = gm.addressing();

                forAll(mapAddressing, i)
                {
                    const labelList& localAddrs = mapAddressing[i];

                    if (!localAddrs.size())
                    {
                        f[i] = pif[i];
                    }
                }
            }
        }
    }

    mapper(*this, ptf);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(ptf.defaultGIBValue_)
{}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_),
    oldTimeFieldPtr_(nullptr),
    defaultGIBValue_(ptf.defaultGIBValue_)
{

    const polyPatch& pthis = this->patch().patch();
    if (isA<indirectPolyPatch>(pthis))
    {
        if (ptf.oldTimeFieldPtr_.valid())
        {
            oldTimeFieldPtr_.reset(new Field<Type>(ptf.oldTimeFieldPtr_()));
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::objectRegistry& Foam::fvPatchField<Type>::db() const
{
    if (notNull(internalField_))
    {
        return internalField_.db();
    }
    else
    {
        return patch_.boundaryMesh().mesh();
    }
}


template<class Type>
void Foam::fvPatchField<Type>::check(const fvPatchField<Type>& ptf) const
{
    if (&patch_ != &(ptf.patch_))
    {
        FatalErrorInFunction
            << "different patches for fvPatchField<Type>s"
            << abort(FatalError);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvPatchField<Type>::snGrad() const
{
    return patch_.deltaCoeffs()*(*this - patchInternalField());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fvPatchField<Type>::patchInternalField() const
{
    return patch_.patchInternalField(internalField_);
}


template<class Type>
void Foam::fvPatchField<Type>::patchInternalField(Field<Type>& pif) const
{
    patch_.patchInternalField(internalField_, pif);
}


template<class Type>
void Foam::fvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mapper(*this, *this);
}


template<class Type>
void Foam::fvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::fvPatchField<Type>::updateGIB()
{
    const polyPatch& pthis = this->patch().patch();
    if (isA<indirectPolyPatch>(pthis))
    {
        storeGIB();
        const indirectPolyPatch& dpp =
        refCast<const indirectPolyPatch>(pthis);
        this->setSize(dpp.size());
        fvPatchField<Type>::forceAssign(pTraits<Type>::zero);
    }
}


template<class Type>
void Foam::fvPatchField<Type>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{}


template<class Type>
void Foam::fvPatchField<Type>::storeGIB()
{
    const polyPatch& pthis = this->patch().patch();
    if (isA<indirectPolyPatch>(pthis))
    {
        oldTimeFieldPtr_.reset(new Field<Type>(*this));
    }
}


template<class Type>
void Foam::fvPatchField<Type>::updateCoeffs()
{
    updated_ = true;
}


template<class Type>
void Foam::fvPatchField<Type>::resetUpdate()
{
    updated_ = false;
}


template<class Type>
Foam::tmp<Foam::CoeffField<Type>>
Foam::fvPatchField<Type>::gradientInternalBCoeffs() const
{
    //- Return segregated Coeffs in block structure
    typedef CoeffField<Type> TypeCoeffField;
    typedef typename TypeCoeffField::linearTypeField
        linearTypeField;

    tmp<TypeCoeffField> bct
    (
        new TypeCoeffField(this->size())
    );
    TypeCoeffField& bc = bct.ref();

    linearTypeField& bcLinear = bc.asLinear();

    tmp<Field<Type>> gbc = gradientInternalCoeffs();

    bcLinear = gbc();

    return bct;
}


template<class Type>
Foam::tmp<Foam::CoeffField<Type>>
Foam::fvPatchField<Type>::gradientInternalBCoeffs
(
    const scalarField& dCoeffs
) const
{
    //- Return segregated Coeffs in block structure
    typedef CoeffField<Type> TypeCoeffField;
    typedef typename TypeCoeffField::linearTypeField
        linearTypeField;

    tmp<TypeCoeffField> bct
    (
        new TypeCoeffField(this->size())
    );
    TypeCoeffField& bc = bct.ref();

    linearTypeField& bcLinear = bc.asLinear();

    tmp<Field<Type>> gbc = gradientInternalCoeffs(dCoeffs);

    bcLinear = gbc();

    return bct;
}


template<class Type>
void Foam::fvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated_)
    {
        updateCoeffs();
    }

    updated_ = false;
    manipulatedMatrix_ = false;
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix(fvMatrix<Type>& matrix)
{
    manipulatedMatrix_ = true;
}


template<class Type>
void Foam::fvPatchField<Type>::boundaryRelaxMatrix
(
    fvBlockMatrix<Type>& bEq
) const
{}


template<class Type>
void Foam::fvPatchField<Type>::addMomentumGradPCoupledBC
(
    fvBlockMatrix<Type>& bEq,
    const volScalarField& p
)
{}


template<class Type>
void Foam::fvPatchField<Type>::addContinuityCoupledBC
(
    BlockLduSystem<Foam::vector, Foam::scalar>& bEq,
    const volScalarField& p,
    const surfaceTensorField& rDf
) const
{}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fvPatchField<Type>::boundarySources() const
{
    Field<Type> df(this->size(), pTraits<Type>::zero);
    return boundarySources(*this, *this, df);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fvPatchField<Type>::boundarySources
(
    const Field<Type>& pf,
    Field<Type>& df
) const
{
    return boundarySources(*this, pf, df);
}


// Partial specialisation for VectorNs (see fvPatchFields.C for general version)

template<class Type>
template<int N>
Foam::tmp<Foam::Field<Foam::VectorN<Foam::scalar,N>>>
Foam::fvPatchField<Type>::boundarySources
(
    const fvPatchField<VectorN<scalar,N>>& pf,
    const Field<VectorN<scalar,N>>& f,
    Field<VectorN<scalar,N>>& df
)
{
    NotImplemented;
}

template<class Type>
template<int N>
Foam::tmp<Foam::Field<Foam::TensorN<Foam::scalar,N>>>
Foam::fvPatchField<Type>::boundarySources
(
    const fvPatchField<TensorN<scalar,N>>& pf,
    const Field<TensorN<scalar,N>>& f,
    Field<TensorN<scalar,N>>& df
)
{
    NotImplemented;
}

template<class Type>
template<int N>
Foam::tmp<Foam::Field<Foam::DiagTensorN<Foam::scalar,N>>>
Foam::fvPatchField<Type>::boundarySources
(
    const fvPatchField<DiagTensorN<scalar,N>>& pf,
    const Field<DiagTensorN<scalar,N>>& f,
    Field<DiagTensorN<scalar,N>>& df
)
{
    NotImplemented;
}

template<class Type>
template<int N>
Foam::tmp<Foam::Field<Foam::SphericalTensorN<Foam::scalar,N>>>
Foam::fvPatchField<Type>::boundarySources
(
    const fvPatchField<SphericalTensorN<scalar,N>>& pf,
    const Field<SphericalTensorN<scalar,N>>& f,
    Field<SphericalTensorN<scalar,N>>& df
)
{
    NotImplemented;
}


template<class Type>
void Foam::fvPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());

    if (overridesConstraint())
    {
        os.writeEntry("patchType", patch().type());
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorInFunction
            << "incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


template<class Type>
void Foam::fvPatchField<Type>::forceAssign
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::forceAssign
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::forceAssign
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/fvPatchFields/fvPatchField/fvPatchFieldNew.C"

// ************************************************************************* //

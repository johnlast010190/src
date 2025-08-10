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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/IOobject/IOobject.H"
#include "db/dictionary/dictionary.H"
#include "faMesh/faMesh.H"
#include "fields/faPatchFields/faPatchField/faPatchFieldMapper.H"
#include "areaMesh/areaMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
faPatchField<Type>::faPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    Field<Type>(p.size(), pTraits<Type>::zero),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(word::null)
{}

template<class Type>
Foam::faPatchField<Type>::faPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const word& patchType
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(patchType)
{}


template<class Type>
faPatchField<Type>::faPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_()
{}


template<class Type>
faPatchField<Type>::faPatchField
(
    const faPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    Field<Type>(mapper(ptf)),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        Field<Type>& f = *this;

        const generalFaPatchFieldMapper& gm =
            dynamic_cast<const generalFaPatchFieldMapper&>(mapper);

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

    mapper(*this, ptf);
}


template<class Type>
faPatchField<Type>::faPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(dict.lookupOrDefault<word>("patchType", word::null))
{
    if (dict.found("value"))
    {
        faPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else if (!valueRequired)
    {
        faPatchField<Type>::operator=(pTraits<Type>::zero);
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


template<class Type>
faPatchField<Type>::faPatchField
(
    const faPatchField<Type>& ptf
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_)
{}


template<class Type>
faPatchField<Type>::faPatchField
(
    const faPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const objectRegistry& faPatchField<Type>::db() const
{
    return patch_.boundaryMesh().mesh().mesh();
}


template<class Type>
void faPatchField<Type>::check(const faPatchField<Type>& ptf) const
{
    if (&patch_ != &(ptf.patch_))
    {
        FatalErrorInFunction
            << "different patches for faPatchField<Type>s"
            << abort(FatalError);
    }
}


// Return gradient at boundary
template<class Type>
tmp<Field<Type>> faPatchField<Type>::snGrad() const
{
    return patch_.deltaCoeffs()*(*this - patchInternalField());
}

// Return internal field next to patch as patch field
template<class Type>
tmp<Field<Type>> faPatchField<Type>::patchInternalField() const
{
    return patch_.patchInternalField(internalField_);
}

template<class Type>
void Foam::faPatchField<Type>::patchInternalField(Field<Type>& pif) const
{
    patch_.patchInternalField(internalField_, pif);
}

template<class Type>
void faPatchField<Type>::autoMap
(
    const faPatchFieldMapper& mapper
)
{
    mapper(*this, *this);
}


template<class Type>
void faPatchField<Type>::rmap
(
    const faPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}

template<class Type>
void Foam::faPatchField<Type>::updateCoeffs()
{
    updated_ = true;
}


template<class Type>
void Foam::faPatchField<Type>::updateCoeffs(const scalarField& weights)
{
    if (!updated_)
    {
        updateCoeffs();

        Field<Type>& fld = *this;
        fld *= weights;

        updated_ = true;
    }
}


template<class Type>
void faPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated_)
    {
        updateCoeffs();
    }

    updated_ = false;
    manipulatedMatrix_ = false;
}

template<class Type>
void Foam::faPatchField<Type>::manipulateMatrix(faMatrix<Type>& matrix)
{
    manipulatedMatrix_ = true;
}


template<class Type>
void Foam::faPatchField<Type>::manipulateMatrix
(
    faMatrix<Type>& matrix,
    const scalarField& weights
)
{
    manipulatedMatrix_ = true;
}


template<class Type>
void faPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());
    if (patchType_.size())
    {
        os.writeEntry("patchType", patchType_);
    }
}

template<class Type>
template<class EntryType>
void Foam::faPatchField<Type>::writeEntryIfDifferent
(
    Ostream& os,
    const word& entryName,
    const EntryType& value1,
    const EntryType& value2
) const
{
    if (value1 != value2)
    {
        os.writeEntry(entryName, value2);
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void faPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void faPatchField<Type>::operator=
(
    const faPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void faPatchField<Type>::operator+=
(
    const faPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void faPatchField<Type>::operator-=
(
    const faPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void faPatchField<Type>::operator*=
(
    const faPatchField<scalar>& ptf
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
void faPatchField<Type>::operator/=
(
    const faPatchField<scalar>& ptf
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
void faPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void faPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void faPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void faPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void faPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void faPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void faPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void faPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void faPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


// Force an assignment, overriding fixedValue status
template<class Type>
void faPatchField<Type>::forceAssign
(
    const faPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void faPatchField<Type>::forceAssign
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void faPatchField<Type>::forceAssign
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const faPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check("Ostream& operator<<(Ostream&, const faPatchField<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/faPatchFields/faPatchField/faPatchFieldNew.C"

// ************************************************************************* //

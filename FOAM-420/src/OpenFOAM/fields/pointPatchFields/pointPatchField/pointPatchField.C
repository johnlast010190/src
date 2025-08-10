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
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "meshes/pointMesh/pointMesh.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    patch_(p),
    internalField_(iF),
    updated_(false),
    patchType_(word::null)
{}


template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    patch_(p),
    internalField_(iF),
    updated_(false),
    patchType_(dict.lookupOrDefault<word>("patchType", word::null))
{}


template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper&
)
:
    patch_(p),
    internalField_(iF),
    updated_(false),
    patchType_(ptf.patchType_)
{}


template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatchField<Type>& ptf
)
:
    patch_(ptf.patch_),
    internalField_(ptf.internalField_),
    updated_(false),
    patchType_(ptf.patchType_)
{}


template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    patch_(ptf.patch_),
    internalField_(iF),
    updated_(false),
    patchType_(ptf.patchType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::objectRegistry& Foam::pointPatchField<Type>::db() const
{
    return patch_.boundaryMesh().mesh()();
}


template<class Type>
void Foam::pointPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());

    if (patchType_.size())
    {
        os.writeEntry("patchType", patchType_);
    }
}


template<class Type>
template<class EntryType>
void Foam::pointPatchField<Type>::writeEntryIfDifferent
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


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::pointPatchField<Type>::patchInternalField() const
{
    return patchInternalField(primitiveField());
}


template<class Type>
template<class Type1>
Foam::tmp<Foam::Field<Type1>>
Foam::pointPatchField<Type>::patchInternalField
(
    const Field<Type1>& iF,
    const labelList& meshPoints
) const
{
    // Check size
    if (iF.size() != primitiveField().size())
    {
        FatalErrorInFunction
            << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << primitiveField().size()
            << abort(FatalError);
    }

    return tmp<Field<Type1>>(new Field<Type1>(iF, meshPoints));
}


template<class Type>
template<class Type1>
Foam::tmp<Foam::Field<Type1>>
Foam::pointPatchField<Type>::patchInternalField
(
    const Field<Type1>& iF
) const
{
    return patchInternalField(iF, patch().meshPoints());
}


template<class Type>
template<class Type1>
void Foam::pointPatchField<Type>::addToInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF
) const
{
    // Check size
    if (iF.size() != primitiveField().size())
    {
        FatalErrorInFunction
            << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << primitiveField().size()
            << abort(FatalError);
    }

    if (pF.size() != size())
    {
        FatalErrorInFunction
            << "given patch field does not correspond to the mesh. "
            << "Field size: " << pF.size()
            << " mesh size: " << size()
            << abort(FatalError);
    }

    // Get the addressing
    const labelList& mp = patch().meshPoints();

    forAll(mp, pointi)
    {
        iF[mp[pointi]] += pF[pointi];
    }
}


template<class Type>
template<class Type1>
void Foam::pointPatchField<Type>::addToInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF,
    const labelList& points
) const
{
    // Check size
    if (iF.size() != primitiveField().size())
    {
        FatalErrorInFunction
            << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << primitiveField().size()
            << abort(FatalError);
    }

    if (pF.size() != size())
    {
        FatalErrorInFunction
            << "given patch field does not correspond to the mesh. "
            << "Field size: " << pF.size()
            << " mesh size: " << size()
            << abort(FatalError);
    }

    // Get the addressing
    const labelList& mp = patch().meshPoints();

    forAll(points, i)
    {
        label pointi = points[i];
        iF[mp[pointi]] += pF[pointi];
    }
}


template<class Type>
template<class Type1>
void Foam::pointPatchField<Type>::setInInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF,
    const labelList& meshPoints
) const
{
    // Check size
    if (iF.size() != primitiveField().size())
    {
        FatalErrorInFunction
            << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << primitiveField().size()
            << abort(FatalError);
    }

    if (pF.size() != meshPoints.size())
    {
        FatalErrorInFunction
            << "given patch field does not correspond to the meshPoints. "
            << "Field size: " << pF.size()
            << " meshPoints size: " << size()
            << abort(FatalError);
    }

    forAll(meshPoints, pointi)
    {
        iF[meshPoints[pointi]] = pF[pointi];
    }
}


template<class Type>
template<class Type1>
void Foam::pointPatchField<Type>::setInInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF
) const
{
    setInInternalField(iF, pF, patch().meshPoints());
}

template<class Type>
void Foam::pointPatchField<Type>::updateFromInternalField()
{
    //nothing to do; inserted for completeness.
}

template<class Type>
void Foam::pointPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated_)
    {
        updateCoeffs();
    }

    updated_ = false;
}

template<class Type>
void Foam::pointPatchField<Type>::manipulateMatrix
(
    const Field<Type>&,
    labelList&,
    Field<Type>&
)
{
    // do nothing
}

template<class Type>
void Foam::pointPatchField<Type>::manipulateMatrix
(
    labelList&,
    Field<Type>&
)
{
    // do nothing
}

template<class Type>
void Foam::pointPatchField<Type>::setField(const Field<Type>&, Field<Type>&)
{
    // do nothing
}

template<class Type>
void Foam::pointPatchField<Type>::setField(Field<Type>&)
{
    // do nothing
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const pointPatchField<Type>& ptf
)
{
    ptf.write(os);

    os.check(FUNCTION_NAME);

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/pointPatchFields/pointPatchField/pointPatchFieldNew.C"

// ************************************************************************* //

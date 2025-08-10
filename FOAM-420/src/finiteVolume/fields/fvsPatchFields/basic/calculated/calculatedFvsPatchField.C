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

\*---------------------------------------------------------------------------*/

#include "fields/fvsPatchFields/basic/calculated/calculatedFvsPatchField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "fields/GeometricFields/GeometricField/GeometricFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::fvsPatchField<Type>::calculatedType()
{
    return calculatedFvsPatchField<Type>::typeName;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(p, iF)
{}


template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fvsPatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const calculatedFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvsPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const calculatedFvsPatchField<Type>& ptf
)
:
    fvsPatchField<Type>(ptf)
{}


template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const calculatedFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(ptf, iF)
{}


template<class Type>
Foam::tmp<Foam::fvsPatchField<Type>>
Foam::fvsPatchField<Type>::NewCalculatedType
(
    const fvPatch& p
)
{
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTable_().find(p.type());

    if (patchTypeCstrIter != patchConstructorTable_().end())
    {
        return patchTypeCstrIter->second
        (
            p,
            DimensionedField<Type, surfaceMesh>::null()
        );
    }
    else
    {
        return tmp<fvsPatchField<Type>>
        (
            new calculatedFvsPatchField<Type>
            (
                p,
                DimensionedField<Type, surfaceMesh>::null()
            )
        );
    }
}


template<class Type>
template<class Type2>
Foam::tmp<Foam::fvsPatchField<Type>>
Foam::fvsPatchField<Type>::NewCalculatedType
(
    const fvsPatchField<Type2>& pf
)
{
    return NewCalculatedType(pf.patch());
}

// ************************************************************************* //

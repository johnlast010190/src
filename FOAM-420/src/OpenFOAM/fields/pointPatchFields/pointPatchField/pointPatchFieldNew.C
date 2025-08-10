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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
{
    if (debug)
    {
        InfoInFunction << "Constructing pointPatchField<Type>" << endl;
    }

    const auto ctor = ctorTableLookup("patchFieldType type", pointPatchConstructorTable_(), patchFieldType);
    autoPtr<pointPatchField<Type>> pfPtr(ctor(p, iF));

    if
    (
        actualPatchType == word::null
     || actualPatchType != p.type()
    )
    {
        if (pfPtr().constraintType() != p.constraintType())
        {
            // Use default constraint type
            typename pointPatchConstructorTable::iterator patchTypeCstrIter =
                pointPatchConstructorTable_().find(p.type());

            if (patchTypeCstrIter == pointPatchConstructorTable_().end())
            {
                FatalErrorInFunction
                    << "inconsistent patch and patchField types for \n"
                    << "    patch type " << p.type()
                    << " and patchField type " << patchFieldType
                    << exit(FatalError);
            }

            return patchTypeCstrIter->second(p, iF);
        }
    }
    else
    {
        if (pointPatchConstructorTable_().count(p.type()) > 0)
        {
            pfPtr().patchType() = actualPatchType;
        }
    }

    return pfPtr;
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const word& patchFieldType,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing pointPatchField<Type>" << endl;
    }

    word patchFieldType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTable_().find(patchFieldType);

    if (cstrIter == dictionaryConstructorTable_().end())
    {
        if (!disallowGenericPointPatchField)
        {
            cstrIter = dictionaryConstructorTable_().find("generic");
        }

        if (cstrIter == dictionaryConstructorTable_().end())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl << nl
                << exit(FatalIOError);
        }
    }

    // Construct (but not necesarily returned)
    autoPtr<pointPatchField<Type>> pfPtr(cstrIter->second(p, iF, dict));

    if
    (
       !dict.found("patchType")
     || word(dict.lookup("patchType")) != p.type()
    )
    {
        if (pfPtr().constraintType() == p.constraintType())
        {
            // Compatible (constraint-wise) with the patch type
            return pfPtr;
        }
        else
        {
            // Use default constraint type
            typename dictionaryConstructorTable::iterator patchTypeCstrIter
                = dictionaryConstructorTable_().find(p.type());

            if (patchTypeCstrIter == dictionaryConstructorTable_().end())
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "inconsistent patch and patchField types for \n"
                    << "    patch type " << p.type()
                    << " and patchField type " << patchFieldType
                    << exit(FatalIOError);
            }

            return patchTypeCstrIter->second(p, iF, dict);
        }
    }

    return cstrIter->second(p, iF, dict);
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const pointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        InfoInFunction << "Constructing pointPatchField<Type>" << endl;
    }

    const auto ctor = ctorTableLookup("patchField type", patchMapperConstructorTable_(), ptf.type());
    return ctor(ptf, p, iF, pfMapper);
}


// ************************************************************************* //

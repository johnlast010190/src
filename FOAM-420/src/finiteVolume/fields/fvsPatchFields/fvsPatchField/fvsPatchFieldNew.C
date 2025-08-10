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
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvsPatchField<Type>> Foam::fvsPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvsPatchField<Type>" << endl;
    }

    const auto ctor = ctorTableLookup("patchField type", patchConstructorTable_(), patchFieldType);
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTable_().find(p.type());

    if
    (
        actualPatchType == word::null
     || actualPatchType != p.type()
    )
    {
        if (patchTypeCstrIter != patchConstructorTable_().end())
        {
            return patchTypeCstrIter->second(p, iF);
        }
        else
        {
            return ctor(p, iF);
        }
    }
    else
    {
        tmp<fvsPatchField<Type>> tfvp = ctor(p, iF);

        // Check if constraint type override and store patchType if so
        if ((patchTypeCstrIter != patchConstructorTable_().end()))
        {
            tfvp.ref().patchType() = actualPatchType;
        }
        return tfvp;
    }
}


template<class Type>
Foam::tmp<Foam::fvsPatchField<Type>> Foam::fvsPatchField<Type>::New
(
    const word& patchFieldType,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::tmp<Foam::fvsPatchField<Type>> Foam::fvsPatchField<Type>::New
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvsPatchField<Type>" << endl;
    }

    const word patchFieldType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTable_().find(patchFieldType);

    if (cstrIter == dictionaryConstructorTable_().end())
    {
        if (!disallowGenericFvsPatchField)
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

    if
    (
        !dict.found("patchType")
     || word(dict.lookup("patchType")) != p.type()
    )
    {
        typename dictionaryConstructorTable::iterator patchTypeCstrIter
            = dictionaryConstructorTable_().find(p.type());

        if (patchTypeCstrIter != dictionaryConstructorTable_().end() &&
            patchTypeCstrIter->second != cstrIter->second)
        {
            tmp<fvsPatchField<Type>> tpf(cstrIter->second(p, iF, dict));
            typename dictionaryConstructorTable::iterator constraintTypeCstrIter
                = dictionaryConstructorTable_().find(tpf->constraintType());

            if (patchTypeCstrIter == constraintTypeCstrIter)
            {
                return tpf;
            }

            FatalIOErrorInFunction
            (
                dict
            )   << "inconsistent patch and patchField types for \n"
                   "    patch type " << p.type()
                << " and patchField type " << patchFieldType
                << exit(FatalIOError);
        }
    }

    return cstrIter->second(p, iF, dict);
}


template<class Type>
Foam::tmp<Foam::fvsPatchField<Type>> Foam::fvsPatchField<Type>::New
(
    const fvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvsPatchField<Type>" << endl;
    }

    typename patchMapperConstructorTable::iterator cstrIter =
        patchMapperConstructorTable_().find(ptf.type());

    if (cstrIter == patchMapperConstructorTable_().end())
    {
        FatalErrorInFunction
            << "Unknown patchField type " << ptf.type() << nl << nl
            << exit(FatalError);
    }

    typename patchMapperConstructorTable::iterator
        patchTypeCstrIter = patchMapperConstructorTable_().find(p.type());

    if (patchTypeCstrIter != patchMapperConstructorTable_().end())
    {
        return patchTypeCstrIter->second(ptf, p, iF, pfMapper);
    }
    else
    {
        return cstrIter->second(ptf, p, iF, pfMapper);
    }
}


// ************************************************************************* //

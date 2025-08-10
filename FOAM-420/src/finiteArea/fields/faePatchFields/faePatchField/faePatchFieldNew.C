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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faePatchField<Type>> faePatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF
)
{
    if (debug)
    {
        Info<< "faePatchField<Type>::New(const word&, const faPatch&, "
               "const DimensionedField<Type, faEdgeMesh>&) : "
               "constructing faePatchField<Type>"
            << endl;
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
        tmp<faePatchField<Type>> tfap = ctor(p, iF);

        // Check if constraint type override and store patchType if so
        if ((patchTypeCstrIter != patchConstructorTable_().end()))
        {
            tfap.ref().patchType() = actualPatchType;
        }
        return tfap;
    }

}

template<class Type>
tmp<faePatchField<Type>> faePatchField<Type>::New
(
    const word& patchFieldType,
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}

template<class Type>
tmp<faePatchField<Type>> faePatchField<Type>::New
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const dictionary& dict
)
{
    if (debug)
    {
        Info<< "faePatchField<Type>::New(const faPatch&, "
               "const DimensionedField<Type, faEdgeMesh>&, "
               "const dictionary&) : "
               "constructing faePatchField<Type>"
            << endl;
    }

    word patchFieldType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTable_().find(patchFieldType);

    if (cstrIter == dictionaryConstructorTable_().end())
    {
        if (!disallowGenericFaePatchField)
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

    typename dictionaryConstructorTable::iterator patchTypeCstrIter
        = dictionaryConstructorTable_().find(p.type());

    if
    (
        patchTypeCstrIter != dictionaryConstructorTable_().end()
     && *patchTypeCstrIter != *cstrIter
    )
    {
        tmp<faePatchField<Type>> tpf(cstrIter->second(p, iF, dict));
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

    return cstrIter->second(p, iF, dict);
}


// Return a pointer to a new patch created on freestore from
// a given faePatchField<Type> mapped onto a new patch

template<class Type>
tmp<faePatchField<Type>> faePatchField<Type>::New
(
    const faePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const faPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        Info<< "faePatchField<Type>::New(const faePatchField<Type>&,"
               " const faPatch&, const DimensionedField<Type, faEdgeMesh>&, "
               "const faPatchFieldMapper&) : "
               "constructing faePatchField<Type>"
            << endl;
    }

    typename patchMapperConstructorTable::iterator cstrIter =
        patchMapperConstructorTable_().find(ptf.type());

    if (cstrIter == patchMapperConstructorTable_().end())
    {
        FatalErrorInFunction
            << "Unknown patchField type " << ptf.type() << nl << nl
            << "Valid patchField types are :" << endl
            << patchMapperConstructorTable_().sortedToc()
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

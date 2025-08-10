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
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faPatchField<Type>> faPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
{
    if (debug)
    {
        Info<< "faPatchField<Type>::New(const word&, const word&, "
               "const faPatch&, const DimensionedField<Type, areaMesh>&) :"
               " patchFieldType="
            << patchFieldType
            << " : " << p.type()
            << endl;
    }

    typename patchConstructorTable::iterator cstrIter =
        patchConstructorTable_().find(patchFieldType);

    if (cstrIter == patchConstructorTable_().end())
    {
        FatalErrorInFunction
            << "Unknown patchField type "
            << patchFieldType << nl << nl
            << exit(FatalError);
    }

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
            return cstrIter->second(p, iF);
        }
    }
    else
    {
        tmp<faPatchField<Type>> tfap = cstrIter->second(p, iF);

        // Check if constraint type override and store patchType if so
        if ((patchTypeCstrIter != patchConstructorTable_().end()))
        {
            tfap.ref().patchType() = actualPatchType;
        }
        return tfap;
    }

}

template<class Type>
tmp<faPatchField<Type>> faPatchField<Type>::New
(
    const word& patchFieldType,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}



template<class Type>
tmp<faPatchField<Type>> faPatchField<Type>::New
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
{
    if (debug)
    {
        Info<< "faPatchField<Type>::New(const faPatch&, "
               "const DimensionedField<Type, areaMesh>&, const dictionary&) : "
               "constructing faPatchField<Type>"
            << endl;
    }

    word patchFieldType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTable_().find(patchFieldType);

    if (cstrIter == dictionaryConstructorTable_().end())
    {
        if (!disallowGenericFaPatchField)
        {
            cstrIter = dictionaryConstructorTable_().find("default");
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

        if
        (
            patchTypeCstrIter != dictionaryConstructorTable_().end()
         && patchTypeCstrIter->second != cstrIter->second
        )
        {
            tmp<faPatchField<Type>> tpf(cstrIter->second(p, iF, dict));
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
tmp<faPatchField<Type>> faPatchField<Type>::New
(
    const faPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        Info<< "faPatchField<Type>::New(const faPatchField<Type>&,"
               " const faPatch&, const DimensionedField<Type, areaMesh>&, "
               "const faPatchFieldMapper&) : "
               "constructing faPatchField<Type>"
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

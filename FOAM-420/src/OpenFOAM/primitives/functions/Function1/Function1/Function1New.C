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
    (v) 2018-2023 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Constant/Constant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::New
(
    const word& entryName,
    const entry* eptr,
    const dictionary& dict,
    const word& redirectType,
    const bool mandatory
)
{
    word modelType(redirectType);

    const dictionary* coeffs = (eptr ? eptr->dictPtr() : nullptr);

    if (coeffs)
    {
        // Dictionary entry

        DebugInFunction
            << "For " << entryName << " with dictionary entries: "
            << flatOutput(coeffs->toc()) << nl;

        // The "type" entry - mandatory if no redirectType provided
        coeffs->readEntry
        (
            "type",
            modelType,
            keyType::LITERAL,
            (
                modelType.empty() ? true : false
            )
        );

        // Fallthrough
    }
    else if (eptr)
    {
        // Primitive entry
        // - word : the modelType
        // - non-word : value for constant function

        DebugInFunction
            << "For " << entryName << " with primitive entry" << nl;

        ITstream& is = eptr->stream();

        if (is.peek().isWord())
        {
            modelType = is.peek().wordToken();
        }
        else
        {
            // A value - compatibility for reading constant

            const Type constValue = pTraits<Type>(is);

            return autoPtr<Function1<Type>>
            (
                new Function1Types::Constant<Type>
                (
                    entryName,
                    constValue
                )
            );
        }

        // Fallthrough
    }


    if (modelType.empty())
    {
        // Entry missing

        if (mandatory)
        {
            FatalIOErrorInFunction(dict)
                << "Missing or invalid Function1 entry: "
                << entryName << nl
                << exit(FatalIOError);
        }

        return nullptr;
    }
    else if (!coeffs)
    {
        // Primitive entry. Coeffs dictionary is optional.

        const word& kw =
        (
            eptr
          ? eptr->keyword()  // Could be a compatibility lookup
          : entryName
        );

        coeffs = &dict.optionalSubDict(kw + "Coeffs");
    }


    const auto ctor = ctorTableLookup
        (
            "Function1 type", dictionaryConstructorTable_(), modelType
        );

    return ctor(entryName, *coeffs);
}



template<class Type>
Foam::autoPtr<Foam::Function1<Type>> Foam::Function1<Type>::New
(
    const word& entryName,
    const dictionary& dict
)
{
    if (dict.isDict(entryName))
    {
        const dictionary& coeffsDict(dict.subDict(entryName));

        const word Function1Type(coeffsDict.lookup("type"));

        const auto ctor = ctorTableLookup("Function1 type", dictionaryConstructorTable_(), Function1Type);
        return ctor(entryName, coeffsDict);
    }
    else
    {
        Istream& is(dict.lookup(entryName, false));

        token firstToken(is);
        word Function1Type;

        if (!firstToken.isWord())
        {
            is.putBack(firstToken);
            return autoPtr<Function1<Type>>
            (
                new Function1Types::Constant<Type>(entryName, is)
            );
        }
        else
        {
            Function1Type = firstToken.wordToken();
        }

        const auto ctor = ctorTableLookup("Function1 type", dictionaryConstructorTable_(), Function1Type);
        return ctor(
            entryName,
            dict.found(entryName + "Coeffs")
          ? dict.subDict(entryName + "Coeffs")
          : dict
        );
    }
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::NewIfPresent
(
    const word& entryName,
    const dictionary& dict
)
{
    // mandatory = false
    return Function1<Type>::New(entryName, dict);
}


template<class Type>
Foam::refPtr<Foam::Function1<Type>>
Foam::Function1<Type>::New
(
    HashPtrTable<Function1<Type>>& cache,
    const word& entryName,
    const dictionary& dict,
    enum keyType::option matchOpt,
    const bool mandatory
)
{
    // Use the dictionary to find the keyword (allowing wildcards).
    // Alternative would be to have
    // a HashTable where the key type uses a wildcard match

    refPtr<Function1<Type>> fref;

    // Try for direct cache hit
    fref.cref(cache.get(entryName));

    if (fref)
    {
        return fref;
    }

    // Lookup from dictionary
    const entry* eptr = dict.lookupEntryPtr(entryName, true, true);

    if (eptr)
    {
        // Use keyword (potentially a wildcard) instead of entry name
        const auto& kw = eptr->keyword();

        // Try for a cache hit
        fref.cref(cache.get(kw));

        if (!fref)
        {
            // Create new entry
            auto fauto
            (
                Function1<Type>::New
                (
                    kw,
                    eptr,  // Already resolved
                    dict,
                    word::null,
                    mandatory
                )
            );

            if (fauto)
            {
                // Cache the newly created function
                fref.cref(fauto.rawPtr());
                cache.set(kw, fauto);
            }
        }
    }

    if (mandatory && !fref)
    {
        FatalIOErrorInFunction(dict)
            << "No match for " << entryName << nl
            << exit(FatalIOError);
    }

    return fref;
}


// ************************************************************************* //

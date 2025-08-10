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

Description
    Macros to ease declaration of run-time selection tables.

    declareRunTimeSelectionTable is used to create a run-time selection table
    for a base-class which holds constructor pointers on the table.

    declareRunTimeNewSelectionTable is used to create a run-time selection
    table for a derived-class which holds "New" pointers on the table.

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/token/token.H"

#ifndef runTimeSelectionTables_H
#define runTimeSelectionTables_H

#include "memory/autoPtr/autoPtr.H"
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Declare a run-time selection
#define declareRunTimeSelectionTable(autoPtr,baseType,argNames,argList,parList)\
                                                                               \
    /* Construct from argList function pointer type */                         \
    typedef autoPtr<baseType> (*argNames##ConstructorPtr)argList;              \
                                                                               \
    /* Construct from argList function table type */                           \
    typedef std::unordered_map<word, argNames##ConstructorPtr, string::hash>   \
        argNames##ConstructorTable;                                            \
                                                                               \
    /* Construct from argList function pointer table pointer */                \
    static argNames##ConstructorTable& argNames##ConstructorTable_();          \
                                                                               \
    /* Class to add constructor from argList to table */                       \
    template<class baseType##Type>                                             \
    class add##argNames##ConstructorToTable                                    \
    {                                                                          \
    public:                                                                    \
        static autoPtr<baseType> New argList                                   \
        {                                                                      \
            return autoPtr<baseType>(new baseType##Type parList);              \
        }                                                                      \
                                                                               \
        add##argNames##ConstructorToTable                                      \
        (                                                                      \
            const word& lookup = baseType##Type::typeName                      \
        )                                                                      \
        {                                                                      \
            if (!argNames##ConstructorTable_().insert(                 \
                std::make_pair(lookup, New)                                    \
            ).second) {                                                        \
                std::cerr<< "Duplicate entry " << lookup                       \
                    << " in runtime selection table " << #baseType             \
                    << std::endl;                                              \
                error::safePrintStack(std::cerr);                              \
            }                                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    /* Class to add constructor from argList to table */                       \
    /* Remove only the entry (not the table) upon destruction */               \
    template<class baseType##Type>                                             \
    class addRemovable##argNames##ConstructorToTable                           \
    {                                                                          \
        /* retain lookup name for later removal */                             \
        const word& lookup_;                                                   \
                                                                               \
    public:                                                                    \
                                                                               \
        static autoPtr<baseType> New argList                                   \
        {                                                                      \
            return autoPtr<baseType>(new baseType##Type parList);              \
        }                                                                      \
                                                                               \
        addRemovable##argNames##ConstructorToTable                             \
        (                                                                      \
            const word& lookup = baseType##Type::typeName                      \
        )                                                                      \
        :                                                                      \
            lookup_(lookup)                                                    \
        {                                                                      \
            argNames##ConstructorTable_().insert_or_assign(lookup, New);       \
        }                                                                      \
                                                                               \
        ~addRemovable##argNames##ConstructorToTable()                          \
        {                                                                      \
            argNames##ConstructorTable_().erase(lookup_);                      \
        }                                                                      \
    };



//- Declare a run-time selection for derived classes
#define declareRunTimeNewSelectionTable(                                       \
    autoPtr,baseType,argNames,argList,parList)                                 \
                                                                               \
    /* Construct from argList function pointer type */                         \
    typedef autoPtr<baseType> (*argNames##ConstructorPtr)argList;              \
                                                                               \
    /* Construct from argList function table type */                           \
    typedef std::unordered_map<word, argNames##ConstructorPtr, string::hash>   \
        argNames##ConstructorTable;                                            \
                                                                               \
    /* Construct from argList function pointer table pointer */                \
    static argNames##ConstructorTable& argNames##ConstructorTable_();          \
                                                                               \
    /* Class to add constructor from argList to table */                       \
    template<class baseType##Type>                                             \
    class add##argNames##ConstructorToTable                                    \
    {                                                                          \
    public:                                                                    \
                                                                               \
        static autoPtr<baseType> New##baseType argList                         \
        {                                                                      \
            return autoPtr<baseType>(baseType##Type::New parList.ptr());       \
        }                                                                      \
                                                                               \
        add##argNames##ConstructorToTable                                      \
        (                                                                      \
            const word& lookup = baseType##Type::typeName                      \
        )                                                                      \
        {                                                                      \
            auto p = argNames##ConstructorTable_().insert(                     \
                std::make_pair(lookup, New##baseType)                          \
            );                                                                 \
            if (!p.second) {                                                   \
                std::cerr<< "Duplicate entry " << lookup                       \
                    << " in runtime selection table " << #baseType             \
                    << std::endl;                                              \
                error::safePrintStack(std::cerr);                              \
            }                                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    /* Class to add constructor from argList to table */                       \
    template<class baseType##Type>                                             \
    class addRemovable##argNames##ConstructorToTable                           \
    {                                                                          \
        /* retain lookup name for later removal */                             \
        const word& lookup_;                                                   \
                                                                               \
    public:                                                                    \
                                                                               \
        static autoPtr<baseType> New##baseType argList                         \
        {                                                                      \
            return autoPtr<baseType>(baseType##Type::New parList.ptr());       \
        }                                                                      \
                                                                               \
        addRemovable##argNames##ConstructorToTable                             \
        (                                                                      \
            const word& lookup = baseType##Type::typeName                      \
        )                                                                      \
        :                                                                      \
            lookup_(lookup)                                                    \
        {                                                                      \
            argNames##ConstructorTable_().insert_or_assign                     \
            (                                                                  \
                lookup,                                                        \
                New##baseType                                                  \
            );                                                                 \
        }                                                                      \
                                                                               \
        ~addRemovable##argNames##ConstructorToTable()                          \
        {                                                                      \
            argNames##ConstructorTable_().erase(lookup_);                      \
        }                                                                      \
    };

#define defineRunTimeSelectionTablePtr(baseType,argNames)                      \
                                                                               \
    /* Define the constructor function table */                                \
    baseType::argNames##ConstructorTable& baseType::argNames##ConstructorTable_() { \
        static baseType::argNames##ConstructorTable table;                     \
        return table;                                                          \
    }



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define run-time selection table
#define defineRunTimeSelectionTable(baseType,argNames)                         \
                                                                               \
    defineRunTimeSelectionTablePtr(baseType,argNames);


//- Define run-time selection table for template classes
//  use when baseType doesn't need a template argument (eg, is a typedef)
#define defineTemplateRunTimeSelectionTable(baseType,argNames)                 \
                                                                               \
    template<>                                                                 \
    defineRunTimeSelectionTablePtr(baseType,argNames);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- Create pointer to hash-table of functions
//  use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionTablePtr(baseType,argNames,Targ)        \
                                                                               \
    /* Define the constructor function table */                                \
    baseType<Targ>::argNames##ConstructorTable&                                \
        baseType<Targ>::argNames##ConstructorTable_() {                        \
        static baseType<Targ>::argNames##ConstructorTable table;               \
        return table;                                                          \
    };


//- Define run-time selection table for template classes
//  use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionTable(baseType,argNames,Targ)           \
                                                                               \
    template<>                                                                 \
    defineTemplatedRunTimeSelectionTablePtr(baseType,argNames,Targ);


// Pretty-print the allowed keys of a constructor table
template<typename T>
void printCtorTableKeys(Foam::OSstream& s, const T& table)
{
    // Sorted keyset of an unordered map is "fun".
    std::vector<std::string> keys;
    for (const auto& p : table)
    {
        keys.push_back(p.first);
    }
    std::sort(keys.begin(), keys.end());

    for (const auto& k : keys)
    {
        s << "    ";
        s.writeQuoted(k, false);
        s << Foam::nl;
    }
}


//- Look something up in a constructor table, providing a pretty error message
//  if unsuccessful.
template<typename T>
auto ctorTableLookup
(
    const std::string& flavour, const T& table, const std::string& key
)
{
    auto cstrIter = table.find(key);
    if (cstrIter == table.end())
    {
        std::stringstream ss;
        Foam::OSstream& s(FatalErrorInFunction);
        s  << "Unknown " << flavour << ": " << key << '.' << Foam::nl
           << "Valid options are:" << Foam::nl;
        printCtorTableKeys(s, table);
        s << exit(Foam::FatalError);
    }

    return cstrIter->second;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

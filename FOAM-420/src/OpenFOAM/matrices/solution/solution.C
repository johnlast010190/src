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
    (c) 2019-2022 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "matrices/solution/solution.H"
#include "db/Time/Time.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitchWithName(solution, "solution", 0);
}

// List of sub-dictionaries to rewrite
static const Foam::List<Foam::word> subDictNames
(
    Foam::IStringStream("(preconditioner smoother)")()
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solution::update()
{
    bool doOverwrite = localSolution_.lookupOrDefault("overwrite", false);

    if (found("select"))
    {
        subDict(word(lookup("select"))).merge(localSolution_, doOverwrite);
    }
    else
    {
        this->merge(localSolution_, doOverwrite);
    }
}

void Foam::solution::read(const dictionary& dict)
{
    if (localSolution_.size())
    {
        update();
    }
    if (dict.found("cache"))
    {
        cache_ = dict.subDict("cache");
        caching_ = cache_.lookupOrDefault("active", true);
    }

    if (dict.found("relaxationFactors"))
    {
        const dictionary& relaxDict(dict.subDict("relaxationFactors"));
        bool needsCompat = true;

        if (relaxDict.found("fields") || relaxDict.found("equations"))
        {
            if (relaxDict.found("fields"))
            {
                needsCompat = false;
                fieldRelaxDict_ = relaxDict.subDict("fields");
                fieldRelaxCache_.clear();
            }

            if (relaxDict.found("equations"))
            {
                needsCompat = false;
                eqnRelaxDict_ = relaxDict.subDict("equations");
                eqnRelaxCache_.clear();
            }
        }
        if (needsCompat)
        {
            // backwards compatibility
            fieldRelaxDict_.clear();
            fieldRelaxCache_.clear();

            const wordList entryNames(relaxDict.toc());
            forAll(entryNames, i)
            {
                const word& e = entryNames[i];
                scalar value = readScalar(relaxDict.lookup(e));

                if (e(0, 1) == "p")
                {
                    fieldRelaxDict_.add(e, value);
                }
                else if (e.length() >= 3)
                {
                    if (e(0, 3) == "rho")
                    {
                        fieldRelaxDict_.add(e, value);
                    }
                }

            }

            eqnRelaxDict_ = relaxDict;
            eqnRelaxCache_.clear();
        }

        if (!fieldRelaxDefault_)
        {
            fieldRelaxDefault_.reset
            (
                new Function1Types::Constant<scalar>("default", 0)
            );
        }
        if (!eqnRelaxDefault_)
        {
            eqnRelaxDefault_.reset
            (
                new Function1Types::Constant<scalar>("default", 0)
            );
        }

        if (debug)
        {
            Info<< "Relaxation factors:" << nl
                << "fields: " << fieldRelaxDict_ << nl
                << "equations: " << eqnRelaxDict_ << endl;
        }
    }


    if (dict.found("solvers"))
    {
        solvers_ = dict.subDict("solvers");
        upgradeSolverDict(solvers_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solution::solution
(
    const objectRegistry& obr,
    const fileName& dictName
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            (
                obr.readOpt() == IOobject::MUST_READ
             || obr.readOpt() == IOobject::READ_IF_PRESENT
              ? IOobject::MUST_READ_IF_MODIFIED
              : obr.readOpt()
            ),
            IOobject::NO_WRITE
        )
    ),
    cache_(dictionary::null),
    caching_(false),
    fieldRelaxDict_(dictionary::null),
    eqnRelaxDict_(dictionary::null),
    fieldRelaxDefault_(),
    eqnRelaxDefault_(),
    solvers_(dictionary::null),
    localSolution_(dictionary())
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read(dict());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::solution::upgradeSolverDict
(
    dictionary& dict,
    const bool verbose
)
{
    label nChanged = 0;

    // backward compatibility:
    // recast primitive entries into dictionary entries
    forAllIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            Istream& is = iter().stream();
            word name(is);
            dictionary subdict;

            subdict.add("solver", name);
            subdict <<= dictionary(is);

            // preconditioner and smoother entries can be
            // 1) primitiveEntry w/o settings,
            // 2) or a dictionaryEntry.
            // transform primitiveEntry with settings -> dictionaryEntry
            forAll(subDictNames, dictI)
            {
                const word& dictName = subDictNames[dictI];
                entry* ePtr = subdict.lookupEntryPtr(dictName,false,false);

                if (ePtr && !ePtr->isDict())
                {
                    Istream& is = ePtr->stream();
                    is >> name;

                    if (!is.eof())
                    {
                        dictionary newDict;
                        newDict.add(dictName, name);
                        newDict <<= dictionary(is);

                        subdict.set(dictName, newDict);
                    }
                }
            }

            // write out information to help people adjust to the new syntax
            if (verbose && Pstream::master())
            {
                Info<< "// using new solver syntax:\n"
                    << iter().keyword() << subdict << endl;
            }

            // overwrite with dictionary entry
            dict.set(iter().keyword(), subdict);

            nChanged++;
        }
    }

    return nChanged;
}

void Foam::solution::setLocalSolutionDict(const dictionary& input)
{
    localSolution_.merge(input);

    read(dict());
}


void Foam::solution::resetReadOpt(const readOption& option)
{
    readOpt() = option;
}


bool Foam::solution::cache(const word& name) const
{
    if (caching_)
    {
        if (debug)
        {
            Info<< "Cache: find entry for " << name << endl;
        }

        return cache_.found(name);
    }
    else
    {
        return false;
    }
}


bool Foam::solution::relaxField(const word& name) const
{
    if (debug)
    {
        Info<< "Field relaxation factor for " << name
            << " is " << (fieldRelaxDict_.found(name) ? "set" : "unset")
            << endl;
    }

    return fieldRelaxDict_.found(name) || fieldRelaxDict_.found("default");
}

bool Foam::solution::storeField(const word& name) const
{
    if (debug)
    {
        Info<< "Field relaxation factor for " << name
            << " is " << fieldRelaxationFactor(name)
            << ". Return true for storing it."
            << endl;
    }
    return fieldRelaxationFactor(name) != 1.0;
}


bool Foam::solution::relaxEquation(const word& name) const
{
    if (debug)
    {
        Info<< "Find equation relaxation factor for " << name << endl;
    }

    return eqnRelaxDict_.found(name) || eqnRelaxDict_.found("default");
}


Foam::scalar Foam::solution::fieldRelaxationFactor(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup variable relaxation factor for " << name << endl;
    }

    if (fieldRelaxDict_.found(name))
    {
        return Function1<scalar>::New
        (
            fieldRelaxCache_,  // cache
            name,
            fieldRelaxDict_,
            keyType::REGEX
        )().value(time().timeOutputValue());
    }
    else if (fieldRelaxDefault_)
    {
        scalar value = fieldRelaxDefault_().value(time().timeOutputValue());
        if (value>SMALL)
        {
            return value;
        }
        else
        {
            FatalIOErrorInFunction(fieldRelaxDict_)
                << "Cannot find variable relaxation factor for '" << name
                << "' and default is not set." << nl
                << exit(FatalIOError);
        }
    }

    FatalIOErrorInFunction(fieldRelaxDict_)
        << "Cannot find variable relaxation factor for '" << name
        << "' or a suitable default value." << nl
        << exit(FatalIOError);

    return 0;
}


Foam::scalar Foam::solution::equationRelaxationFactor(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup equation relaxation factor for " << name << endl;
    }

    if (eqnRelaxDict_.found(name))
    {
        return Function1<scalar>::New
        (
            eqnRelaxCache_,  // cache
            name,
            eqnRelaxDict_,
            keyType::REGEX
        )().value(time().timeOutputValue());
    }
    else if (eqnRelaxDefault_)
    {
        scalar value = eqnRelaxDefault_().value(time().timeOutputValue());
        if (value>SMALL)
        {
            return value;
        }
        else
        {
            FatalIOErrorInFunction(eqnRelaxDict_)
                << "Cannot find equation relaxation factor for '" << name
                << "' and default is not set." << nl
                << exit(FatalIOError);
        }
    }

     FatalIOErrorInFunction
     (
         eqnRelaxDict_
     )   << "Cannot find equation relaxation factor for '" << name
         << "' or a suitable default value."
         << exit(FatalIOError);

    return 0;
}


const Foam::dictionary& Foam::solution::dict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }
    else
    {
        return *this;
    }
}


const Foam::dictionary& Foam::solution::solverDict(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup solver for " << name << endl;
    }

    return solvers_.subDict(name);
}


const Foam::dictionary& Foam::solution::solver(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup solver for " << name << endl;
    }

    return solvers_.subDict(name);
}


bool Foam::solution::solverDictExists(const word& name) const
{
    return solvers_.found(name);
}


bool Foam::solution::read()
{
    if (regIOobject::read())
    {
        read(dict());

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

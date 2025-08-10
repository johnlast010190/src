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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Pstreams/Pstream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Should issue warning if there is +ve versioning (+ve version number)
// and the this version number is older than the current OpenFOAM version
// as conveyed by the OPENFOAM compiler define.

static inline constexpr bool shouldWarnVersion(const int version)
{
    return (version > 0 && version < FOAM_API);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary::const_searcher Foam::dictionary::csearchCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    const_searcher finder(csearch(keyword, matchOpt));

    if (finder.found())
    {
        return finder;
    }

    for (const std::pair<const char*,int>& iter : compat)
    {
        finder = csearch(word::validated(iter.first), matchOpt);

        if (finder.found())
        {
            // Only want a single warning (on master), but guard with a
            // parRun check to avoid Pstream::master() when Pstream has not
            // yet been initialized
            if
            (
                shouldWarnVersion(iter.second)
             && (Pstream::parRun() ? Pstream::master() : true)
            )
            {
                DeprecationIOWarningInFunction
                (
                    *this, iter.first, "keyword", iter.second
                )
                    << "Should be '"
                    << keyword.c_str() << "' in dictionary \""
                    << name().c_str() << "\" "
                    << nl
                    << endl;
            }

            break;
        }
    }

    return finder;
}


bool Foam::dictionary::foundCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    return csearchCompat(keyword, compat, matchOpt).found();
}


const Foam::entry* Foam::dictionary::findCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    return csearchCompat(keyword, compat, matchOpt).ptr();
}


const Foam::entry& Foam::dictionary::lookupEntryCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    const const_searcher finder(csearchCompat(keyword, compat, matchOpt));

    if (!finder.found())
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return finder.ref();
}


Foam::ITstream& Foam::dictionary::lookupCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    return lookupEntryCompat(keyword, compat, matchOpt).stream();
}


// ************************************************************************* //

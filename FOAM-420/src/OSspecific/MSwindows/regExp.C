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
    (c) 2016 blueCAPE Ltda.
    (c) 2011-2016 OpenFOAM Foundation

Modifications
    This file is based on the original version for POSIX:
        OpenFOAM/src/OSspecific/POSIX/

    This file has been created by blueCAPE's unofficial mingw patches for
    OpenFOAM.
    For more information about these patches, visit:
        http://bluecfd.com/Core

    Modifications made:
      - Only file renaming.

\*---------------------------------------------------------------------------*/

#include "regExp.H"
#include "primitives/strings/string/string.H"
#include "containers/Lists/List/List.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class StringType>
bool Foam::regExp::matchGrouping
(
    const std::string& str,
    List<StringType>& groups
) const
{
    if (preg_ && str.size())
    {
        size_t nmatch = ngroups() + 1;
        regmatch_t pmatch[nmatch];

        // Also verify that the entire string was matched.
        // pmatch[0] is the entire match
        // pmatch[1..] are the (...) sub-groups
        if
        (
            regexec(preg_, str.c_str(), nmatch, pmatch, 0) == 0
         && (pmatch[0].rm_so == 0 && pmatch[0].rm_eo == label(str.size()))
        )
        {
            groups.setSize(ngroups());
            label groupI = 0;

            for (size_t matchI = 1; matchI < nmatch; matchI++)
            {
                if (pmatch[matchI].rm_so != -1 && pmatch[matchI].rm_eo != -1)
                {
                    groups[groupI] = str.substr
                    (
                        pmatch[matchI].rm_so,
                        pmatch[matchI].rm_eo - pmatch[matchI].rm_so
                    );
                }
                else
                {
                    groups[groupI].clear();
                }
                groupI++;
            }

            return true;
        }
    }

    groups.clear();
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regExp::regExp()
:
    preg_(0)
{}


Foam::regExp::regExp(const char* pattern)
:
    preg_(0)
{
    set(pattern, false);
}


Foam::regExp::regExp(const std::string& pattern)
:
    preg_(0)
{
    set(pattern.c_str(), false);
}

Foam::regExp::regExp(const char* pattern, bool ignoreCase)
:
    preg_(0)
{
    set(pattern, ignoreCase);
}

Foam::regExp::regExp(const std::string& pattern, bool ignoreCase)
:
    preg_(0)
{
    set(pattern, ignoreCase);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regExp::~regExp()
{
    clear();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regExp::set(const char* pattern, const bool ignoreCase) const
{
    clear();

    // Avoid nullptr and zero-length patterns
    if (pattern && *pattern)
    {
        int cflags = REG_EXTENDED;
        if (ignoreCase)
        {
            cflags |= REG_ICASE;
        }

        const char* pat = pattern;

        // Check for embedded prefix for ignore-case
        // this is the only embedded prefix we support
        // - a simple check is sufficient
        if (!strncmp(pattern, "(?i)", 4))
        {
            cflags |= REG_ICASE;
            pat += 4;

            // avoid zero-length patterns
            if (!*pat)
            {
                return;
            }
        }

        preg_ = new regex_t;
        int err = regcomp(preg_, pat, cflags);

        if (err != 0)
        {
            char errbuf[200];
            regerror(err, preg_, errbuf, sizeof(errbuf));

            FatalErrorInFunction
                << "Failed to compile regular expression '" << pattern << "'"
                << nl << errbuf
                << exit(FatalError);
        }
    }
}


void Foam::regExp::set(const std::string& pattern, const bool ignoreCase) const
{
    return set(pattern.c_str(), ignoreCase);
}


bool Foam::regExp::clear() const
{
    if (preg_)
    {
        regfree(preg_);
        delete preg_;
        preg_ = 0;

        return true;
    }

    return false;
}


std::string::size_type Foam::regExp::find(const std::string& str) const
{
    if (preg_ && str.size())
    {
        size_t nmatch = 1;
        regmatch_t pmatch[1];

        if (regexec(preg_, str.c_str(), nmatch, pmatch, 0) == 0)
        {
            return pmatch[0].rm_so;
        }
    }

    return string::npos;
}


bool Foam::regExp::match(const std::string& str) const
{
    if (preg_ && str.size())
    {
        size_t nmatch = 1;
        regmatch_t pmatch[1];

        // Also verify that the entire string was matched
        // pmatch[0] is the entire match
        if
        (
            regexec(preg_, str.c_str(), nmatch, pmatch, 0) == 0
         && (pmatch[0].rm_so == 0 && pmatch[0].rm_eo == label(str.size()))
        )
        {
            return true;
        }
    }

    return false;
}


bool Foam::regExp::match
(
    const std::string& str,
    List<std::string>& groups
) const
{
    return matchGrouping(str, groups);
}


bool Foam::regExp::match
(
    const std::string& str,
    List<Foam::string>& groups
) const
{
    return matchGrouping(str, groups);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::regExp::operator=(const char* pat)
{
    set(pat);
}


void Foam::regExp::operator=(const std::string& pat)
{
    set(pat);
}


// ************************************************************************* //

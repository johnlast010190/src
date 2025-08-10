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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/strings/word/word.H"
#include "global/debug/debug.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::word::typeName = "word";
int Foam::word::debug(Foam::debug::debugSwitch(word::typeName, 0));
const Foam::word Foam::word::null;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::word::validated(const std::string& s)
{
    std::string::size_type count = 0;
    bool prefix = false;

    // Count number of valid characters and detect if the first character
    // happens to be a digit, which we'd like to avoid having since this
    // will cause parse issues when read back later.
    for (std::string::const_iterator it = s.cbegin(); it != s.cend(); ++it)
    {
        const char c = *it;

        if (word::valid(c))
        {
            if (!count && isdigit(c))
            {
                // First valid character was a digit - prefix with '_'
                prefix = true;
                ++count;
            }

            ++count;
        }
    }

    if (count == s.size() && !prefix)
    {
        return word(s, false);  // Already checked, can just return as word
    }

    word out;
    out.resize(count);
    count = 0;

    // Copy valid content.
    if (prefix)
    {
        out[count++] = '_';
    }

    for (std::string::const_iterator it = s.cbegin(); it != s.cend(); ++it)
    {
        const char c = *it;

        if (word::valid(c))
        {
            out[count++] = c;
        }
    }

    out.resize(count);

    return out;
}

Foam::word Foam::word::validate(const std::string& s, const bool prefix)
{
    word out;
    out.resize(s.size() + (prefix ? 1 : 0));

    std::string::size_type len = 0;

    // As per validate, but optionally detect if the first character
    // is a digit, which we'd like to avoid having since this will
    // cause parse issues when read back later.
    for (auto iter = s.cbegin(); iter != s.cend(); ++iter)
    {
        const char c = *iter;

        if (word::valid(c))
        {
            if (!len && prefix && isdigit(c))
            {
                // First valid character was a digit - prefix with '_'
                out[len++] = '_';
            }

            out[len++] = c;
        }
    }

    out.resize(len);

    return out;
}

Foam::word Foam::word::validate
(
    const char* first,
    const char* last,
    const bool prefix
)
{
    std::string::size_type len = (last - first) + (prefix ? 1 : 0);

    word out;
    out.resize(len);

    for (len=0; first != last; ++first)
    {
        const char c = *first;

        if (word::valid(c))
        {
            if (!len && prefix && isdigit(c))
            {
                // First valid character was a digit - prefix with '_'
                out[len++] = '_';
            }

            out[len++] = c;
        }
    }

    out.resize(len);

    return out;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::word::lessExt() const
{
    const size_type i = find_ext();

    if (i == npos)
    {
        return *this;
    }
    else
    {
        return substr(0, i);
    }
}


Foam::word Foam::word::ext() const
{
    return string::ext();
}


Foam::word& Foam::word::ext(const word& ending)
{
    string::ext(ending);
    return *this;
}


bool Foam::word::hasExt(const word& ending) const
{
    return string::hasExt(ending);
}


bool Foam::word::hasExt(const wordRe& ending) const
{
    return string::hasExt(ending);
}


// ************************************************************************* //

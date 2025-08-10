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
    (c) 2017-2018 OpenCFD Ltd.
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline bool Foam::wordRes::found_match
(
    const UList<wordRe>& patterns,
    const std::string& text,
    bool literal
)
{
    for (const wordRe& select : patterns)
    {
        if (select.match(text, literal))
        {
            return true;
        }
    }

    return false;
}


inline Foam::wordRe::compOption Foam::wordRes::found_matched
(
    const UList<wordRe>& patterns,
    const std::string& text
)
{
    auto retval(wordRe::compOption::DETECT);

    for (const wordRe& select : patterns)
    {
        if (!select.isPattern())
        {
            if (select.match(text, true))
            {
                return wordRe::compOption::LITERAL;
            }
        }
        else if (wordRe::compOption::DETECT == retval)
        {
            if (select.match(text, false))
            {
                retval = wordRe::compOption::REGEX;
            }
        }
    }

    return retval;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::wordRes::match(const std::string& text, bool literal) const
{
    return found_match(*this, text, literal);
}


inline Foam::wordRe::compOption
Foam::wordRes::matched(const std::string& text) const
{
    return found_matched(*this, text);
}


template<class StringType>
inline Foam::labelList Foam::wordRes::matching
(
    const UList<StringType>& input,
    const bool invert
) const
{
    const label len = input.size();

    labelList indices(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        if (match(input[i]) ? !invert : invert)
        {
            indices[count] = i;
            ++count;
        }
    }
    indices.resize(count);

    return indices;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wordRes::operator()(const std::string& text) const
{
    return found_match(*this, text);
}


// ************************************************************************* //

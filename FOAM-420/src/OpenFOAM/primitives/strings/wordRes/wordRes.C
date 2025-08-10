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
    (c) 2016-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/strings/wordRes/wordRes.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


Foam::wordRes Foam::wordRes::uniq(const UList<wordRe>& input)
{
    wordRes output(input.size());

    // Use linear List search instead of HashSet, since the lists are
    // normally fairly small and mostly just have unique entries
    // anyhow. This reduces the overall overhead.

    List<bool> duplicate(input.size(), false);  // Track duplicates

    label count = 0;

    forAll(input, i)
    {
        const wordRe& val = input[i];

        const label next = input.find(val, i+1);

        if (next > i)
        {
            duplicate[next] = true;  // Duplicate
        }

        if (!duplicate[i])
        {
            output[count] = val;
            ++count;
        }
    }

    output.resize(count);

    return output;
}


void Foam::wordRes::inplaceUniq(wordReList& input)
{
    wordHashSet uniqWord;

    label count = 0;
    forAll(input, i)
    {
        const auto& select = input[i];

        if (select.isPattern() || uniqWord.insert(select))
        {
            if (count != i)
            {
                input[count] = input[i];
            }
            ++count;
        }
    }

    input.setSize(count);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wordRes::uniq()
{
    List<wordRe> input = *this;

    wordRes& output = *this;

    // Use linear List search instead of HashSet, since the lists are
    // normally fairly small and mostly just have unique entries
    // anyhow. This reduces the overall overhead.

    List<bool> duplicate(input.size(), false);  // Track duplicates

    label count = 0;

    forAll(input, i)
    {
        wordRe& val = input[i];

        const label next = input.find(val, i+1);

        if (next > i)
        {
            duplicate[next] = true;  // Duplicate
        }

        if (!duplicate[i])
        {
            output[count] = std::move(val);
            ++count;
        }
    }

    output.resize(count);
}


// ************************************************************************* //

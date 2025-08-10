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
    (c) 2011 OpenFOAM Foundation
    (c) 2016-2019 OpenCFD Ltd.
    (c) 2020 Esi Ltd.

Class
    Foam::wordRes

Description
    A List of wordRe with additional matching capabilities.

SourceFiles
    wordResI.H
    wordRes.C

\*---------------------------------------------------------------------------*/

#ifndef wordRes_H
#define wordRes_H

#include "primitives/strings/lists/wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wordRes Declaration
\*---------------------------------------------------------------------------*/

class wordRes
:
    public List<wordRe>
{
    // Private Methods

        //- Smart match as literal or regex, stopping on the first match.
        inline static bool found_match
        (
            const UList<wordRe>& patterns,
            const std::string& text,
            bool literal=false
        );

        //- Smart match across entire list, returning the match type.
        //  Stops on the first literal match, or continues to examine
        //  if a regex match occurs.
        //  \return wordRe::LITERAL, wordRe::REGEX on match and
        //      wordRe::UNKNOWN otherwise.
        inline static wordRe::compOption found_matched
        (
            const UList<wordRe>& patterns,
            const std::string& text
        );


public:

    // Public Classes

        //- Functor wrapper for matching against a List of wordRe
        struct matcher
        {
            const UList<wordRe>& values;

            matcher(const UList<wordRe>& list)
            :
                values(list)
            {}

            //- Return true if string matches ANY of the regular expressions
            //  Allows use as a predicate.
            bool operator()(const std::string& text) const
            {
                return found_match(values, text);
            }
        };


    // Factory Methods

        //- Return a wordRes with duplicate entries filtered out.
        //  No distinction made between literals or regular expressions.
        static wordRes uniq(const UList<wordRe>& input);

        // Define inplaceUniq for backwards compatability with FOAM
        static void inplaceUniq(wordReList& input);


    // Constructors

        //- Inherit constructors from List of wordRe
        using List<wordRe>::List;

        // For backwards-compatability with FOAM, construct from UList<wordRe>
        //- Construct from a list of wordRe
        inline wordRes(const UList<wordRe>& list) : List(list) {};

        // Older versions of gcc appear to explicitly require this delegation
        // Specifically g++ version 5.4.0 on Ubuntu 16.04
        inline wordRes() : List() {};


    //- Destructor
    ~wordRes() = default;


    // Member Functions

        //- Filter out duplicate entries (inplace).
        //  No distinction made between literals or regular expressions.
        void uniq();

        //- Smart match as literal or regex, stopping on the first match.
        //
        //  \param literal Force literal match only.
        //  \return True if text matches ANY of the entries.
        inline bool match(const std::string& text, bool literal=false) const;

        //- Smart match in the list of matchers, returning the match type.
        //  It stops if there is a literal match, or continues to examine
        //  other regexs.
        //  \return LITERAL if a lteral match was found, REGEX if a regex
        //      match was found and UNKNOWN otherwise.
        inline wordRe::compOption matched(const std::string& text) const;

        //- Extract list indices for all matches.
        //
        //  \param input  A list of string inputs to match against
        //  \param invert invert the matching logic
        //  \return The locations (indices) in the input list where match()
        //      is true
        template<class StringType>
        inline labelList matching
        (
            const UList<StringType>& input,
            const bool invert=false
        ) const;


    // Member Operators

        //- Identical to match(), for use as a predicate.
        inline bool operator()(const std::string& text) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitives/strings/wordRes/wordResI.H"

#endif

// ************************************************************************* //

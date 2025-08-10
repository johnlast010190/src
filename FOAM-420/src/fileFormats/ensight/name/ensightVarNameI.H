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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include <cctype>


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ensight::VarName::VarName(const VarName& vn)
:
    word(vn, false)
{}


inline Foam::ensight::VarName::VarName(const char* s)
:
    word(s, false)
{
    stripInvalid();
}


inline Foam::ensight::VarName::VarName(const std::string& s)
:
    word(s, false)
{
    stripInvalid();
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::ensight::VarName::stripInvalid()
{
    string::stripInvalid<VarName>(*this);

    if (empty())
    {
        FatalErrorInFunction
            << "ensight::VarName empty after stripping" << nl
            << exit(FatalError);
    }

    // prefix with '_' to avoid starting with leading digits
    std::string::iterator iter = begin();
    if (isdigit(*iter))
    {
        insert(iter, '_');
    }
}


inline bool Foam::ensight::VarName::valid(char c)
{
    return
    (
        word::valid(c) // includes space, quotes, /\;{}
     && c != '!'
     && c != '#'
     && c != '%'
     && c != '('
     && c != ')'
     && c != '*'
     && c != '+'
     && c != ','
     && c != '-'
     && c != '.'
     && c != '@'
     && c != '['
     && c != ']'
     && c != '^'
    );
}


// ************************************************************************* //

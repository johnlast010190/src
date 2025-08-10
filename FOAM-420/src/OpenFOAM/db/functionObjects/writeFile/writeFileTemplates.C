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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/StringStreams/StringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::writeFile::writeHeaderValue
(
    Ostream& os,
    const string& property,
    const Type& value
) const
{
    os  << setw(1) << '#' << setw(1) << ' '
        << setf(ios_base::left) << setw(charWidth() - 2) << property.c_str()
        << setw(1) << ':' << setw(1) << ' ' << value << nl;
}

template<typename T1>
void Foam::functionObjects::writeFile::writeDelimited
(
    Ostream& os,
    const T1& value
) const
{
    unsigned short width = 0;
    if (setWidth_)
    {
        width = charWidth();
    }

    // setw() messes with anything that streams into os directly.  For example,
    // anything that derives from VectorSpace will stream in a single
    // token::BEGIN_LIST, then a single space, then the first value, etc...
    // setw() will operate on the token::BEGIN_LIST, leaving a huge space
    // between the open bracket and the rest of the object you want to print.
    //
    // For example, take a Point.  The space gets put after the first bracket,
    // not at the end, so ( x y z ) becomes (          x y z).
    //
    // This probably isn't the most elegant fix (probably involves making
    // copies, etc...), but should be pretty robust.  This won't ever be a
    // bottleneck, so in this case robustness > elegance.
    OStringStream buf;
    buf << value;

    os  << delimiter_ << std::string(spacesAfterDelimiter_, ' ').c_str()
        << setw(width) << buf.str().c_str();
}

// For backwards-compatibility, don't print "" around strings
template<>
void Foam::functionObjects::writeFile::writeDelimited<Foam::string>
(
    Ostream& os,
    const string& str
) const
{
    writeDelimited(os, str.c_str());
}

// For backwards-compatibility, don't print "" around strings
template<>
void Foam::functionObjects::writeFile::writeDelimited<std::string>
(
    Ostream& os,
    const std::string& str
) const
{
    writeDelimited(os, str.c_str());
}

template<typename T1, typename... Args>
void Foam::functionObjects::writeFile::writeDelimited
(
    Ostream& os,
    const T1& value,
    const Args&... args
) const
{
    writeDelimited(os, value);
    writeDelimited(os, args...);
}

template<typename... Args>
void Foam::functionObjects::writeFile::writeDelimitedComment(
    Ostream& os,
    const string& str,
    const Args&... args
) const
{
    writeDelimitedComment(os, str.c_str(), args...);
}


template<typename T1, typename... Args>
void Foam::functionObjects::writeFile::writeDelimitedComment
(
    Ostream& os,
    const T1& value,
    const Args&... args
) const
{
    // Write a comment and the first argument (in order to get the fixed
    // width size right, if applicable)
    unsigned short width = 0;
    if (setWidth_)
    {
        width = charWidth();
    }
    os  << setw(1) << "#" << setw(1) << ' '
        << setf(ios_base::left) << setw(width - 2) << value;

    // Now write everything else as normal
    writeDelimited(os, args...);
}

// ************************************************************************* //

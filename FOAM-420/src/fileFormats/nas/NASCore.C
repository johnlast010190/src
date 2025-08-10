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

\*---------------------------------------------------------------------------*/

#include "nas/NASCore.H"
#include "db/IOstreams/StringStreams/IStringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::NASCore::NASCore()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::fileFormats::NASCore::parseNASCoord
(
#if defined(WIN64) || defined(WIN32)
    const Foam::string& s
    #else
    const string& s
    #endif
)
{
    size_t expSign = s.find_last_of("+-");

    if (expSign != string::npos && expSign > 0 && !isspace(s[expSign-1]))
    {
        scalar mantissa = readScalar(IStringStream(s.substr(0, expSign))());
        scalar exponent = readScalar(IStringStream(s.substr(expSign+1))());

        if (s[expSign] == '-')
        {
            exponent = -exponent;
        }
        return mantissa * pow(10, exponent);
    }
    else
    {
        return readScalar(IStringStream(s)());
    }
}


// ************************************************************************* //

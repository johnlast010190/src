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
    (c) 2014-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/ints/int32/int32.H"
#include "primitives/strings/stringOps/stringOps.H"
#include "db/IOstreams/IOstreams.H"

#include <inttypes.h>
#include <sstream>
#include <cerrno>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const char* fmt, const int32_t val)
{
    return stringOps::name(fmt, val);
}


Foam::word Foam::name(const std::string& fmt, const int32_t val)
{
    return stringOps::name(fmt, val);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, int32_t& i)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        i = int32_t(t.labelToken());
    }
    else
    {
        is.setBad();
        FatalIOErrorInFunction(is)
            << "wrong token type - expected int32_t, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


int32_t Foam::readInt32(Istream& is)
{
    int32_t val;
    is >> val;

    return val;
}


bool Foam::read(const char* buf, int32_t& s)
{
    char *endptr = nullptr;
    errno = 0;
    intmax_t l = strtoimax(buf, &endptr, 10);
    s = int32_t(l);
    return
        (*endptr == 0) && (errno == 0)
     && (l >= INT32_MIN) && (l <= INT32_MAX);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const int32_t i)
{
    os.write(label(i));
    os.check(FUNCTION_NAME);
    return os;
}


//  32-bit operating systems are no longer supported by FOAM
/*
#if FOAM_ARCH_OPTION == 32
Foam::Istream& Foam::operator>>(Istream& is, long& i)
{
    return operator>>(is, reinterpret_cast<int32_t&>(i));
}

Foam::Ostream& Foam::operator<<(Ostream& os, const long i)
{
    os << int32_t(i);
    return os;
}
#endif*/


// ************************************************************************* //

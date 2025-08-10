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
    (c) 2019-2021 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "simpleObjectRegistry.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Tstreams/ITstream.H"
#include "db/IOstreams/StringStreams/StringStream.H"
#include "primitives/ints/int/int.H"
#include "primitives/Scalar/floatScalar/floatScalar.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleObjectRegistry::setNamedValue
(
    std::string name,
    int val,
    bool report
)
{
    // Report enables output, but respect DetailInfo state as well.
    // The local log variable captures this logic.

    const bool log = (report /*&& Foam::infoDetailLevel > 0*/);

    token tok(static_cast<label>(val));

    // Handle name=value
    const auto eq = name.find('=');

    if (eq != std::string::npos)
    {
        std::string strval(name.substr(eq+1));
        name.erase(eq);  // Truncate the name

        scalar fvalue = Foam::readScalar(IStringStream(strval)());
        int ival = Foam::readInt(IStringStream(strval)());

        if (ival)
        {
            tok = static_cast<label>(ival);
        }
        else if (fvalue > 0.0)
        {
            tok = fvalue;
        }
        // Treat 'name=' like 'name' (ie, no value parameter)
        // silently ignore 'name=junk', but could warn
    }


    simpleObjectRegistryEntry* objPtr = this->lookupPtr(name.c_str());

    if (objPtr)
    {
        // The generic interface requires an Istream.
        ITstream is("", tokenList(1, tok));

        Log << name.c_str() << '=' << tok << nl;

        const List<simpleRegIOobject*>& objects = *objPtr;

        for (simpleRegIOobject* obj : objects)
        {
            is.rewind();
            obj->readData(is);
        }
    }
    else
    {
        Log << name.c_str() << " (unregistered)" << nl;
    }
}


// ************************************************************************* //

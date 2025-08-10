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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2023 Esi Ltd.

Description
    Writes the header description of the File to the stream
    associated with the File.

\*---------------------------------------------------------------------------*/

#include "db/IOobject/IOobject.H"
#include "db/objectRegistry/objectRegistry.H"
#include "global/foamVersion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::IOobject::writeBanner(Ostream& os, bool noHint)
{
    static bool spacesSet(false);
    static char spaces[40];

/*
    const char* cmd =
        "/sbin/ifconfig eth0 | grep -o -E '([[:xdigit:]]{1,2}:){5}[[:xdigit:]]{1,2}'";
    FILE* pipe = popen(cmd, "r");
    if (!pipe)
    {
        std::exit(0);
    }

    char buffer[128];
    Foam::string result("");

    while (!feof(pipe))
    {
        if (fgets(buffer, 128, pipe) != nullptr)
        {
            result += buffer;
        }
    }
    pclose(pipe);
    string mac("00:25:90:4C:C9:7E");

    if (string(result).substr(0, 17) != mac || Foam::hostName() != "gpu02")
    {
        std::exit(0);
    }
*/

    // Capitalise "dev" version
    const char* FOAMversion =
        (strcmp(FOAMversion, "dev") == 0)
      ? "Dev"
      : Foam::FOAMversion;

    if (!spacesSet)
    {
        memset(spaces, ' ', 40);

        size_t len = strlen(FOAMversion);
        if (len < 38)
        {
            spaces[38 - len] = '\0';
        }
        else
        {
            spaces[0] = '\0';
        }
        spacesSet = true;
    }

    if (noHint)
    {
        os  <<
            "/*--------------------------------------"
            "-------------------------------------*\\\n";
    }
    else
    {
        os  <<
            "/*--------------------------------*- C++ "
            "-*----------------------------------*\\\n";
    }
    os  <<
        "|       o        |                                                            |\n"
        "|    o     o     |  FOAM (R) : Open-source CFD for Enterprise                |\n"
        "|   o   O   o    |  Version : " << FOAMversion << spaces << "          |\n"
        "|    o     o     |  ESI Ltd. <http://esi.com/>                            |\n"
        "|       o        |                                                            |\n"
        "\\*---------------------------------------------------------------------------*/\n";

    return os;
}


Foam::Ostream& Foam::IOobject::writeDivider(Ostream& os)
{
    os  <<
        "// * * * * * * * * * * * * * * * * * "
        "* * * * * * * * * * * * * * * * * * * * //\n";

    return os;
}

Foam::Ostream& Foam::IOobject::writeEndDivider(Ostream& os)
{
    os  << "\n\n"
        "// *****************************************"
        "******************************** //\n";

    return os;
}


bool Foam::IOobject::writeHeader(Ostream& os, const word& type) const
{
    if (!os.good())
    {
        InfoInFunction
            << "No stream open for write" << nl
            << os.info() << endl;

        return false;
    }

    writeBanner(os)
        << "FoamFile\n{\n"
        << "    version     " << os.version() << ";\n"
        << "    format      " << os.format() << ";\n"
        << "    class       " << type << ";\n";

    if (os.format() == IOstream::BINARY)
    {
        os  << "    arch        " << Foam::FOAMbuildArch << ";\n";
    }

    if (!note().empty())
    {
        os  << "    note        " << note() << ";\n";
    }

    os  << "    object      " << name() << ";\n"
        << "}" << nl;

    writeDivider(os) << nl;

    return true;
}


bool Foam::IOobject::writeHeader(Ostream& os) const
{
    return writeHeader(os, type());
}


// ************************************************************************* //

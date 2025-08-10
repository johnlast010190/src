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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "printDictionary/printDictionary.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(printDictionary, 0);
    addToRunTimeSelectionTable(functionObject, printDictionary, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::printDictionary::printDictionary
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::printDictionary::~printDictionary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::printDictionary::execute()
{
    // Do nothing - only valid on write
    return true;
}


bool Foam::functionObjects::printDictionary::write()
{
    // Do nothing - only valid on write
    return true;
}


bool Foam::functionObjects::printDictionary::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    //read dictionaries and print them to screen
    Log << type() << " " << name() <<  " read:" << nl;

    Log << nl << "    OpenFOAM dictionary settings:" << nl
         << "    #############################" << nl << endl;

    //will only check system and constant for dictionaries

    List<word> dicts(dict.lookup("dictionaries"));

    Log << "    System:" << nl << endl;
    forAll(dicts, dI)
    {
        IOobject header
        (
            dicts[dI],
            obr_.time().system(),
            obr_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (header.typeHeaderOk<IOdictionary>(true))
        {
            //found dict in system, read and print to screen
            IOdictionary cDict(header);

            Log << "    " << dicts[dI];
            Log << cDict << endl;
            Log << endl;
        }
    }

    Log << "    Constant:" << nl << endl;
    forAll(dicts, dI)
    {
        IOobject header
        (
            dicts[dI],
            obr_.time().constant(),
            obr_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (header.typeHeaderOk<IOdictionary>(true))
        {
            //found dict in system, read and print to screen
            IOdictionary cDict(header);

            Log << "    " << cDict << endl;

            Log << endl;
        }
    }

    Log << endl;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

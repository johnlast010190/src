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
    (c) 2015 OpenCFD Ltd.
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cloudInfo/cloudInfo.H"
#include "clouds/baseClasses/kinematicCloud/kinematicCloud.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Pstreams/PstreamReduceOps.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudInfo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        cloudInfo,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::cloudInfo::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Cloud information");
    writeCommented(os, "Time");
    writeDelimited(os, "nParcels");
    writeDelimited(os, "mass");
    writeDelimited(os, "Dmax");
    writeDelimited(os, "D10");
    writeDelimited(os, "D32");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudInfo::cloudInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    logFiles(obr_, name, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudInfo::~cloudInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::cloudInfo::read(const dictionary& dict)
{
    if (regionFunctionObject::read(dict) && logFiles::read(dict))
    {
        logFiles::resetNames(dict.lookup("clouds"));

        Info<< type() << " " << name() << ": ";
        if (writeToFile() && names().size())
        {
            Info<< "applying to clouds:" << nl;
            forAll(names(), i)
            {
                Info<< "    " << names()[i] << nl;
                writeFileHeader(files(i));
            }
            Info<< endl;
        }
        else
        {
            Info<< "no clouds to be processed" << nl << endl;
        }

        return true;
    }

    return true;
}


bool Foam::functionObjects::cloudInfo::execute()
{
    return true;
}


bool Foam::functionObjects::cloudInfo::write()
{
    forAll(names(), i)
    {
        const word& cloudName = names()[i];

        const kinematicCloud& cloud =
            obr_.lookupObject<kinematicCloud>(cloudName);

        label nParcels = returnReduce(cloud.nParcels(), sumOp<label>());
        scalar massInSystem =
            returnReduce(cloud.massInSystem(), sumOp<scalar>());

        scalar Dmax = cloud.Dmax();
        scalar D10 = cloud.Dij(1, 0);
        scalar D32 = cloud.Dij(3, 2);

        Log << type() << " " << name() <<  " write:" << nl
            << "    number of parcels : " << nParcels << nl
            << "    mass in system    : " << massInSystem << nl
            << "    maximum diameter  : " << Dmax << nl
            << "    D10 diameter      : " << D10 << nl
            << "    D32 diameter      : " << D32 << nl
            << endl;

        if (writeToFile())
        {
            writeTime(files(i));
            files(i)
                << token::TAB
                << nParcels << token::TAB
                << massInSystem << token::TAB
                << Dmax << token::TAB
                << D10 << token::TAB
                << D32 << token::TAB
                << endl;
        }
    }

    return true;
}


// ************************************************************************* //

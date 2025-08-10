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
    (c) 2013-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "writeDictionary/writeDictionary.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeDictionary, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeDictionary,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::writeDictionary::tryDirectory
(
    const label dicti,
    const word& location,
    bool& firstDict
)
{
    IOobject dictIO
    (
        dictNames_[dicti],
        location,
        obr_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (dictIO.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(dictIO);

        if (dict.digest() != digests_[dicti])
        {
            if (firstDict)
            {
                Info<< type() << " " << name() << " write:" << nl << endl;

                IOobject::writeDivider(Info);
                Info<< endl;
                firstDict = false;
            }

            Info<< dict.dictName() << dict << nl;

            IOobject::writeDivider(Info);

            digests_[dicti] = dict.digest();
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeDictionary::writeDictionary
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    dictNames_(),
    digests_()
{
    read(dict);
    execute();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeDictionary::~writeDictionary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeDictionary::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);

    wordList dictNames(dict.lookup("dictNames"));
    HashSet<word> uniqueNames(dictNames);
    dictNames_ = uniqueNames.toc();

    digests_.setSize(dictNames_.size(), SHA1Digest());

    Info<< type() << " " << name() << ": monitoring dictionaries:" << nl;
    if (dictNames_.size())
    {
        forAll(dictNames_, i)
        {
            Info<< "    " << dictNames_[i] << endl;
        }
    }
    else
    {
        Info<< "    none" << nl;
    }
    Info<< endl;

    return true;
}


bool Foam::functionObjects::writeDictionary::execute()
{
    return true;
}


bool Foam::functionObjects::writeDictionary::write()
{
    bool firstDict = true;
    forAll(dictNames_, i)
    {
        if (obr_.foundObject<dictionary>(dictNames_[i]))
        {
            const dictionary& dict =
                obr_.lookupObject<dictionary>(dictNames_[i]);

            if (dict.digest() != digests_[i])
            {
                if (firstDict)
                {
                    Info<< type() << " " << name() << " write:" << nl << endl;

                    IOobject::writeDivider(Info);
                    Info<< endl;
                    firstDict = false;
                }

                digests_[i] = dict.digest();

                Info<< dict.dictName() << dict << nl;
                IOobject::writeDivider(Info);
                Info<< endl;
            }
        }
        else
        {
            bool processed = tryDirectory(i, obr_.time().timeName(), firstDict);

            if (!processed)
            {
                processed = tryDirectory(i, obr_.time().constant(), firstDict);
            }

            if (!processed)
            {
                processed = tryDirectory(i, obr_.time().system(), firstDict);
            }

            if (!processed)
            {
                Info<< "    Unable to locate dictionary " << dictNames_[i]
                    << nl << endl;
            }
            else
            {
                Info<< endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //

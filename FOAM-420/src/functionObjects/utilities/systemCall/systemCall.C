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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "systemCall/systemCall.H"
#include "db/Time/Time.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCode.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(systemCall, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        systemCall,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::systemCall::systemCall
(
    const word& name,
    const Time&,
    const dictionary& dict
)
:
    functionObject(name),
    executeCalls_(),
    endCalls_(),
    writeCalls_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::systemCall::~systemCall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::systemCall::read(const dictionary& dict)
{
    functionObject::read(dict);

    dict.readIfPresent("executeCalls", executeCalls_);
    dict.readIfPresent("endCalls", endCalls_);
    dict.readIfPresent("writeCalls", writeCalls_);
    masterOnly_ = dict.lookupOrDefault<Switch>("masterOnly", false);

    if (executeCalls_.empty() && endCalls_.empty() && writeCalls_.empty())
    {
        WarningInFunction
            << "no executeCalls, endCalls or writeCalls defined."
            << endl;
    }
    else if (!dynamicCode::allowSystemOperations)
    {
        FatalErrorInFunction
            << "Executing user-supplied system calls is not enabled by "
            << "default because of " << nl
            << "security issues.  If you trust the case you can enable this "
            << "facility by " << nl
            << "adding to the InfoSwitches setting in the system controlDict:"
            << nl << nl
            << "    allowSystemOperations 1" << nl << nl
            << "The system controlDict is either" << nl << nl
            << "    ~/.FOAMcore/$FOAM_PROJECT_VERSION/controlDict" << nl << nl
            << "or" << nl << nl
            << "    $FOAM_PROJECT_DIR/etc/controlDict" << nl << nl
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::systemCall::execute()
{
    forAll(executeCalls_, calli)
    {
	    if (masterOnly_)
	    {
		    if (Pstream::master())
		    {
			    Foam::system(executeCalls_[calli]);
		    }
	    }
	    else
	    {
		    Foam::system(executeCalls_[calli]);
	    }
    }

    return true;
}


bool Foam::functionObjects::systemCall::end()
{
    forAll(endCalls_, calli)
    {
        Foam::system(endCalls_[calli]);
    }

    return true;
}


bool Foam::functionObjects::systemCall::write()
{
    forAll(writeCalls_, calli)
    {
        Foam::system(writeCalls_[calli]);
    }

    return true;
}


// ************************************************************************* //

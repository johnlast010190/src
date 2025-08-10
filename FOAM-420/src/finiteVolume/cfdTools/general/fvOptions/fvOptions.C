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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2019-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/fvOptions/fvOptions.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(options, 0);
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::fv::options::createIOobject
(
    const fvMesh& mesh
) const
{
    return createIOobject(mesh, mesh);
}

Foam::IOobject Foam::fv::options::createIOobject
(
    const fvMesh& mesh,
    const objectRegistry& obr
) const
{
    IOobject io
    (
        typeName,
        mesh.time().constant(),
        obr,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorInFunction
            << "Found " << io.name() << " file in constant folder. "
            << "Please move it to " << mesh.time().system()
            << " or remove the file! "
            << exit(FatalError);
        return io;
    }
    else
    {
        // Check if the fvOptions file is in system
        io.instance() = mesh.time().system();

        if (io.typeHeaderOk<IOdictionary>(true))
        {
            Info<< "Creating finite volume options from "
                << io.instance()/io.name() << nl
                << endl;

            io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
            return io;
        }
        else
        {
            Info<< "No finite volume options present" << nl << endl;

            io.readOpt() = IOobject::NO_READ;
            return io;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::options::options
(
    const fvMesh& mesh,
    const dictionary& injectDict,
    bool init
)
:
    IOdictionary(createIOobject(mesh)),
    optionList
    (
        mesh,
        (
            // Non-recursive, no-overwrite merge
            *this |= injectDict,
            *this
        )
    ),
    injectDict_(injectDict)
{
    if (init)
    {
        this->initialise();
    }
}

Foam::fv::options::options
(
    const fvMesh& mesh,
    const objectRegistry& obr,
    const dictionary& injectDict,
    bool init
)
:
    IOdictionary(createIOobject(mesh, obr)),
    optionList
    (
        obr,
        (
            // Non-recursive, no-overwrite merge
            *this |= injectDict,
            *this
        )
    ),
    injectDict_(injectDict)
{
    if (init)
    {
        this->initialise();
    }
}

Foam::fv::options& Foam::fv::options::New
(
    const fvMesh& mesh,
    const dictionary& injectDict,
    bool init
)
{
    if (mesh.thisDb().foundObject<options>(typeName))
    {
        return const_cast<options&>
        (
            mesh.lookupObject<options>(typeName)
        );
    }
    else
    {
        if (debug)
        {
            InfoInFunction
                << "Constructing " << typeName
                << " for region " << mesh.name() << endl;
        }

        options* objectPtr = new options(mesh, injectDict, init);
        regIOobject::store(objectPtr);
        return *objectPtr;
    }
}

Foam::fv::options& Foam::fv::options::New
(
    const fvMesh& mesh,
    const objectRegistry& obr,
    const dictionary& injectDict,
    bool init
)
{
    if (obr.foundObject<options>(typeName))
    {
        return const_cast<options&>
        (
            obr.lookupObject<options>(typeName)
        );
    }
    else
    {
        if (debug)
        {
            InfoInFunction
                << "Constructing " << typeName
                << " for region " << mesh.name() << endl;
        }

        options* objectPtr = new options(mesh, obr, injectDict, init);
        regIOobject::store(objectPtr);
        return *objectPtr;
    }
}

bool Foam::fv::options::read()
{
    if (IOdictionary::regIOobject::read())
    {
        // Non-recursive, no-overwrite merge
        *this |= injectDict_;
        optionList::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

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
    (c) 2012-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/porosityModel/porosityModel/IOporosityModelList.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::IOporosityModelList::createIOobject
(
    const fvMesh& mesh,
    bool readFromFvOptions
) const
{

    if (readFromFvOptions)
    {
        IOobject io
        (
            "fvOptions",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (io.typeHeaderOk<IOdictionary>(true))
        {
            Info<< "Creating porosity model list from " << io.name() << nl << endl;

            io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
            return io;
        }
        else
        {
            Info<< "No porosity models present" << nl << endl;

            io.readOpt() = IOobject::NO_READ;
            return io;
        }
    }
    else
    {
        IOobject io
        (
            "porosityProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (io.typeHeaderOk<IOdictionary>(true))
        {
            Info<< "Creating porosity model list from " << io.name() << nl << endl;

            io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
            return io;
        }
        else
        {
            Info<< "No porosity models present" << nl << endl;

            io.readOpt() = IOobject::NO_READ;
            return io;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOporosityModelList::IOporosityModelList
(
    const fvMesh& mesh,
    bool readFromFvOptions
)
:
    IOdictionary(createIOobject(mesh, readFromFvOptions)),
    porosityModelList(mesh, mesh, *this, readFromFvOptions)
{}

Foam::IOporosityModelList::IOporosityModelList
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    bool readFromFvOptions
)
:
    IOdictionary(createIOobject(mesh, readFromFvOptions)),
    porosityModelList(obr, mesh, *this, readFromFvOptions)
{}


bool Foam::IOporosityModelList::read()
{
    if (regIOobject::read())
    {
        porosityModelList::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

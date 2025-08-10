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

#include "engineTime/engineTime.H"
#include "ignition/ignition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ignition::ignition
(
    const dictionary& combustionProperties,
    const Time& db,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    ignite_(combustionProperties.lookup("ignite")),
    ignSites_
    (
        combustionProperties.lookup("ignitionSites"),
        ignitionSite::iNew(db, mesh)
    )
{
    if (ignite_)
    {
        Info<< "\nIgnition on" << endl;
    }
    else
    {
        Info<< "\nIgnition switched off" << endl;
    }
}


Foam::ignition::ignition
(
    const dictionary& combustionProperties,
    const engineTime& edb,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    ignite_(combustionProperties.lookup("ignite")),
    ignSites_
    (
        combustionProperties.lookup("ignitionSites"),
        ignitionSite::iNew(edb, mesh)
    )
{
    if (ignite_)
    {
        Info<< "\nIgnition on" << endl;
    }
    else
    {
        Info<< "\nIgnition switched off" << endl;
    }
}


// ************************************************************************* //

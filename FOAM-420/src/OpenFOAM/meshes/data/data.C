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

\*---------------------------------------------------------------------------*/

#include "meshes/data/data.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::data::debug(Foam::debug::debugSwitch("data", false));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::data::data(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "data",
            obr.time().system(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    prevTimeIndices_(),
    prevTimeBlockIndices_()
{
    set("solverPerformance", dictionary());
    set("blockSolverPerformance", dictionary());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::data::solverPerformanceDict() const
{
    return subDict("solverPerformance");
}

const Foam::dictionary& Foam::data::blockSolverPerformanceDict() const
{
    return subDict("blockSolverPerformance");
}


// ************************************************************************* //

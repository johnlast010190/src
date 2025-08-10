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
    (c) 2018-2022 OpenCFD Ltd.
\*---------------------------------------------------------------------------*/

#include "db/Time/instant/instant.H"
#include "db/Time/Time.H"
#include <cstdlib>  // std::atof
#include <utility>  // std::move

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::instant::typeName = "instant";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::instant::instant(scalar timeValue)
:
    Instant<word>(timeValue, Time::timeName(timeValue))
{}


Foam::instant::instant(const word& timeName)
:
    Instant<word>(0, timeName)
{
    value() = std::atof(name().c_str());
}


Foam::instant::instant(word&& timeName)
:
    Instant<word>(0, std::move(timeName))
{
    value() = std::atof(name().c_str());
}


// ************************************************************************* //

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

#include "FSD/reactionRateFlameAreaModels/reactionRateFlameArea/reactionRateFlameArea.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::reactionRateFlameArea> Foam::reactionRateFlameArea::New
(
    const dictionary& dict,
    const objectRegistry& obr,
    const combustionModel& combModel
)
{
    word reactionRateFlameAreaType
    (
        dict.lookup("reactionRateFlameArea")
    );

    Info<< "Selecting reaction rate flame area correlation "
        << reactionRateFlameAreaType << endl;

    const auto ctor = ctorTableLookup("reactionRateFlameArea type", dictionaryConstructorTable_(), reactionRateFlameAreaType);

    const label tempOpen = reactionRateFlameAreaType.find('<');

    const word className = reactionRateFlameAreaType(0, tempOpen);

    return autoPtr<reactionRateFlameArea>
        (ctor(className, dict, obr, combModel));
}


// ************************************************************************* //

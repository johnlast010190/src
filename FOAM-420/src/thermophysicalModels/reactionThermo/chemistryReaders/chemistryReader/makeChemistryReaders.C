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
    (c) 2017-2022 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "makeReactionThermo.H"
#include "include/thermoPhysicsTypes.H"
#include "include/dummyThermo.H"

#include "include/solidThermoPhysicsTypes.H"

#include "chemistryReaders/chemistryReader/chemistryReader.H"
#include "chemistryReaders/foamChemistryReader/foamChemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Solid chemistry readers based on sensibleEnthalpy

makeChemistryReader(constGasHThermoPhysics);
makeChemistryReader(gasHThermoPhysics);
makeChemistryReader(constIncompressibleGasHThermoPhysics);
makeChemistryReader(incompressibleGasHThermoPhysics);
makeChemistryReader(constIcoFluidHThermoPhysics);
makeChemistryReader(icoPoly8HThermoPhysics);
makeChemistryReader(tabulatedGasHThermoPhysics);

makeChemistryReaderType(foamChemistryReader, constGasHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, gasHThermoPhysics);
makeChemistryReaderType
(
    foamChemistryReader,
    constIncompressibleGasHThermoPhysics
);
makeChemistryReaderType(foamChemistryReader, incompressibleGasHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, constIcoFluidHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, icoPoly8HThermoPhysics);
makeChemistryReaderType(foamChemistryReader, tabulatedGasHThermoPhysics);

// Solid chemistry readers based on sensibleInternalEnergy

makeChemistryReader(constGasEThermoPhysics);
makeChemistryReader(gasEThermoPhysics);
makeChemistryReader(constIncompressibleGasEThermoPhysics);
makeChemistryReader(incompressibleGasEThermoPhysics);
makeChemistryReader(constIcoFluidEThermoPhysics);
makeChemistryReader(icoPoly8EThermoPhysics);
makeChemistryReader(tabulatedGasEThermoPhysics);

makeChemistryReaderType(foamChemistryReader, constGasEThermoPhysics);
makeChemistryReaderType(foamChemistryReader, gasEThermoPhysics);
makeChemistryReaderType
(
    foamChemistryReader,
    constIncompressibleGasEThermoPhysics
);
makeChemistryReaderType(foamChemistryReader, incompressibleGasEThermoPhysics);
makeChemistryReaderType(foamChemistryReader, constIcoFluidEThermoPhysics);
makeChemistryReaderType(foamChemistryReader, icoPoly8EThermoPhysics);
makeChemistryReaderType(foamChemistryReader, tabulatedGasEThermoPhysics);

// Solid chemistry readers for solids based on sensibleInternalEnergy

makeChemistryReader(hConstSolidThermoPhysics);
makeChemistryReader(hPowerSolidThermoPhysics);
makeChemistryReader(hExpKappaConstSolidThermoPhysics);

makeChemistryReaderType(foamChemistryReader, hConstSolidThermoPhysics);
makeChemistryReaderType(foamChemistryReader, hPowerSolidThermoPhysics);
makeChemistryReaderType(foamChemistryReader, hExpKappaConstSolidThermoPhysics);


// Material models
makeChemistryReader(dummyThermo);
makeChemistryReaderType(foamChemistryReader, dummyThermo);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

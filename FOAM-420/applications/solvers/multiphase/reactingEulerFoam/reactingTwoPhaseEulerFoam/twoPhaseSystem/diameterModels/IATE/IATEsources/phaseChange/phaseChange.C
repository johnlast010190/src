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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "phaseChange.H"
#include "twoPhaseSystem.H"
#include "phaseSystem/phaseSystem.H"
#include "PhaseSystems/ThermalPhaseChangePhaseSystem/ThermalPhaseChangePhaseSystem.H"
#include "PhaseSystems/MomentumTransferPhaseSystem/MomentumTransferPhaseSystem.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{
    defineTypeNameAndDebug(phaseChange, 0);
    addToRunTimeSelectionTable(IATEsource, phaseChange, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATEsources::phaseChange::phaseChange
(
    const IATE& iate,
    const dictionary& dict
)
:
    IATEsource(iate)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::IATEsources::phaseChange::R
(
    const volScalarField& alphai,
    volScalarField& kappai
) const
{
    const ThermalPhaseChangePhaseSystem
    <
        MomentumTransferPhaseSystem
        <
            twoPhaseSystem
        >
    >& phaseChangeFluid = refCast
    <
        const ThermalPhaseChangePhaseSystem
        <
            MomentumTransferPhaseSystem<twoPhaseSystem>
        >
    >(fluid());

    return -fvm::SuSp
    (
        (1.0/3.0)
       *phaseChangeFluid.iDmdt(phase())()()
       /(alphai()*phase().rho()()),
        kappai
    );
}


// ************************************************************************* //

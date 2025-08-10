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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "SolverFunction.H"

#include "turbulence/turbulenceSolver.H"
#include "concentrationTransport/concentrationTransport.H"
#include "disperseEulerian/disperseEulerian.H"
#include "nonParticipatingRadiation/nonParticipatingRadiation.H"
#include "passiveScalarTransport/passiveScalarTransport.H"
#include "scalarTransport/scalarTransport.H"
#include "windDrivenRain/windDrivenRain.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeSolverFunction(turbulenceSolver);
makeSolverFunction(concentrationTransport);
makeSolverFunction(disperseEulerian);
makeSolverFunction(nonParticipatingRadiation);
makeSolverFunction(passiveScalarTransport);
makeSolverFunction(scalarTransport);
makeSolverFunction(windDrivenRain);

// ************************************************************************* //

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

#include "clouds/derived/basicReactingCloud/basicReactingCloud.H"

#include "parcels/include/makeParcelCloudFunctionObjects.H"

// Kinematic
#include "parcels/include/makeThermoParcelForces.H" // thermo variant
#include "parcels/include/makeParcelDispersionModels.H"
#include "parcels/include/makeReactingParcelInjectionModels.H" // Reacting variant
#include "parcels/include/makeParcelPatchInteractionModels.H"
#include "parcels/include/makeParcelStochasticCollisionModels.H"
#include "parcels/include/makeReactingParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "parcels/include/makeParcelHeatTransferModels.H"

// Reacting
#include "parcels/include/makeReactingParcelCompositionModels.H"
#include "parcels/include/makeReactingParcelPhaseChangeModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(basicReactingCloud);

// Kinematic sub-models
makeThermoParcelForces(basicReactingCloud);
makeParcelDispersionModels(basicReactingCloud);
makeReactingParcelInjectionModels(basicReactingCloud);
makeParcelPatchInteractionModels(basicReactingCloud);
makeParcelStochasticCollisionModels(basicReactingCloud);
makeReactingParcelSurfaceFilmModels(basicReactingCloud);

// Thermo sub-models
makeParcelHeatTransferModels(basicReactingCloud);

// Reacting sub-models
makeReactingParcelCompositionModels(basicReactingCloud);
makeReactingParcelPhaseChangeModels(basicReactingCloud);


// ************************************************************************* //

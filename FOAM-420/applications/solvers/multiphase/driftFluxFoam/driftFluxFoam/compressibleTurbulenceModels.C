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
    (c) 2014-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "CompressibleTurbulenceModel/CompressibleTurbulenceModel.H"
#include "incompressibleTwoPhaseInteractingMixture.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "laminar/laminarModel/laminarModel.H"
#include "RAS/RASModel/RASModel.H"
#include "LES/LESModel/LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeTurbulenceModelTypes
(
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    CompressibleTurbulenceModel,
    incompressibleTwoPhaseInteractingMixture
);

makeBaseTurbulenceModel
(
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    CompressibleTurbulenceModel,
    incompressibleTwoPhaseInteractingMixture
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedTurbulenceModel                                               \
    (                                                                          \
        incompressibleTwoPhaseInteractingMixtureCompressibleTurbulenceModel,   \
        laminar,                                                               \
        Type                                                                   \
    )

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (                                                                          \
        incompressibleTwoPhaseInteractingMixtureCompressibleTurbulenceModel,   \
        RAS,                                                                   \
        Type                                                                   \
    )

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (                                                                          \
        incompressibleTwoPhaseInteractingMixtureCompressibleTurbulenceModel,   \
        LES,                                                                   \
        Type                                                                   \
    )

#include "laminar/Stokes/Stokes.H"
makeLaminarModel(Stokes);

#include "RAS/kEpsilon/kEpsilon.H"
makeRASModel(kEpsilon);

#include "RAS/buoyantKEpsilon/buoyantKEpsilon.H"
makeRASModel(buoyantKEpsilon);

#include "LES/Smagorinsky/Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "LES/kEqn/kEqn.H"
makeLESModel(kEqn);


// ************************************************************************* //

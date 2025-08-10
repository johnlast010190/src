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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#ifndef makeReactingParcelInjectionModels_H
#define makeReactingParcelInjectionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Kinematic/InjectionModel/CellZoneInjection/CellZoneInjection.H"
#include "submodels/Kinematic/InjectionModel/ConeInjection/ConeInjection.H"
#include "submodels/Kinematic/InjectionModel/ConeNozzleInjection/ConeNozzleInjection.H"
#include "submodels/Kinematic/InjectionModel/FieldActivatedInjection/FieldActivatedInjection.H"
#include "submodels/Kinematic/InjectionModel/ManualInjection/ManualInjection.H"
#include "submodels/Kinematic/InjectionModel/NoInjection/NoInjection.H"
#include "submodels/Kinematic/InjectionModel/PatchInjection/PatchInjection.H"
#include "submodels/Kinematic/InjectionModel/PatchFlowRateInjection/PatchFlowRateInjection.H"
#include "submodels/Reacting/InjectionModel/ReactingLookupTableInjection/ReactingLookupTableInjection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingParcelInjectionModels(CloudType)                           \
                                                                               \
    makeInjectionModel(CloudType);                                             \
    makeInjectionModelType(CellZoneInjection, CloudType);                      \
    makeInjectionModelType(ConeInjection, CloudType);                          \
    makeInjectionModelType(ConeNozzleInjection, CloudType);                    \
    makeInjectionModelType(FieldActivatedInjection, CloudType);                \
    makeInjectionModelType(ManualInjection, CloudType);                        \
    makeInjectionModelType(NoInjection, CloudType);                            \
    makeInjectionModelType(PatchInjection, CloudType);                         \
    makeInjectionModelType(PatchFlowRateInjection, CloudType);                 \
    makeInjectionModelType(ReactingLookupTableInjection, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

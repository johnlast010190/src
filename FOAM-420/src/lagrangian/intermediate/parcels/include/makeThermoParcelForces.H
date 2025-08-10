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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#ifndef makeThermoParcelForces_H
#define makeThermoParcelForces_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Kinematic/ParticleForces/Drag/SphereDrag/SphereDragForce.H"
#include "submodels/Kinematic/ParticleForces/Drag/NonSphereDrag/NonSphereDragForce.H"

#include "submodels/Kinematic/ParticleForces/Lift/SaffmanMeiLift/SaffmanMeiLiftForce.H"
#include "submodels/Kinematic/ParticleForces/Lift/TomiyamaLift/TomiyamaLiftForce.H"

#include "submodels/Kinematic/ParticleForces/Gravity/GravityForce.H"
#include "submodels/Kinematic/ParticleForces/NonInertialFrame/NonInertialFrameForce.H"
#include "submodels/Kinematic/ParticleForces/Paramagnetic/ParamagneticForce.H"
#include "submodels/Kinematic/ParticleForces/PressureGradient/PressureGradientForce.H"
#include "submodels/Kinematic/ParticleForces/SRF/SRFForce.H"
#include "submodels/Kinematic/ParticleForces/VirtualMass/VirtualMassForce.H"
#include "submodels/Kinematic/ParticleForces/ElectricField/ElectricFieldForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermoParcelForces(CloudType)                                      \
                                                                               \
    makeParticleForceModel(CloudType);                                         \
    makeParticleForceModelType(SphereDragForce, CloudType);                    \
    makeParticleForceModelType(NonSphereDragForce, CloudType);                 \
    makeParticleForceModelType(SaffmanMeiLiftForce, CloudType);                \
    makeParticleForceModelType(TomiyamaLiftForce, CloudType);                  \
    makeParticleForceModelType(GravityForce, CloudType);                       \
    makeParticleForceModelType(NonInertialFrameForce, CloudType);              \
    makeParticleForceModelType(ParamagneticForce, CloudType);                  \
    makeParticleForceModelType(PressureGradientForce, CloudType);              \
    makeParticleForceModelType(SRFForce, CloudType);                           \
    makeParticleForceModelType(VirtualMassForce, CloudType);                   \
    makeParticleForceModelType(ElectricFieldForce, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

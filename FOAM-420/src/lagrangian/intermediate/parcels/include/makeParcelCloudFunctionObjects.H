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
    (c) 2010-2015 Esi Ltd.

\*---------------------------------------------------------------------------*/

#ifndef makeParcelCloudFunctionObjects_H
#define makeParcelCloudFunctionObjects_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/CloudFunctionObjects/CloudToVTK/CloudToVTK.H"
#include "submodels/CloudFunctionObjects/FacePostProcessing/FacePostProcessing.H"
#include "submodels/CloudFunctionObjects/ParticleCollector/ParticleCollector.H"
#include "submodels/CloudFunctionObjects/ParticleErosion/ParticleErosion.H"
#include "submodels/CloudFunctionObjects/ParticleTracks/ParticleTracks.H"
#include "submodels/CloudFunctionObjects/ParticleTrap/ParticleTrap.H"
#include "submodels/CloudFunctionObjects/PatchPostProcessing/PatchPostProcessing.H"
#include "submodels/CloudFunctionObjects/VoidFraction/VoidFraction.H"
#include "submodels/CloudFunctionObjects/areaAveragedForce/areaAveragedForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelCloudFunctionObjects(CloudType)                              \
                                                                               \
    makeCloudFunctionObject(CloudType);                                        \
                                                                               \
    makeCloudFunctionObjectType(FacePostProcessing, CloudType);                \
    makeCloudFunctionObjectType(ParticleCollector, CloudType);                 \
    makeCloudFunctionObjectType(ParticleErosion, CloudType);                   \
    makeCloudFunctionObjectType(ParticleTracks, CloudType);                    \
    makeCloudFunctionObjectType(ParticleTrap, CloudType);                      \
    makeCloudFunctionObjectType(PatchPostProcessing, CloudType);               \
    makeCloudFunctionObjectType(VoidFraction, CloudType);                      \
    makeCloudFunctionObjectType(areaAveragedForce, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

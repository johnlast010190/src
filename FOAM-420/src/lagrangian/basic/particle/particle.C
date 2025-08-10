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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "particle/particle.H"
#include "primitives/transform/transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::particle::particleCount_ = 0;

const Foam::scalar Foam::particle::trackingCorrectionTol = 1e-5;

const Foam::scalar Foam::particle::lambdaDistanceToleranceCoeff = 1e3*SMALL;

const Foam::scalar Foam::particle::minStepFractionTol = 1e5*SMALL;

namespace Foam
{
    defineTypeNameAndDebug(particle, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particle::particle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    mesh_(mesh),
    position_(position),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(tetFacei),
    tetPti_(tetPti),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID()),
    noBasePt_(false)
{}


Foam::particle::particle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    bool doCellFacePt
)
:
    mesh_(mesh),
    position_(position),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(-1),
    tetPti_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID()),
    noBasePt_(false)
{
    if (doCellFacePt)
    {
        initCellFacePt();
    }
}


Foam::particle::particle(const particle& p)
:
    mesh_(p.mesh_),
    position_(p.position_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    origProc_(p.origProc_),
    origId_(p.origId_),
    noBasePt_(p.noBasePt_)
{}


Foam::particle::particle(const particle& p, const polyMesh& mesh)
:
    mesh_(mesh),
    position_(p.position_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    origProc_(p.origProc_),
    origId_(p.origId_),
    noBasePt_(p.noBasePt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::particle::transformProperties(const transformer&)
{}


Foam::scalar Foam::particle::wallImpactDistance(const vector&) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

bool Foam::operator==(const particle& pA, const particle& pB)
{
    return (pA.origProc() == pB.origProc() && pA.origId() == pB.origId());
}


bool Foam::operator!=(const particle& pA, const particle& pB)
{
    return !(pA == pB);
}


// ************************************************************************* //

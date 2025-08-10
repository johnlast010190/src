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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/axisRotationMotion/axisRotationMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/unitConversion/unitConversion.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(axisRotationMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        axisRotationMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        axisRotationMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::axisRotationMotion::axisRotationMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    origin_(SBMFCoeffs_.lookup("origin"))
{
    axisRotationMotion::read(SBMFCoeffs);
}


Foam::solidBodyMotionFunctions::axisRotationMotion::axisRotationMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    origin_(CofR0())
{
    axisRotationMotion::read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::axisRotationMotion::~axisRotationMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::vectorTuple
Foam::solidBodyMotionFunctions::axisRotationMotion::calcInfo() const
{
    const vector omega
    (
        degToRad(radialVelocity_.x()),
        degToRad(radialVelocity_.y()),
        degToRad(radialVelocity_.z())
    );

    if (mag(omega) > VSMALL)
    {
        vector axis(normalised(omega));

        if (hasFrame() && isIncrementalMotion())
        {
            axis = frame().coorSys().transform(axis);
            origin_ = frame().coorSys().origin();
        }

        return vectorTuple(omega, axis);
    }

    return vectorTuple(Zero, Zero);
}


Foam::septernion
Foam::solidBodyMotionFunctions::axisRotationMotion::transformation
(
    const scalar t0,
    const scalar t
) const
{
    const scalar dt = t - t0;
    const vectorTuple omegaAxis(calcInfo());
    const vector angles(omegaAxis.first()*dt);
    scalar angle = mag(angles);

    septernion TR = septernion::I;

    if (angle > VSMALL)
    {
        const quaternion R(omegaAxis.second(), angle);
        TR = septernion(-origin_)*R*septernion(origin_);
    }

    DebugInFunction << "Time = " << t << ", t0 = " << t0
                    << ", transformation: " << TR << endl;

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::axisRotationMotion::velocity() const
{
    const vectorTuple omegaAxis(calcInfo());
    return vectorTuple(Zero, omegaAxis.second()*mag(omegaAxis.first()));
}


bool Foam::solidBodyMotionFunctions::axisRotationMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);
    radialVelocity_ = SBMFCoeffs_.lookup<vector>("radialVelocity");

    return true;
}


// ************************************************************************* //

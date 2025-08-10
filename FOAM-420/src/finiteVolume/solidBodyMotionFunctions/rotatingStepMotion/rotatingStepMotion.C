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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/rotatingStepMotion/rotatingStepMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/unitConversion/unitConversion.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvMesh.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(rotatingStepMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingStepMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingStepMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingStepMotion::rotatingStepMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            SBMFCoeffs_.lookupOrDefault<word>("region", "region0")
        )
    ),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis")),
    theta_(SBMFCoeffs_.lookup<scalar>("theta")),
    theta0_(SBMFCoeffs_.lookupOrDefault<scalar>("theta0", 0.0)),
    stepPeriod_(readLabel(SBMFCoeffs_.lookup("period"))),
    stO_(septernion::zero),
    UName_(SBMFCoeffs_.lookupOrDefault<word>("U", "U")),
    lastUpdate_(-1)
{}


Foam::solidBodyMotionFunctions::rotatingStepMotion::rotatingStepMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    obr_(obr),
    origin_(CofR0()),
    axis_(frame().axis0()),
    theta_(SBMFCoeffs_.lookup<scalar>("theta")),
    theta0_(SBMFCoeffs_.lookupOrDefault<scalar>("theta0", 0.0)),
    stepPeriod_(readLabel(SBMFCoeffs_.lookup("period"))),
    stO_(septernion::zero),
    UName_(SBMFCoeffs_.lookupOrDefault<word>("U", "U")),
    lastUpdate_(-1)
{
    if (isIncrementalMotion())
    {
        FatalErrorInFunction
            << "The rotatingStepMotion function "
            << "does not support incremental motion"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingStepMotion::~rotatingStepMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::rotatingStepMotion::transformation
(
    const scalar t0,
    const scalar t1
) const
{
    label t = time_.timeIndex() - 1;
    //-1 added so mesh changes after end of period

    label nSteps = t / stepPeriod_;
    label rem = t % stepPeriod_;

    // Stepwise rotation around axis
    scalar angle = degToRad(theta0_ + nSteps*theta_);

    quaternion R(axis_, angle);

    stO_ = (septernion(-origin_)*R*septernion(origin_));

    if (rem == 0 && lastUpdate_ != t && nSteps > 0)
    {
        lastUpdate_ = t;

        Info<< "solidBodyMotionFunctions::rotatingStepMotion"
            << "::transformation(): Time = " << t
            << " transformation: " << stO_ << endl;

    }

    return stO_;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::rotatingStepMotion::velocity() const
{
    // Dummy steady-state motion, so return zero velocity!
    return vectorTuple(Zero, Zero);
}


bool Foam::solidBodyMotionFunctions::rotatingStepMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    return true;
}

// ************************************************************************* //

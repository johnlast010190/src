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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/rotatingMotion/rotatingMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(rotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingMotion::rotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis")),
    omega_(Function1<scalar>::New("omega", SBMFCoeffs_))
{}


Foam::solidBodyMotionFunctions::rotatingMotion::rotatingMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    origin_(CofR0()),
    axis_(frame().axis0()/mag(frame().axis0())),
    omega_(Function1<scalar>::New("omega", SBMFCoeffs_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingMotion::~rotatingMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::solidBodyMotionFunctions::rotatingMotion::updateOriginAndAxis() const
{
    // Update origin and axis if relative motion
    if (hasFrame() && isIncrementalMotion())
    {
        origin_ = frame().coorSys().origin();
        axis_ = normalised(frame().axis());
    }
}


Foam::septernion
Foam::solidBodyMotionFunctions::rotatingMotion::transformation
(
    const scalar t1,
    const scalar t2
) const
{
    updateOriginAndAxis();

    // Rotation around axis
    const scalar angle = omega_->integrate(t1, t2);
    const quaternion R(axis_, angle);
    const septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction << "Time = " << t2 << ", t1 = " << t1
                    << ", transformation: " << TR << endl;

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::rotatingMotion::velocity() const
{
    return vectorTuple(Zero, axis_*omega_->value(time_.value()));
}


bool Foam::solidBodyMotionFunctions::rotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);
    omega_.reset(Function1<scalar>::New("omega", SBMFCoeffs_).ptr());
    return true;
}


// ************************************************************************* //

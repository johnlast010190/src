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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/linearMotion/linearMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(linearMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        linearMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        linearMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearMotion::linearMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    velocity_(nullptr)
{
    linearMotion::read(SBMFCoeffs);
}


Foam::solidBodyMotionFunctions::linearMotion::linearMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    velocity_(nullptr)
{
    linearMotion::read(SBMFCoeffs);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearMotion::~linearMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::linearMotion::transformation
(
    const scalar t0,
    const scalar t
) const
{
    // Translation of centre of gravity with constant velocity
    vector displacement = velocity_->integrate(t0, t);

    if (hasFrame())
    {
        if (isIncrementalMotion() && frame().coorSys().uniform())
        {
            displacement = frame().coorSys().transform(displacement);
        }
        else if (frame().coorSys0().uniform())
        {
            displacement = frame().coorSys0().transform(displacement);
        }
        else
        {
            FatalErrorInFunction
                << "linearMotion allow only uniform coordinate systems"
                << exit(FatalError);
        }
    }

    quaternion R(1);
    septernion TR(septernion(-displacement)*R);

    DebugInFunction << "Time = " << t << ", t0 = " << t0
                    << ", transformation: " << TR << endl;

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::linearMotion::velocity() const
{
    vector vel(velocity_->value(time_.value()));

    if (hasFrame())
    {
        // Change direction of velocity vector by nested frames
        if (isIncrementalMotion() && frame().coorSys().uniform())
        {
            vel = frame().coorSys().transform(vel);
        }
        if (frame().coorSys0().uniform())
        {
            vel = frame().coorSys0().transform(vel);
        }
        else
        {
            FatalErrorInFunction
                << "linearMotion allow only uniform coordinate systems"
                << exit(FatalError);
        }
    }

    return vectorTuple(vel, Zero);
}


bool Foam::solidBodyMotionFunctions::linearMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    velocity_.reset
    (
        Function1<vector>::New("velocity", SBMFCoeffs_).ptr()
    );

    return true;
}


// ************************************************************************* //

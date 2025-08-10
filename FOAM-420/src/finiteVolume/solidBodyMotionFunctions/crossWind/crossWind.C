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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/crossWind/crossWind.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(crossWind, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        crossWind,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::crossWind::crossWind
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    velocity_(nullptr),
    frameForRotation_
    (
        coordinateFrame::New
        (
            dynamic_cast<const fvMesh&>(obr),
            SBMFCoeffs.lookup<word>("frameForRotation")
        )
    )
{
    crossWind::read(SBMFCoeffs);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::crossWind::~crossWind()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion Foam::solidBodyMotionFunctions::crossWind::transformation
(
    const scalar t0,
    const scalar t
) const
{
    frameForRotation_.updateState();

    septernion TR(septernion::I);

    forAll(frameForRotation_.parents(), i)
    {
        if (!frameForRotation_.parents()[i].isDynamic())
        {
            TR *=
                inv(frameForRotation_.parents()[i].decoupledTransformation());
        }
    }

    DebugInFunction << "Time = " << t << ", t0 = " << t0
                    << ", transformation: " << TR << endl;

    return TR;
}


Foam::vectorTuple Foam::solidBodyMotionFunctions::crossWind::velocity() const
{
    vector vel(velocity_->value(time_.value()));

    // Change direction of velocity vector by nested frames
    return vectorTuple(frame().coorSys().transform(vel), Zero);
}


bool Foam::solidBodyMotionFunctions::crossWind::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    velocity_.reset(Function1<vector>::New("velocity", SBMFCoeffs_).ptr());

    if (!frame().coorSys().uniform())
    {
        FatalErrorInFunction
            << "crossWind allow only uniform coordinate systems"
            << exit(FatalError);
    }

    return true;
}


// ************************************************************************* //

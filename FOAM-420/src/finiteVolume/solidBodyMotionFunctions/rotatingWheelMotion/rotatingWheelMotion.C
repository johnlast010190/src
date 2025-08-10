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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/rotatingWheelMotion/rotatingWheelMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(rotatingWheelMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingWheelMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingWheelMotion::rotatingWheelMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    omega0_(),
    omega_(vector::zero),
    vehicleSpeed_()
{
    rotatingWheelMotion::read(SBMFCoeffs);
    setIncrementalMotion(rotatingWheelMotion::typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingWheelMotion::~rotatingWheelMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::rotatingWheelMotion::transformation() const
{
    updateTime();
    updateMotionCoeffs();

    // Approximation rotation around axis
    scalar angle = dt_*mag(omega_);

    if (SBMFCoeffs_.found("omega"))
    {
        angle = omega0_->integrate(tOld_, time_.value());
    }

    quaternion R(frame().axis(), angle);
    const vector& CofR = frame().CofR();
    septernion TR(septernion(-CofR)*R*septernion(CofR));

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::rotatingWheelMotion::velocity() const
{
    updateMotionCoeffs();
    return vectorTuple(Zero, omega_);
}


void
Foam::solidBodyMotionFunctions::rotatingWheelMotion::updateMotionCoeffs() const
{
    const scalar t = time_.value();
    const vector axis = frame().axis();
    const vector CofR = frame().CofR();

    if (!omega0_.valid())
    {
        if (frame().validParentFrame())
        {
            const coordinateFrame& parentFrame = frame().parentFrame();
            const vector parentAxis = parentFrame.axis();
            const vector parentCofR = parentFrame.CofR();

            vector Uwoc = parentFrame.frameVelocity(CofR);
            Uwoc -= axis*(axis & Uwoc);

            if (mag(Uwoc) == 0.0)
            {
                omega_ = 0.0;
            }
            else
            {
                scalar contactRadius;
                // Corrected wheel velocity account for toe-in and camber
                vector contactDir(Uwoc ^ axis);
                if (pos0(contactDir & parentAxis))
                {
                    contactDir = -contactDir;
                }
                contactDir.normalise();

                if (SBMFCoeffs_.found("contactRadius"))
                {
                    contactRadius =
                        SBMFCoeffs_.lookup<scalar>("contactRadius");
                }
                else
                {
                    const scalar zcd(mag((parentAxis & (parentCofR - CofR))));
                    contactRadius = zcd / (contactDir & (-parentAxis));
                }
                omega_ = (axis & ((Uwoc ^ contactDir)/contactRadius))*axis;
            }
        }
        else if (vehicleSpeed_.valid())
        {
            // for wind tunnel cases where vehicleFrame is not defined
            // assume wheel coordinate system is specified correctly:
            // 1) axis is set so omega is positive
            // 2) origin is either specified as distance from ground
            //    or contactRadius is given
            const scalar Uwoc = vehicleSpeed_->value(t);

            if (Uwoc == 0.0)
            {
                omega_ = 0.0;
            }
            else
            {
                scalar contactRadius;

                if (SBMFCoeffs_.found("contactRadius"))
                {
                    contactRadius =
                        SBMFCoeffs_.lookup<scalar>("contactRadius");
                }
                else
                {
                    contactRadius = mag(CofR.z());
                }
                omega_ = Uwoc/contactRadius * axis;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Could not calculate or find omega value for "
                << rotatingWheelMotion::typeName
                << exit(FatalError);
        }
    }
    else
    {
        omega_ = omega0().value(t)*axis;
    }
}


const Foam::Function1<Foam::scalar>&
Foam::solidBodyMotionFunctions::rotatingWheelMotion::omega0() const
{
    if (!omega0_.valid())
    {
        FatalErrorInFunction
            << "Could not calculate or find omega value for "
            << rotatingWheelMotion::typeName
            << exit(FatalError);
    }

    return omega0_();
}


bool Foam::solidBodyMotionFunctions::rotatingWheelMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    if (SBMFCoeffs_.found("omega"))
    {
        omega0_ = Function1<scalar>::New("omega", SBMFCoeffs_);
    }
    if (SBMFCoeffs_.found("vehicleSpeed"))
    {
        vehicleSpeed_ = Function1<scalar>::New("vehicleSpeed", SBMFCoeffs_);
    }

    return true;
}

// ************************************************************************* //

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

#include "solidBodyMotionFunctions/vehicleMotion/vehicleMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(vehicleMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        vehicleMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::vehicleMotion::vehicleMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    CofR_(Zero),
    vehicleSpeed_(),
    turnRadius_(),
    invTurnRadius_(),
    omega_(Zero),
    linearVelocity_(Zero)
{
    vehicleMotion::read(SBMFCoeffs);
    setIncrementalMotion(vehicleMotion::typeName);
    updateMotionCoeffs();
    CofR0_ = CofR_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::vehicleMotion::~vehicleMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::vehicleMotion::transformation() const
{
    updateTime();
    updateMotionCoeffs();

    const scalar t = time_.value();
    septernion TR(septernion::zero);

    if (mag(omega_) > VSMALL)
    {
        // Approximation
        const scalar angle = dt_*mag(omega_);
        quaternion R(frame().axis(), angle);
        TR = septernion(-CofR_)*R*septernion(CofR_);
    }
    else
    {
        const vector displacement = vehicleSpeed_->integrate(tOld_, t);
        quaternion R(1);
        TR = septernion(-displacement)*R;
    }

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::vehicleMotion::velocity() const
{
    updateMotionCoeffs();
    return vectorTuple(linearVelocity_, omega_);
}


void Foam::solidBodyMotionFunctions::vehicleMotion::updateMotionCoeffs() const
{
    const scalar t = time_.value();

    scalar r = 0.0;

    if (turnRadius_.valid())
    {
        r = turnRadius_->value(t);
    }
    else if (invTurnRadius_.valid())
    {
        r = invTurnRadius_->value(t);

        if (r != 0.0)
        {
            r = 1/r;
        }
    }

    CofR_ = frame().coorSys().origin();
    if (r == 0.0)
    {
        omega_ = Zero;
        linearVelocity_ = frame().coorSys().e1()*vehicleSpeed_->value(t);
    }
    else
    {
        CofR_ += frame().coorSys().e2()*r;
        omega_ = vehicleSpeed_->value(t)/r*frame().axis();
        linearVelocity_ = Zero;
    }
}


bool Foam::solidBodyMotionFunctions::vehicleMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    vehicleSpeed_ = Function1<scalar>::New("vehicleSpeed", SBMFCoeffs_);

    if (SBMFCoeffs_.found("turnRadius"))
    {
        turnRadius_ = Function1<scalar>::New("turnRadius", SBMFCoeffs_);
    }

    if (SBMFCoeffs_.found("invTurnRadius"))
    {
        invTurnRadius_ = Function1<scalar>::New("invTurnRadius", SBMFCoeffs_);
    }

    vehicleMotion::updateMotionCoeffs();

    return true;
}


// ************************************************************************* //

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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyMotionFunction, 0);
    defineRunTimeSelectionTable(solidBodyMotionFunction, dictionary);
    defineRunTimeSelectionTable(solidBodyMotionFunction, registry);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunction::solidBodyMotionFunction
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    time_(runTime),
    SBMFCoeffs_
    (
        SBMFCoeffs.optionalSubDict
        (
            word(SBMFCoeffs.lookup("solidBodyMotionFunction")) + "Coeffs"
        )
    ),
    dt_(runTime.deltaT().value()),
    dtOld_(runTime.deltaT().value()),
    // t0_(SBMFCoeffs_.lookupOrDefault<scalar>("t0", 0.0)),
    tOld_
    (
        runTime.value() > 0 ? runTime.value() - runTime.deltaT().value() : 0
    ),
    updateIndex_(0),
    incrementalMotion_
    (
        SBMFCoeffs_.lookupOrDefault<Switch>("incrementalMotion", false)
    ),
    frameName_(frameName),
    coorFramePtr_(nullptr)
{
    if (SBMFCoeffs_.found("coordinateFrame"))
    {
        frameName_ = word(SBMFCoeffs_.lookup("coordinateFrame"));
    }

    // Backward compatibility
    if (!SBMFCoeffs_.found("incrementalMotion"))
    {
        incrementalMotion_ =
            SBMFCoeffs_.lookupOrDefault<Switch>("relativeMotion", false);
    }

    if (SBMFCoeffs_.found("t0"))
    {
        WarningInFunction
            << "The start time isn't supported yet." << endl;
    }
}


Foam::solidBodyMotionFunction::solidBodyMotionFunction
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    time_(obr.time()),
    SBMFCoeffs_
    (
        SBMFCoeffs.optionalSubDict
        (
            word
            (
                SBMFCoeffs.found("type")
              ? SBMFCoeffs.lookup("type")
              : SBMFCoeffs.lookup(solidBodyMotionFunction::typeName)
            ) + "Coeffs"
        )
    ),
    dt_(obr.time().deltaT().value()),
    dtOld_(obr.time().deltaT().value()),
    // t0_(SBMFCoeffs_.lookupOrDefault<scalar>("t0", 0.0)),
    tOld_
    (
        obr.time().value() > 0
      ? obr.time().value() - obr.time().deltaT().value()
      : 0
    ),
    updateIndex_(0),
    incrementalMotion_
    (
        SBMFCoeffs_.lookupOrDefault<Switch>("incrementalMotion", false)
    ),
    frameName_(frameName),
    coorFramePtr_(&obr.lookupObjectRef<coordinateFrame>(frameName_))
{
    if (SBMFCoeffs_.found("t0"))
    {
        WarningInFunction
            << "The start time isn't supported yet." << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunction::~solidBodyMotionFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion Foam::solidBodyMotionFunction::transformation() const
{
    updateTime();

    // Incremental motion
    if (isIncrementalMotion())
    {
        return transformation(tOld_, time_.value());
    }

    // Time 0 motion
    return transformation(0.0, time_.value());
}


Foam::vectorTuple Foam::solidBodyMotionFunction::acceleration() const
{
    updateTime();
    vector linearAccel = Zero;
    vector angularAccel = Zero;

    if (hasFrame())
    {
        const scalar a0 = (dtOld_ + 2.0*dt_)/dt_/(dtOld_ + dt_);
        const scalar a1 = -(dtOld_ + dt_)/dtOld_/dt_;
        const scalar a2 = dt_/dtOld_/(dtOld_ + dt_);
        const vectorTuple v(frame().velocity());
        const vectorTuple vOld(frame().oldTime().velocity());
        const vectorTuple vOldOld(frame().oldTime().oldTime().velocity());

        // Estimate Linear Acceleration
        linearAccel = a0*v.first() + a1*vOld.first() + a2*vOldOld.first();

        // Estimate Angular Acceleration
        angularAccel = a0*v.second() + a1*vOld.second() + a2*vOldOld.second();
    }

    return vectorTuple(linearAccel, angularAccel);
}


const Foam::vector& Foam::solidBodyMotionFunction::CofR() const
{
    return frame().coorSys().origin();
}


const Foam::vector& Foam::solidBodyMotionFunction::CofR0() const
{
    return frame().coorSys0().origin();
}


void Foam::solidBodyMotionFunction::updateTime() const
{
    if (time_.value() > 0 && updateIndex_ != time_.timeIndex())
    {
        tOld_ = time_.time().value() - time_.time().deltaTValue();
        dtOld_ = dt_;
        dt_ = time_.deltaT().value();
    }
}


void Foam::solidBodyMotionFunction::setIncrementalMotion(const word& type)
{
    if
    (
        SBMFCoeffs_.found("incrementalMotion")
     && SBMFCoeffs_.lookup<Switch>("incrementalMotion") == false
    )
    {
        FatalErrorInFunction
            << "The motion function " << type
            << " is only implemented for incremental motion."
            << nl << exit(FatalError);
    }
    incrementalMotion_ = true;
}


Foam::Switch& Foam::solidBodyMotionFunction::isIncrementalMotion() const
{
    return const_cast<Switch&>(incrementalMotion_);
}


bool Foam::solidBodyMotionFunction::read(const dictionary& SBMFCoeffs)
{
    SBMFCoeffs_ =
        SBMFCoeffs.optionalSubDict
        (
            word
            (
                SBMFCoeffs.found("type")
              ? SBMFCoeffs.lookup("type")
              : SBMFCoeffs.lookup(solidBodyMotionFunction::typeName)
            ) + "Coeffs"
        );
    // t0_ = SBMFCoeffs_.lookupOrDefault<scalar>("t0", 0.0);

    incrementalMotion_ =
        SBMFCoeffs_.lookupOrDefault<Switch>("incrementalMotion", false);
    return true;
}


void Foam::solidBodyMotionFunction::writeData(Ostream& os) const
{
    os << SBMFCoeffs_;
}


// ************************************************************************* //

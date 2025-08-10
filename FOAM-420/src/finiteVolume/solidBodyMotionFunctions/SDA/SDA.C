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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solidBodyMotionFunctions/SDA/SDA.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(SDA, 0);
    addToRunTimeSelectionTable(solidBodyMotionFunction, SDA, dictionary);
    addToRunTimeSelectionTable(solidBodyMotionFunction, SDA, registry);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::SDA::SDA
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName),
    CofG_(SBMFCoeffs_.lookup("CofG"))
{
    SDA::read(SBMFCoeffs);
}

Foam::solidBodyMotionFunctions::SDA::SDA
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName),
    CofG_(SBMFCoeffs_.lookup("CofG"))
{
    if (isIncrementalMotion())
    {
        FatalErrorInFunction
            << "The SDA function does not support incremental motion"
            << exit(FatalError);
    }
    SDA::read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::SDA::~SDA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::SDA::transformation
(
    const scalar t0,
    const scalar t
) const
{
    scalar tRel = t - t0;

    scalar Tpi = Tp_ + dTp_*(tRel/dTi_);   // Current roll period [sec]
    scalar wr = twoPi/Tpi; // Current Freq [/sec]

    // Current Phase for roll [rad]
    scalar r = dTp_/dTi_;
    scalar u = Tp_ + r*tRel;
    scalar phr = twoPi*((Tp_/u - 1) + log(mag(u)) - log(Tp_))/r;

    // Current Phase for Sway [rad]
    scalar phs = phr + pi;

    // Current Phase for Heave [rad]
    scalar phh = phr + piByTwo;

    scalar rollA = max(rollAmax_*exp(-sqr(Tpi - Tpn_)/(2*Q_)), rollAmin_);

    vector T
    (
        0,
        swayA_*(sin(wr*tRel + phs) - sin(phs)),
        heaveA_*(sin(wr*tRel + phh) - sin(phh))
    );
    quaternion R(quaternion::XYZ, vector(rollA*sin(wr*tRel + phr), 0, 0));
    septernion TR(septernion(-CofG_ - T)*R*septernion(CofG_));

    DebugInFunction << "Time = " << t << ", t0 = " << t0
                    << ", transformation: " << TR << endl;

    return TR;
}


Foam::vectorTuple
Foam::solidBodyMotionFunctions::SDA::velocity() const
{
    scalar tRel = time_.deltaTValue();

    scalar Tpi = Tp_ + dTp_*(tRel/dTi_);   // Current roll period [sec]
    scalar wr = twoPi/Tpi; // Current Freq [/sec]

    // Current Phase for roll [rad]
    scalar r = dTp_/dTi_;
    scalar u = Tp_ + r*tRel;
    scalar phr = twoPi*((Tp_/u - 1) + log(mag(u)) - log(Tp_))/r;

    // Current Phase for Sway [rad]
    scalar phs = phr + pi;

    // Current Phase for Heave [rad]
    scalar phh = phr + piByTwo;

    scalar rollA = max(rollAmax_*exp(-sqr(Tpi - Tpn_)/(2*Q_)), rollAmin_);

    vector T
    (
        0,
        swayA_*(sin(wr*tRel + phs) - sin(phs)),
        heaveA_*(sin(wr*tRel + phh) - sin(phh))
    );

    return vectorTuple(T/tRel, vector(rollA*sin(wr*tRel + phr), 0, 0)/tRel);
}


bool Foam::solidBodyMotionFunctions::SDA::read(const dictionary& SBMFCoeffs)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("CofG") >> CofG_;
    SBMFCoeffs_.lookup("lamda") >> lamda_;
    SBMFCoeffs_.lookup("rollAmax") >> rollAmax_;
    SBMFCoeffs_.lookup("rollAmin") >> rollAmin_;
    SBMFCoeffs_.lookup("heaveA") >> heaveA_;
    SBMFCoeffs_.lookup("swayA") >> swayA_;
    SBMFCoeffs_.lookup("Q") >> Q_;
    SBMFCoeffs_.lookup("Tp") >> Tp_;
    SBMFCoeffs_.lookup("Tpn") >> Tpn_;
    SBMFCoeffs_.lookup("dTi") >> dTi_;
    SBMFCoeffs_.lookup("dTp") >> dTp_;

    // Rescale parameters according to the given scale parameter
    if (lamda_ > 1 + SMALL)
    {
        heaveA_ /= lamda_;
        swayA_ /= lamda_;
        Tp_ /= sqrt(lamda_);
        Tpn_ /= sqrt(lamda_);
        dTi_ /= sqrt(lamda_);
        dTp_ /= sqrt(lamda_);
    }

    return true;
}


// ************************************************************************* //

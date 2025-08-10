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

#include "controllers.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::controllers::coeffsMode,
        2
    >::names[] =
    {
        "fixed",
        "calculated"
    };
    template<>
    const char* Foam::NamedEnum
    <
        Foam::controllers::tuningMode,
        8
    >::names[] =
    {
        "P",
        "PI",
        "PID",
        "classicPID",
        "pessenIntegralRule",
        "someOvershoot",
        "noOvershoot",
        "nuutinen"
    };

    const NamedEnum
    <
    Foam::controllers::coeffsMode,
    2
    >  Foam::controllers::coeffsModeNames_;

    const NamedEnum
    <
    Foam::controllers::tuningMode,
    8
    >  Foam::controllers::tuningModeNames_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::controllers::controllers
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    coeffsMode_(calculated),
    tuningMode_(PID),
    mesh_(mesh),
    P_(1.0),
    I_(0.0),
    D_(0.0),
    error_(0.0),
    errorIntegral_(0.0),
    oldError_(0.0),
    oldErrorIntegral_(0.0),
    controlParameter_(0.0),
    startTime_(0.0),
    timeIndex_(mesh_.time().timeIndex()),
    setParams_(false),
    printControl_(false)
{
    if (dict.found("controller"))
    {
        controllerDict_ = dict.subDict("controller");

        coeffsMode_ =
            coeffsModeNames_.read(controllerDict_.lookup("coeffsType"));

        switch (coeffsMode_)
        {
            case fixed:
                controllerDict_.lookup("P") >> P_;
                controllerDict_.lookup("I") >> I_;
                D_ = controllerDict_.lookupOrDefault<scalar>("D", 0.0);
                setParams_ = true;
                break;
            case calculated:
                tuningMode_ =
                    tuningModeNames_.read(controllerDict_.lookup("tuningMode"));
                break;
        }

        controllerDict_.lookup("startTime") >> startTime_;
        printControl_ = controllerDict_.lookupOrDefault<bool>("printControl", false);

        Info<< "Controller " << controllerDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::controllers::~controllers()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::controllers::zieglerNichols
(
    const scalar& ultimateGain,
    const scalar& ultimatePeriod
)
{
    switch (tuningMode_)
    {
        case P:
            P_ = 0.5*ultimateGain;
            break;

        case PI:
            P_ = 0.45*ultimateGain;
            I_ = 0.54*ultimateGain/ultimatePeriod;
            break;

        case PID:
            P_ = 0.8*ultimateGain;
            D_ = 0.1*ultimateGain*ultimatePeriod;
            break;

        case classicPID:
            P_ = 0.6*ultimateGain;
            I_ = 1.2*ultimateGain/ultimatePeriod;
            D_ = 0.075*ultimateGain*ultimatePeriod;
            break;

        case pessenIntegralRule:
            P_ = 0.7*ultimateGain;
            I_ = 1.75*ultimateGain/ultimatePeriod;
            D_ = 0.105*ultimateGain*ultimatePeriod;
            break;

        case someOvershoot:
            P_ = 0.33*ultimateGain;
            I_ = 0.66*ultimateGain/ultimatePeriod;
            D_ = 0.11*ultimateGain*ultimatePeriod;
            break;

        case noOvershoot:
            P_ = 0.2*ultimateGain;
            I_ = 0.4*ultimateGain/ultimatePeriod;
            D_ = 0.066*ultimateGain*ultimatePeriod;
            break;

        case nuutinen:
            const scalar T = readScalar(controllerDict_.lookup("estimatedThrust"));
            const scalar n = readScalar(controllerDict_.lookup("estimatedRPS"));

            P_ = 0.5*pow(n, 2)/T;
            I_ = 0.1*P_;
            break;
    }
}

void Foam::controllers::calcParameters
(
    const scalar& inSetPoint,
    const scalar& finSetPoint
)
{
    const scalar t = mesh_.time().value();

    // based on the response curve after the step change (approximation)
    // We do not actually apply the step change (approx method)
    if (inSetPoint >= 0.632*finSetPoint && !setParams_)
    {
        const scalar Ku = 1.27323*P_; // = 4*P_/pi
        const scalar Tu = 0.833333*t; // = t/1.2

        zieglerNichols(Ku, Tu);
        setParams_ = true;
    }
}

Foam::scalar Foam::controllers::executePID
(
    scalar& inSetPoint,
    scalar& finSetPoint
)
{
    const scalar t = mesh_.time().value();
    const scalar dt = mesh_.time().deltaTValue();

    if (t >= startTime_)
    {
        calcParameters(inSetPoint, finSetPoint);

        // Update the old-time quantities
        if (timeIndex_ != mesh_.time().timeIndex())
        {
            timeIndex_ = mesh_.time().timeIndex();
            oldError_ = error_;
            oldErrorIntegral_ = errorIntegral_;
        }

        error_ = inSetPoint - finSetPoint;
        errorIntegral_ += error_*dt;
    }

    const scalar errorDifferential = (error_ - oldError_)/dt;

    if (printControl_)
    {
        Info<< "inSetPoint " << inSetPoint << endl;
        Info<< "finSetPoint " << finSetPoint << endl;
        Info<< "error " << error_ << endl;
        Info<< "P = " << P_ << " I = " << I_ << " D = " << D_ << endl;
        Info<< "Prop = " << P_*error_ << " Int = " << I_*errorIntegral_ << " Der = " << D_*errorDifferential << endl;
    }

    return (P_*error_ + I_*errorIntegral_ + D_*errorDifferential)*dt;
}


// ************************************************************************* //

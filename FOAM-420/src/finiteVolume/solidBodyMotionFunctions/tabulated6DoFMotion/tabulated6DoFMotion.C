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

#include "solidBodyMotionFunctions/tabulated6DoFMotion/tabulated6DoFMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "primitives/Tuple2/Tuple2.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "interpolations/interpolateSplineXY/interpolateSplineXY.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(tabulated6DoFMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        tabulated6DoFMotion,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        tabulated6DoFMotion,
        registry
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tabulated6DoFMotion::tabulated6DoFMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName)
{
    tabulated6DoFMotion::read(SBMFCoeffs);
}

Foam::solidBodyMotionFunctions::tabulated6DoFMotion::tabulated6DoFMotion
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
:
    solidBodyMotionFunction(obr, SBMFCoeffs, frameName)
{
    tabulated6DoFMotion::read(SBMFCoeffs);

    if (isIncrementalMotion())
    {
        FatalErrorInFunction
            << "The tabulated6DoFMotion function "
            << "does not support incremental motion"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tabulated6DoFMotion::~tabulated6DoFMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::tabulated6DoFMotion::transformation
(
    const scalar t0,
    const scalar t
) const
{
    scalar tRel = t - t0;

    if (tRel < times_[0])
    {
        FatalErrorInFunction
            << "current time (" << tRel
            << ") is less than the minimum in the data table ("
            << times_[0] << ')'
            << exit(FatalError);
    }

    if (tRel > times_.last())
    {
        FatalErrorInFunction
            << "current time (" << tRel
            << ") is greater than the maximum in the data table ("
            << times_.last() << ')'
            << exit(FatalError);
    }

    translationRotationVectors TRV = interpolateSplineXY
    (
        tRel,
        times_,
        values_
    );

    // Convert the rotational motion from deg to rad
    TRV[1] *= pi/180.0;

    quaternion R(quaternion::XYZ, TRV[1]);
    septernion TR(septernion(-CofG_ + -TRV[0])*R*septernion(CofG_));

    DebugInFunction << "Time = " << t << ", t0 = " << t0
                    << ", transformation: " << TR << endl;

    return TR;
}

Foam::vectorTuple
Foam::solidBodyMotionFunctions::tabulated6DoFMotion::velocity() const
{
    scalar tRel = time_.deltaTValue();

    if (tRel < times_[0])
    {
        FatalErrorInFunction
            << "current time (" << tRel
            << ") is less than the minimum in the data table ("
            << times_[0] << ')'
            << exit(FatalError);
    }

    if (tRel > times_.last())
    {
        FatalErrorInFunction
            << "current time (" << tRel
            << ") is greater than the maximum in the data table ("
            << times_.last() << ')'
            << exit(FatalError);
    }

    translationRotationVectors TRV =
        interpolateSplineXY(tRel, times_, values_);

    // Convert the rotational motion from deg to rad
    TRV[1] *= pi/180.0;

    return vectorTuple(TRV[0]/tRel, TRV[1]/tRel);
}


bool Foam::solidBodyMotionFunctions::tabulated6DoFMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    // If the timeDataFileName has changed read the file

    fileName newTimeDataFileName
    (
        fileName(SBMFCoeffs_.lookup("timeDataFileName")).expand()
    );

    if (newTimeDataFileName != timeDataFileName_)
    {
        timeDataFileName_ = newTimeDataFileName;

        IFstream dataStream(timeDataFileName_);

        if (dataStream.good())
        {
            List<Tuple2<scalar, translationRotationVectors>> timeValues
            (
                dataStream
            );

            times_.setSize(timeValues.size());
            values_.setSize(timeValues.size());

            forAll(timeValues, i)
            {
                times_[i] = timeValues[i].first();
                values_[i] = timeValues[i].second();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot open time data file " << timeDataFileName_
                << exit(FatalError);
        }
    }

    SBMFCoeffs_.lookup("CofG") >> CofG_;

    return true;
}


// ************************************************************************* //

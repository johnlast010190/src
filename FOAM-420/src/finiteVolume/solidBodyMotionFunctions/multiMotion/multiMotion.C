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

#include "solidBodyMotionFunctions/multiMotion/multiMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(multiMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        multiMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::multiMotion::multiMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime, frameName)
{
    multiMotion::read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::multiMotion::~multiMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::multiMotion::transformation
(
    const scalar t0,
    const scalar t
) const
{
    const septernion TR = multiMotionTransform(0.0, t);
    if (isIncrementalMotion())
    {
        // check all motions
        for (label i = 0; i < SBMFs_.size(); i++)
        {
            if (SBMFs_[i].isIncrementalMotion())
            {
                FatalErrorInFunction
                    << "Inconsistent multiMotion setup. Associated motion"
                    << nl << "functions must all be absolute."
                    << exit(FatalError);
            }
        }

        return TR/multiMotionTransform(0.0, t - time_.deltaT().value());
    }

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::multiMotion::multiMotionTransform
(
    const scalar t0,
    const scalar t
) const
{
    septernion TR = SBMFs_[0].transformation(t0, t);

    for (label i = 1; i < SBMFs_.size(); i++)
    {
        TR *= SBMFs_[i].transformation(t0, t);
    }

    DebugInFunction << "Time = " << t << ", t0 = " << t0
                    << ", transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::multiMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    label i = 0;
    SBMFs_.setSize(SBMFCoeffs_.size());

    forAllConstIter(IDLList<entry>, SBMFCoeffs_, iter)
    {
        if (iter().isDict())
        {
            SBMFs_.set
            (
                i,
                solidBodyMotionFunction::New(iter().dict(), time_)
            );

            Info<< "Constructed SBMF " << i << " : "
                << iter().keyword() << " of type "
                << SBMFs_[i].type() << endl;

            i++;
        }
    }
    SBMFs_.setSize(i);

    return true;
}


// ************************************************************************* //

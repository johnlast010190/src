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

Class
    Foam::controllers

Description
    Controllers
    It implements PID controllers

SourceFiles
    controllers.C

\*---------------------------------------------------------------------------*/

#ifndef controllers_H
#define controllers_H

#include "db/dictionary/dictionary.H"
#include "forces/forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class controllers Declaration
\*---------------------------------------------------------------------------*/

class controllers
{
public:

    // Public data

        //- Coeffs mode for PID, if input-based or calculated (tuningMode)
        enum coeffsMode
        {
            fixed,
            calculated
        };

        //- Tuning mode to calculate PID parameters
        enum tuningMode
        {
            P,
            PI,
            PID,
            classicPID,
            pessenIntegralRule,
            someOvershoot,
            noOvershoot,
            nuutinen
        };

        //- Coeffs mode source type names
        static const NamedEnum<coeffsMode, 2> coeffsModeNames_;

        //- Tuning mode source type names
        static const NamedEnum<tuningMode, 8> tuningModeNames_;

private:

    // Private Data

        //- Controller mode type
        coeffsMode coeffsMode_;

        //- Controller mode type
        tuningMode tuningMode_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Controllers dictionary
        dictionary controllerDict_;

        //- Proportional gain
        scalar P_;

        //- Integral gain
        scalar I_;

        //- Derivative gain
        scalar D_;

        //- Error
        scalar error_;

        //- Error integral
        scalar errorIntegral_;

        //- Old error
        scalar oldError_;

        //- Old error integral
        scalar oldErrorIntegral_;

        //- Control parameter
        scalar controlParameter_;

        //- Old control parameter
        scalar oldControlParameter_;

        //- Ramp time for the controller
        scalar startTime_;

        //- Time index of the last update
        label timeIndex_;

        //- Set PID parameters
        bool setParams_;

        //- Print controller info
        bool printControl_;

    // Private Member Functions

        //- Disallow copy construct
        controllers(const controllers&);

        //- Disallow default bitwise assignment
        void operator=(const controllers&);

public:

    // Constructors

        //- Construct from components
        controllers
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~controllers();


    // Member Functions

        //- Ziegler-Nichols method for calculating PID parameters (approx method)
        void zieglerNichols(const scalar& ultimateGain, const scalar& ultimatePeriod);

        //- Method for setting PID parameters using Ziegler-Nichols method
        void calcParameters(const scalar& inSetPoint, const scalar& finSetPoint);

        //- PID controller applied to a target parameter: finSetPoint
        scalar executePID(scalar& inSetPoint, scalar& finSetPoint);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

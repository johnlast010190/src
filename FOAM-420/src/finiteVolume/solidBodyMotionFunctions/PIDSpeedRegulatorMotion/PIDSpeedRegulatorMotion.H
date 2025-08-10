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
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::solidBodyMotionFunctions::PIDSpeedRegulatorMotion

Description
    Propeller speed adjusted using PID controller

SourceFiles
    PIDSpeedRegulatorMotion.C

\*---------------------------------------------------------------------------*/

#ifndef PIDSpeedRegulatorMotion_H
#define PIDSpeedRegulatorMotion_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class PIDSpeedRegulatorMotion Declaration
\*---------------------------------------------------------------------------*/

class PIDSpeedRegulatorMotion
:
    public solidBodyMotionFunction
{

        class regulatorItem
        {
        public:

            scalar steeringAngle_;

            scalar propellerSpeed_;

            scalar propellerTorque_;

            scalar torqueLimitAtSpeed_;

            scalar regulatedSpeed_;

            regulatorItem()
            :
                steeringAngle_(0),
                propellerSpeed_(0),
                propellerTorque_(0),
                torqueLimitAtSpeed_(0),
                regulatedSpeed_(0)
            {};

            regulatorItem
            (
                const scalar steeringAngle,
                const scalar propellerSpeed,
                const scalar propellerTorque,
                const scalar torqueLimitAtSpeed,
                const scalar regulatedSpeed
            )
            :
                steeringAngle_(steeringAngle),
                propellerSpeed_(propellerSpeed),
                propellerTorque_(propellerTorque),
                torqueLimitAtSpeed_(torqueLimitAtSpeed),
                regulatedSpeed_(regulatedSpeed)
            {}

            scalar steeringAngle() const {return steeringAngle_;}
            scalar propellerSpeed() const {return propellerSpeed_;}
            scalar propellerTorque() const {return propellerTorque_;}
            scalar torqueLimitAtSpeed()const {return torqueLimitAtSpeed_;}
            scalar regulatedSpeed() const {return regulatedSpeed_;}
        };

        //- angle in radians
        mutable scalar angle_;

        //- Propeller speed (rad/s)
        mutable scalar omega_;

        //- Origin of the axis
        mutable vector origin_;

        //- Axis vector
        mutable vector axis_;

        //- Forces function object name
        word functionName_;

        //- Settings for lower and upper speed limits
        Pair<scalar> speedLimits_;

        //- Regulator parameter
        scalar kP_;

        //- Regulator parameter
        scalar tN_;

        //- Regulator parameter
        scalar tV_;

        //- Steering angle (omega)
        scalar steeringSpeed_;

        //- Design point for omega and torque
        Pair<scalar> designPoint_;

        //- Torque limit curve
        autoPtr<Function1<scalar>> limitCurve_;

        //- Time to start regulating
        scalar startTime_;

        mutable DynamicList<regulatorItem> results_;

    // Private Member Functions

        //- Disallow copy construct
        PIDSpeedRegulatorMotion(const PIDSpeedRegulatorMotion&);

        //- Disallow default bitwise assignment
        void operator=(const PIDSpeedRegulatorMotion&);

        //- Return lookup of propeller torque
        scalar torque() const;

public:

    //- Runtime type information
    TypeName("PIDSpeedRegulatorMotion");


    // Constructors

        //- Construct from components
        PIDSpeedRegulatorMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        //- Construct from components
        PIDSpeedRegulatorMotion
        (
            const objectRegistry& obr,
            const dictionary& SBMFCoeffs,
            const word& frameName
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new PIDSpeedRegulatorMotion
                (
                    SBMFCoeffs_,
                    time_,
                    frameName_
                )
            );
        }


    //- Destructor
    virtual ~PIDSpeedRegulatorMotion();


    // Member Functions

        //- Function that updates origin and axis for incemental motion
        void updateOriginAndAxis() const;

        //- Return the solid-body motion transformation septernion from t0 to t
        virtual septernion transformation
        (
            const scalar t0,
            const scalar t
        ) const override;

        //- Return the solid-body motion transformation septernion from t0 to t
        virtual vectorTuple velocity() const override;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

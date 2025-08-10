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
    (c) 2017 Esi Ltd.

Class
    Foam::solidBodyMotionFunctions::rotatingStepMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function.

    The rotation is defined by an origin and axis of rotation and an angular
    speed.

SourceFiles
    rotatingStepMotion.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingStepMotion_H
#define rotatingStepMotion_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class rotatingStepMotion Declaration
\*---------------------------------------------------------------------------*/

class rotatingStepMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- parent object registry
        const objectRegistry& obr_;

        //- Origin of the axis
        const vector origin_;

        //- Axis vector
        const vector axis_;

        //- Step angle (degrees)
        scalar theta_;

        //- Initial angle at index 0
        scalar theta0_;

        //- step period (iterations)
        label stepPeriod_;

        //- store rotation septernion so update can be detected
        mutable septernion stO_;

        //- name of velocity field to rotate
        word UName_;

        //- index of last update
        //- to prevent multiple update calls performing repeat transformations
        mutable label lastUpdate_;


    // Private Member Functions

        //- Disallow copy construct
        rotatingStepMotion(const rotatingStepMotion&);

        //- Disallow default bitwise assignment
        void operator=(const rotatingStepMotion&);


public:

    //- Runtime type information
    TypeName("rotatingStepMotion");


    // Constructors

        //- Construct from components
        rotatingStepMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        rotatingStepMotion
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
                new rotatingStepMotion
                (
                    SBMFCoeffs_,
                    time_,
                    frameName_
                )
            );
        }

    //- Destructor
    virtual ~rotatingStepMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion from t0 to t
        virtual septernion transformation
        (
            const scalar t0,
            const scalar t
        ) const override;

        //- Return the motion velocity
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

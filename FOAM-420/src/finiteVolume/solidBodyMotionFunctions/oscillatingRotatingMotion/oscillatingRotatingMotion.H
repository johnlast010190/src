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

Class
    Foam::solidBodyMotionFunctions::oscillatingRotatingMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating rotation.

SourceFiles
    oscillatingRotatingMotion.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingRotatingMotion_H
#define oscillatingRotatingMotion_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class oscillatingRotatingMotion Declaration
\*---------------------------------------------------------------------------*/

class oscillatingRotatingMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        mutable vector origin_;

        //- Amplitude
        autoPtr<Function1<vector>> amplitude_;

        //- Radial velocity
        autoPtr<Function1<scalar>> omega_;


    // Private Member Functions

        //- Disallow copy construct
        oscillatingRotatingMotion(const oscillatingRotatingMotion&);

        //- Disallow default bitwise assignment
        void operator=(const oscillatingRotatingMotion&);


public:

    //- Runtime type information
    TypeName("oscillatingRotatingMotion");


    // Constructors

        //- Construct from components
        oscillatingRotatingMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        oscillatingRotatingMotion
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
                new oscillatingRotatingMotion
                (
                    SBMFCoeffs_,
                    time_,
                    frameName_
                )
            );
        }


    //- Destructor
    virtual ~oscillatingRotatingMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation
        (
            const scalar t0,
            const scalar t
        ) const override;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;

        //- Return the motion velocity
        virtual vectorTuple velocity() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

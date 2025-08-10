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
    Foam::solidBodyMotionFunctions::multiMotion

Description
    Combination of SolidBodyMotionFvMesh 6DoF motion functions.

SourceFiles
    multiMotion.C

\*---------------------------------------------------------------------------*/

#ifndef multiMotion_H
#define multiMotion_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multiMotion Declaration
\*---------------------------------------------------------------------------*/

class multiMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Motions to combine
        PtrList<solidBodyMotionFunction> SBMFs_;


    // Private Member Functions

        //- Disallow copy construct
        multiMotion(const multiMotion&);

        //- Disallow default bitwise assignment
        void operator=(const multiMotion&);


public:

    //- Runtime type information
    TypeName("multiMotion");


    // Constructors

        //- Construct from components
        multiMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime,
            const word& frameName = "none"
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new multiMotion
                (
                    SBMFCoeffs_,
                    time_,
                    frameName_
                )
            );
        }


    //- Destructor
    virtual ~multiMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion from t0 to t
        virtual septernion transformation
        (
            const scalar t0,
            const scalar t
        ) const override;

        //- Return the solid-body motion transformation septernion from t0 to t
        virtual septernion multiMotionTransform
        (
            const scalar t0,
            const scalar t
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;

        //- Multimotion shouldn't be used to calc velocity
        virtual vectorTuple velocity() const override
        {
            FatalErrorInFunction
                << "No times selected." << exit(FatalError);

            return vectorTuple(Zero, Zero);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

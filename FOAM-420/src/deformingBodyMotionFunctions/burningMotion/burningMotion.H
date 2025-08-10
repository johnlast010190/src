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
    (c) 2018-2022 Esi Ltd.

Class
    Foam::deformingBodyMotionFunctions::burningMotion

Description
    Motion for burning interface.
    This class calculates the interface velocity and velocity boundary values

SourceFiles
    burningMotion.C

\*---------------------------------------------------------------------------*/

#ifndef burningMotion_H
#define burningMotion_H

#include "dynamicGIBFvMesh/deformingBodyMotionFunction/deformingBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"
#include "fields/Fields/vectorField/vectorField.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#ifndef FOAM_SWAK4FOAM_DISABLED
#include "expressionSource.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace deformingBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class burningMotion Declaration
\*---------------------------------------------------------------------------*/

class burningMotion
:
    public deformingBodyMotionFunction
{
    // Private data

            //- type of motion
            word mode_;

            //- Solid density
            scalar rhos_;

            //- constants
            autoPtr<Function1<scalar>> uniformBurnRate_;

                        //- User defined field name
            word fieldName_;

            #ifndef FOAM_SWAK4FOAM_DISABLED
            //- expression entered by user
            mutable autoPtr<expressionSource<scalar>> expression_;
            #endif

            //- Scale factor to accelerate the burn progression without
            //  influencing boundary velocity
            scalar scaleFactor_;


    // Private Member Functions

        //- Return the interface velocity
        tmp<scalarField> burnrate(const label&) const;

        //- Disallow copy construct
        burningMotion(const burningMotion&);

        //- Disallow default bitwise assignment
        void operator=(const burningMotion&);


public:

    //- Runtime type information
    TypeName("burningMotion");


    // Constructors

        //- Construct from components
        burningMotion
        (
            const fvMesh& mesh,
            const dictionary& DBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~burningMotion();


    // Member Functions

        //- Update properties from given dictionary
        virtual bool read(const dictionary&);

        //- update
        virtual void update();

        //- Return the boundary velocity
        virtual tmp<vectorField> boundaryVelocity(const label&) const;

        //- Return the interface velocity
        virtual tmp<vectorField> interfaceVelocity(const label&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace deformingBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

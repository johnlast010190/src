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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::deformingBodyMotionFunctions::adjointMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Constant velocity displacement.

SourceFiles
    adjointMotion.C

\*---------------------------------------------------------------------------*/

#ifndef adjointMotion_H
#define adjointMotion_H

#include "dynamicGIBFvMesh/deformingBodyMotionFunction/deformingBodyMotionFunction.H"
#include "fields/Fields/primitiveFields.H"
#include "meshes/primitiveShapes/point/point.H"
#include "fields/Fields/vectorField/vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace deformingBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class adjointMotion Declaration
\*---------------------------------------------------------------------------*/

class adjointMotion
:
    public deformingBodyMotionFunction
{
    // Private data

    // Private Member Functions

        //- Return the interface velocity
        tmp<scalarField> normalVelocity(const label&) const;

        //- Disallow default bitwise assignment
        void operator=(const adjointMotion&);


public:

    //- Runtime type information
    TypeName("adjointMotion");


    // Constructors

        //- Construct from components
        adjointMotion
        (
            const fvMesh& mesh,
            const dictionary& DBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~adjointMotion();


    // Member Functions

        //- Update properties from given dictionary
        virtual bool read(const dictionary& DBMFCoeffs);

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

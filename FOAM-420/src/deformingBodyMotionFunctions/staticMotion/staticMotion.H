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
    (c) 2018 Esi Ltd.

Class
    Foam::deformingBodyMotionFunctions::staticMotion

Description
    Motion for burning interface.
    This class calculates the interface velocity and velocity boudnary values

SourceFiles
    staticMotion.C

\*---------------------------------------------------------------------------*/

#ifndef staticMotion_H
#define staticMotion_H

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
                          Class staticMotion Declaration
\*---------------------------------------------------------------------------*/

class staticMotion
:
    public deformingBodyMotionFunction
{
    // Private data

    // Private Member Functions

        //- Return the interface velocity
        tmp<scalarField> normalVelocity(const label&) const;

        //- Disallow copy construct
        staticMotion(const staticMotion&);

        //- Disallow default bitwise assignment
        void operator=(const staticMotion&);


public:

    //- Runtime type information
    TypeName("staticMotion");


    // Constructors

        //- Construct from components
        staticMotion
        (
            const fvMesh& mesh,
            const dictionary& DBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~staticMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
//        virtual septernion transformation() const;

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

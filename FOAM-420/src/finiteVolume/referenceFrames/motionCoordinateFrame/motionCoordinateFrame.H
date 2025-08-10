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
    (c) 2022-2023 Esi Ltd.

Class
    Foam::motionCoordinateFrame

Description
    Moving frame of reference based on motion function.

    Derivation of velocity from rotating acceleration in absolute velocity
    formulation:
    https://openfoamwiki.net/index.php/See_the_MRF_development

Usage
    Example of the vehicle frame specification:
    \verbatim
    <referenceFrameName>
    {
        type            motionCoordinateFrame;
        motionFunction
        {
            type vehicleMotion;
            vehicleMotionCoeffs
            {
                vehicleSpeed    constant 50;
                invTurnRadius   table ( ( 0 0.02 ) ( 2000 0.02 ) );
            }
        }
        coordinateSystem
        {
            type            cartesian;
            origin          ( 0 50 -0.05 );
            e1              ( -1 0 0 );
            e2              ( 0 -1 0 );
        }
    }
    \endverbatim
\*---------------------------------------------------------------------------*/

#ifndef motionCoordinateFrame_H
#define motionCoordinateFrame_H

#include "referenceFrames/coordinateFrame.H"
#include "meshes/MeshObject/MeshObject.H"
#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class solidBodyMotionFunction;

/*---------------------------------------------------------------------------*\
                    Class motionCoordinateFrame Declaration
\*---------------------------------------------------------------------------*/
class motionCoordinateFrame
:
    public coordinateFrame
{

protected:

    // Protected Member Data

        //- Access to the motion function
        autoPtr<solidBodyMotionFunction> motionFunction_;

        //- Dimensionality check
        void dimensionCheck
        (
            const dimensionSet& dim1,
            const dimensionSet& dim2
        ) const;

        //- List of pointers to all the parents
        UPtrList<coordinateFrame> parentFrames_;


public:

    //- Runtime type information
    TypeName("motionCoordinateFrame");

    //- Constructors
        motionCoordinateFrame
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& frameName
        );

    //- Destructor
        virtual ~motionCoordinateFrame();


    // Public Member Functions

        //- Access to motion function
        const solidBodyMotionFunction& motion() const
        {
            return motionFunction_();
        }

        //- Center of rotation
        virtual const vector& CofR() const override
        {
            return motionFunction_().CofR();
        }

        //- Initial center of rotation
        virtual const vector& CofR0() const override
        {
            return motionFunction_().CofR0();
        }

        //- Updating state of coordinate frame
        virtual void updateState() const override;

        //- Return frame relative velocity vector from position(s)
        virtual tmp<volVectorField> frameVelocity
        (
            const volVectorField& positions,
            bool addParentFrames = true
        ) const override;

        //- Return velocity surface vector field from position(s)
        virtual tmp<surfaceVectorField> frameVelocity
        (
            const surfaceVectorField& positions,
            bool addParentFrames = true
        ) const override;

        //- Return velocity vector field from position(s)
        virtual tmp<vectorField> frameVelocity
        (
            const vectorField& positions,
            bool addParentFrames = true
        ) const override;

        //- Return frame velocity for specific point
        virtual vector frameVelocity
        (
            const vector& position,
            bool addParentFrames = true
        ) const override;

        //- Data writer
        virtual bool writeData(Ostream& os) const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}//end namespace FOAM

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

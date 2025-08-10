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
    (c) 2017-2022 Esi Ltd.

Class
    Foam::deformingBodyMotionFunctions::deformingTyreTransformation

Description
    Coordinate transformation for deforming tyres
    User inputs:
        -3D point that defines the center of the tyre
    -Vector that defines the axis of the tyre
    -Vector that defines the normal plane to the ground
    -Radius of the tyre (including tread)
    -Height of the tread (radial distance)
    -Width of the tyre
    -Angle of contact with the ground (full angle)

SourceFiles
    deformingTyreTransformation.C

\*---------------------------------------------------------------------------*/

#ifndef deformingTyreTransformation_H
#define deformingTyreTransformation_H

#include "vNurbsDeformation/controlBox/controlBox.H"
#include "dynamicGIBFvMesh/geometryTransformations/geometryTransformation/geometryTransformation.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "triSurface/triSurface.H"
#include "vNurbsDeformation/volumetricNurbs/volumetricNurbs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class deformingTyreTransformation Declaration
\*---------------------------------------------------------------------------*/

class deformingTyreTransformation
:
    public geometryTransformation
{
    //- private data

        // Center of tyre in 3d
        point center_;

        // vector parallel to the axis of the tyre
        vector tNormal_;

        // vector normal to the ground plane
        vector gNormal_;

        // vector that corresponds to the virtual x-axis
        vector vX_;

        // Radius of tyre
        scalar radius_;

        // Height of tread
        scalar hTread_;

        // Width of the tyre
        scalar tWidth_;

        // Angle of contact with the ground
        scalar alphaContact_;

        // stl displacement to touch the ground
        scalar distToGround_;

        // Morpher with zero initialization
        volumetricNurbs *vNurbs;

    //- private member functions

        // Function to initialize the control box based on the original
        // dimensions of the tyre
        void initializeCBox();

        // Function that defines the appropriate control point displacements
        // so that the tyre deforms according to user input
        void createDisplacementFields();

    public:

        //- Runtime type info
        TypeName("deformingTyreTransformation");

        //- Constructor
        deformingTyreTransformation
        (
            const dictionary& dict
        );

        //- Destructor
        ~deformingTyreTransformation()
        {
            delete vNurbs;
            vNurbs = nullptr;
        }

    //- Member functions

        //- Function that reads the dict
        //- Calls initializeCBox() and createDisplacementFields() internally
        virtual bool read(const dictionary&);

        //- Function that transforms a given pointField
        //- Initially the mopher is trained, then the cBox is updated and
        //- consequently a new pointField is returned
        virtual tmp<pointField> transformPoints
        (
            const pointField&
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

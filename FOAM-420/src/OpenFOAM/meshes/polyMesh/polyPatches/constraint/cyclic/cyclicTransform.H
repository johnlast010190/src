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
    (c) 2020 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::cyclicTransform

Description
    Cyclic plane transformation.

SourceFiles
    cyclicTransform.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicTransform_H
#define cyclicTransform_H

#include "primitives/transform/transformer/transformer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cyclicTransform Declaration
\*---------------------------------------------------------------------------*/

class cyclicTransform
{
public:

    // Public enumerations

        //- Transformation types
        enum transformTypes
        {
            UNSPECIFIED,        // Unspecified -> automatic transformation
            NONE,               // No transformation
            ROTATIONAL,         // Rotation around a coordinate axis
            TRANSLATIONAL       // Translation
        };

        //- Transformation type names
        static const NamedEnum<transformTypes, 4> transformTypeNames;


        //- Methods for the specification of the rotation angle
        enum rotationalMethods
        {
            USERSPECIFIED,      // User input
            CENTREOFMASS,       // Based on patches' centres of mass
            FURTHESIPOINT       // Based on furthest points on specified axis
        };

        //- Rotational method names
        static const NamedEnum<rotationalMethods, 3> rotationalMethodNames;


private:

    // Private data

        //- Type of transformation
        transformTypes transformType_;

        //- Rotational method, if a rotational transform is set
        rotationalMethods rotationalMethod_;

        //- Axis of rotation
        vector rotationAxis_;

        //- Centre of rotation
        point rotationCentre_;

        //- Rotation angle
        scalar rotationAngle_;

        //- Separation vector
        vector separation_;

        //- Is the transformation specification complete, or do some parts of
        //  the above data still need to be computed from the patch geometry?
        bool transformComplete_;

        //- Generic transformer object
        transformer transform_;


    // Private Member Functions

        //- Read and return entry for the transformation type
        word readTransform(const dictionary& dict, const bool defaultIsNone);

        //- Re-calculate the transformer object from the other data if possible
        void update();

        //- Set the transform to that supplied. Return true if the
        //  transformations are compatible. Used for copying transformation
        //  data from the neighbour patch.
        bool set
        (
            const cyclicTransform& t,
            const scalar lengthScale,
            const scalar matchTolerance
        );


public:

    //- Runtime type information
    TypeName("cyclicTransform");


    // Constructors

        // Transformation set NONE
        cyclicTransform();

        // Transformation set as UNSPECIFIED or NONE
        cyclicTransform(const bool defaultIsNone);

        // Transformation read from dictionary
        cyclicTransform(const dictionary& dict, const bool defaultIsNone);

        //- Construct from coupled patch data. Copies from the supplied
        //  transform and neighbour transform. Suitable for
        //  geometrically dissimilar patches.
        cyclicTransform
        (
            const word& name,
            const vectorField& areas,
            const cyclicTransform& transform,
            const word& nbrName,
            const cyclicTransform& nbrTransform,
            const scalar matchTolerance
        );

        //- Construct from coupled patch data. Copies from the supplied
        //  transform and neighbour transform, then tries to determine missing
        //  parts of the transformation automatically from the patch geometry.
        //  Suitable for geometrically similar patches only.
        cyclicTransform
        (
            const word& name,
            const pointField& ctrs,
            const vectorField& areas,
            const cyclicTransform& transform,
            const word& nbrName,
            const pointField& nbrCtrs,
            const vectorField& nbrAreas,
            const cyclicTransform& nbrTransform,
            const scalar matchTolerance
        );


    //- Destructor
    virtual ~cyclicTransform();


    // Member Functions

        // Access

            //- Type of transform
            transformTypes transformType() const
            {
                return transformType_;
            }

            //- Rotational method
            rotationalMethods rotationalMethod() const
            {
                return rotationalMethod_;
            }

            //- Rotation angle
            scalar rotationAngle() const
            {
                return rotationAngle_;
            }

            //- Is the transform fully specified?
            bool transformComplete() const
            {
                return transformComplete_;
            }

            //- Return transformation between the coupled patches
            const transformer& transform() const
            {
                if (!transformComplete_)
                {
                    FatalErrorInFunction
                        << "The transformation has not been fully specified or"
                        << " calculated" << exit(FatalError);
                }

                return transform_;
            }


        //- Write the data to a dictionary
        void write(Ostream& os) const;


        // Global Operators

            friend cyclicTransform operator&
            (
                const transformer& t,
                const cyclicTransform& c
            );

            friend cyclicTransform inv(const cyclicTransform& c);
};


cyclicTransform operator&(const transformer& t, const cyclicTransform& c);

cyclicTransform inv(const cyclicTransform& c);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

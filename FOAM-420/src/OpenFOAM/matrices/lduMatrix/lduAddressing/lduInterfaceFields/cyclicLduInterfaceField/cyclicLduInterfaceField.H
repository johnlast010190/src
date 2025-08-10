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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::cyclicLduInterfaceField

Description
    Abstract base class for cyclic coupled interfaces.

SourceFiles
    cyclicLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicLduInterfaceField_H
#define cyclicLduInterfaceField_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "primitives/transform/transformer/transformer.H"
#include "db/typeInfo/typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicLduInterfaceField
{

public:

    //- Runtime type information
    TypeName("cyclicLduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        cyclicLduInterfaceField()
        {}


    //- Destructor
    virtual ~cyclicLduInterfaceField();


    // Member Functions

        // Access

            //- Return the transformation between the coupled patches
            virtual const transformer& transform() const = 0;

            //- Return rank of component for transform
            virtual int rank() const = 0;

            //- Is the transform required
            bool transforms() const
            {
                return transform().transforms() && rank() != 0;
            }


        //- Transform given patch field
        template<class Type>
        void transformCoupleField(Field<Type>& f) const;

        //- Transform given patch component field
        void transformCoupleField
        (
            scalarField& f,
            const direction cmpt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void Foam::cyclicLduInterfaceField::transformCoupleField
(
    Field<Type>& f
) const
{
    if (transforms())
    {
        transform().transform(f, f);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

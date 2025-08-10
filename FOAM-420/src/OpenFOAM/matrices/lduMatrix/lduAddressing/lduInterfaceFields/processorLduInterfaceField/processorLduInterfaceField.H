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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::processorLduInterfaceField

Description
    Abstract base class for processor coupled interfaces.

SourceFiles
    processorLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef processorLduInterfaceField_H
#define processorLduInterfaceField_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "primitives/transform/transformer/transformer.H"
#include "db/typeInfo/typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class processorLduInterfaceField
{

public:

    //- Runtime type information
    TypeName("processorLduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        processorLduInterfaceField()
        {}


    //- Destructor
    virtual ~processorLduInterfaceField();


    // Member Functions

        // Access

            //- Return communicator used for comms
            virtual label comm() const = 0;

            //- Return processor number
            virtual int myProcNo() const = 0;

            //- Return neigbour processor number
            virtual int neighbProcNo() const = 0;

            //- Return transformation between the coupled patches
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
void Foam::processorLduInterfaceField::transformCoupleField
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

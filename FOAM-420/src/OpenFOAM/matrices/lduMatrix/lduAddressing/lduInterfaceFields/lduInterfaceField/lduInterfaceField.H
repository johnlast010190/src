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
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::lduInterfaceField

Description
    An abstract base class for implicitly-coupled interface fields
    e.g. processor and cyclic patch fields.

SourceFiles
    lduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef lduInterfaceField_H
#define lduInterfaceField_H

#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.H"
#include "fields/Fields/primitiveFieldsFwd.H"
#include "db/IOstreams/Pstreams/Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class lduMatrix;

/*---------------------------------------------------------------------------*\
                     Class lduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class lduInterfaceField
{
    // Private data

        //- Reference to the coupled patch this field is defined for
        const lduInterface& interface_;

        //- Update index used so that updateInterfaceMatrix is called only once
        //  during the construction of the matrix
        bool updatedMatrix_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduInterfaceField(const lduInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const lduInterfaceField&);


public:

    //- Runtime type information
    TypeName("lduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        lduInterfaceField(const lduInterface& patch)
        :
            interface_(patch),
            updatedMatrix_(false)
        {}


    //- Destructor
    virtual ~lduInterfaceField();


    // Member Functions

        // Access

            //- Return the interface
            const lduInterface& interface() const
            {
                return interface_;
            }

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return type();
            }


        // Coupled interface matrix update

            //- Whether matrix has been updated
            bool updatedMatrix() const
            {
                return updatedMatrix_;
            }

            //- Whether matrix has been updated
            bool& updatedMatrix()
            {
                return updatedMatrix_;
            }

            //- Is all data available
            virtual bool ready() const
            {
                return true;
            }

            //- Initialise neighbour matrix update. Add
            //  or subtract coupled contributions to matrix
            virtual void initInterfaceMatrixUpdate
            (
                scalarField&,
                const bool add,
                const scalarField&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Update result field based on interface functionality. Add
            //  or subtract coupled contributions to matrix
            virtual void updateInterfaceMatrix
            (
                scalarField&,
                const bool add,
                const scalarField&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const = 0;

            //- Helper: add (or subtract) weighted contributions to internal
            //  field
            template<class Type>
            void addToInternalField
            (
                Field<Type>& result,
                const bool add,
                const scalarField& coeffs,
                const Field<Type>& vals
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matrices/lduMatrix/lduAddressing/lduInterfaceFields/lduInterfaceField/lduInterfaceFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

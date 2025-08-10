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
    (c) 2022 Esi Ltd.

Class
    Foam::referenceFields

SourceFiles
    referenceFieldsI.H
    referenceFields.C

\*---------------------------------------------------------------------------*/

#ifndef referenceFields_H
#define referenceFields_H

#include "memory/autoPtr/autoPtr.H"
#include "db/regIOobject/regIOobject.H"
#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class referenceFields Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class referenceFields
:
    public regIOobject
{
    // Private member data

        //- Type of geometric field
        typedef GeometricField<Type, fvPatchField, volMesh> geoField;

        //- Original field name
        const word fieldName_;

        //- Store field name
        const word refFieldName_;

        //- Field reference value
        dimensioned<Type> fieldRef_;

        //- Is field constant? (default false)
        bool isConst_;

        //- Value of constant field
        dimensioned<Type> constVal_;

        //- Reference to the filed
        const geoField* field_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Class to create consistent boundary access
        class Boundary
        {
        public:
            Boundary(const referenceFields<Type>& base)
            :
                refFiels_(base)
            {};
            const referenceFields& refFiels_;
            ~Boundary(){};
            tmp<Field<Type>> operator[](const label patchi) const
            {
                return refFiels_.patchField(patchi);
            }
        };

        //- Boundary class
        Boundary boundary_;


public:

    // Constructors

        //- Construct form the components
        referenceFields
        (
            const objectRegistry& obr,
            const word& fieldName,
            const dictionary& dict = dictionary::null
        );


    // Member functions

        //- Return primitive field
        tmp<Field<Type>> primitiveField() const;

        //- Return boundary field
        tmp<Field<Type>> patchField(const label patchi) const;

        //- Boundary field reference
        const Boundary& boundaryField() const
        {
            return boundary_;
        }

        //- Return offset from original value
        const dimensioned<Type>& offset() const;

        //- Return offset from original value
        const dimensioned<Type>& constant() const;

        //- Is field constant?
        bool isConst() const;

        //- Make reference constant (even if not specified in mat. dict)
        void makeConst(const Type& value = Zero);

        //- Make reference non-constant (even if not specified in mat. dict)
        void makeNonConst();

        //- Update fields
        bool updateScalarField
        (
            const word& fieldName,
            const volScalarField& volField
        );

        //- Update fields
        bool updateVectorField
        (
            const word& fieldName,
            const volVectorField& volField
        );

        //- Update fields
        bool updateTensorField
        (
            const word& fieldName,
            const volTensorField& volField
        );


    // IO functions

        //- Write for data
        virtual bool writeData(Ostream& os) const override;

        //- Write
        virtual bool write(const bool valid) const override;


    // Operators

        // Single cell access operator
        inline Type operator[](const label celli) const;

        //- Return dimless geometric field
        tmp<geoField> operator()() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "referenceFieldsI.H"

#ifdef NoRepository
    #include "referenceFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

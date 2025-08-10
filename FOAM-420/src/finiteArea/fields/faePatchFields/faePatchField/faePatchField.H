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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Class
    faePatchField

Description
    faePatchField<Type> abstract base class.  This class gives a fat-interface
    to all derived classes covering all possible ways in which they might be
    used.  The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.  The next
    level of derivation covers all the specialised typed with specific
    evaluation proceedures, particularly with respect to specific fields.

SourceFiles
    faePatchField.C
    newPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef faePatchField_H
#define faePatchField_H

#include "faMesh/faPatches/faPatch/faPatch.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "faEdgeMesh/faEdgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

class objectRegistry;
class dictionary;
class faPatchFieldMapper;
class faEdgeMesh;


// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class faePatchField;

template<class Type>
class calculatedFaePatchField;

template<class Type>
Ostream& operator<<(Ostream&, const faePatchField<Type>&);


/*---------------------------------------------------------------------------*\
                        Class faePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faePatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to a patch
        const faPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, faEdgeMesh>& internalField_;

        //- Optional patch type, used to allow specified boundary conditions
        //  to be applied to constraint patches by providing the constraint
        //  patch type as 'patchType'
        word patchType_;

public:

    typedef faPatch Patch;
    typedef calculatedFaePatchField<Type> Calculated;

    //- Runtime type information
    TypeName("faePatchField");

    //- Debug switch to disallow the use of genericFaePatchField
    static int disallowGenericFaePatchField;

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            faePatchField,
            patch,
            (
                const faPatch& p,
                const DimensionedField<Type, faEdgeMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faePatchField,
            patchMapper,
            (
                const faePatchField<Type>& ptf,
                const faPatch& p,
                const DimensionedField<Type, faEdgeMesh>& iF,
                const faPatchFieldMapper& m
            ),
            (dynamic_cast<const faePatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faePatchField,
            dictionary,
            (
                const faPatch& p,
                const DimensionedField<Type, faEdgeMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        faePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct from patch and internal field and patch field
        faePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        faePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given faePatchField onto a new patch
        faePatchField
        (
            const faePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        faePatchField(const faePatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<faePatchField<Type>> clone() const
        {
            return tmp<faePatchField<Type>>(new faePatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        faePatchField
        (
            const faePatchField<Type>&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faePatchField<Type>> clone
        (
            const DimensionedField<Type, faEdgeMesh>& iF
        ) const
        {
            return tmp<faePatchField<Type>>
            (
                new faePatchField<Type>(*this, iF)
            );
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<faePatchField<Type>> New
        (
            const word&,
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        //  Allows override of constraint type
        static tmp<faePatchField<Type>> New
        (
            const word&,
            const word& actualPatchType,
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given faePatchField mapped onto a new patch
        static tmp<faePatchField<Type>> New
        (
            const faePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<faePatchField<Type>> New
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFaePatchField created on
        //  freestore without setting patchField values
        static tmp<faePatchField<Type>> NewCalculatedType
        (
            const faPatch&
        );

        //- Return a pointer to a new calculatedFaePatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<faePatchField<Type>> NewCalculatedType
        (
            const faePatchField<Type2>&
        );


    // Destructor

        virtual ~faePatchField() = default;


    // Member functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const faPatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, faEdgeMesh>&
            internalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& primitiveField() const
            {
                return internalField_;
            }

            //- Optional patch type
            const word& patchType() const
            {
                return patchType_;
            }

            //- Optional patch type
            word& patchType()
            {
                return patchType_;
            }

            virtual word constraintType()
            {
                return word("");
            }

            //- Return the type of the calculated for of faePatchField
            static const word& calculatedType();

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const faPatchFieldMapper&
            );

            //- Reverse map the given faePatchField onto this faePatchField
            virtual void rmap
            (
                const faePatchField<Type>&,
                const labelList&
            );


        //- Write
        virtual void write(Ostream&) const;


        // Check

            //- Check faePatchField<Type> against given faePatchField<Type>
            void check(const faePatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const faePatchField<Type>&);
        virtual void operator+=(const faePatchField<Type>&);
        virtual void operator-=(const faePatchField<Type>&);
        virtual void operator*=(const faePatchField<scalar>&);
        virtual void operator/=(const faePatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void forceAssign(const faePatchField<Type>&);
        virtual void forceAssign(const Field<Type>&);
        virtual void forceAssign(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const faePatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "fields/faePatchFields/faePatchField/faePatchField.C"
#include "fields/faePatchFields/basic/calculated/calculatedFaePatchField.H"
#endif


#define makeFaePatchTypeFieldTypeName(type)                                \
                                                                           \
defineNamedTemplateTypeNameAndDebug(type, 0);

#define makeFaePatchFieldsTypeName(type)                                   \
                                                                           \
makeFaePatchTypeFieldTypeName(type##FaePatchScalarField);                  \
makeFaePatchTypeFieldTypeName(type##FaePatchVectorField);                  \
makeFaePatchTypeFieldTypeName(type##FaePatchSphericalTensorField);         \
makeFaePatchTypeFieldTypeName(type##FaePatchSymmTensorField);              \
makeFaePatchTypeFieldTypeName(type##FaePatchTensorField);

#define makeFaePatchTypeField(PatchTypeField, typePatchTypeField)          \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField, typePatchTypeField, patch                              \
);                                                                         \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField,                                                        \
    typePatchTypeField,                                                    \
    patchMapper                                                            \
);                                                                         \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField, typePatchTypeField, dictionary                         \
);


#define makeFaePatchFields(type)                                           \
                                                                           \
makeFaePatchTypeField(faePatchScalarField, type##FaePatchScalarField);     \
makeFaePatchTypeField(faePatchVectorField, type##FaePatchVectorField);     \
makeFaePatchTypeField                                                      \
(                                                                          \
faePatchSphericalTensorField,                                                 \
type##FaePatchSphericalTensorField                                         \
);                                                                         \
makeFaePatchTypeField(faePatchSymmTensorField, type##FaePatchSymmTensorField); \
makeFaePatchTypeField(faePatchTensorField, type##FaePatchTensorField);


#define makeFaePatchTypeFieldTypedefs(type)                                \
                                                                           \
typedef type##FaePatchField<scalar> type##FaePatchScalarField;             \
typedef type##FaePatchField<vector> type##FaePatchVectorField;             \
typedef type##FaePatchField<sphericalTensor>                               \
type##FaePatchSphericalTensorField;                                        \
typedef type##FaePatchField<symmTensor> type##FaePatchSymmTensorField;     \
typedef type##FaePatchField<tensor> type##FaePatchTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

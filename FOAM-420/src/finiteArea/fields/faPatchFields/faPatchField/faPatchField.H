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
    (c) 2016 Esi Ltd.
    (c) 2016-2017 Wikki Ltd

Class
    Foam::faPatchField

Description
    faPatchField<Type> abstract base class.  This class gives a fat-interface
    to all derived classes covering all possible ways in which they might be
    used.  The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.  The next
    level of derivation covers all the specialised typed with specific
    evaluation proceedures, particularly with respect to specific fields.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faPatchField.C
    newPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef faPatchField_H
#define faPatchField_H

#include "faMesh/faPatches/faPatch/faPatch.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

class objectRegistry;
class dictionary;
class faPatchFieldMapper;
class areaMesh;


// * * * * * * Forward declaration of template friend fuctions * * * * * * * //
template<class Type>
class faMatrix;

template<class Type>
class faPatchField;

template<class Type>
class calculatedFaPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const faPatchField<Type>&);

/*---------------------------------------------------------------------------*\
                         Class faPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faPatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to a patch
        const faPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, areaMesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;

        //- Update index used so that manipulateMatrix is called only once
        //  during the construction of the matrix
        bool manipulatedMatrix_;

        //- Optional patch type, used to allow specified boundary conditions
        //  to be applied to constraint patches by providing the constraint
        //  patch type as 'patchType'
        word patchType_;



public:

    typedef faPatch Patch;
    typedef calculatedFaPatchField<Type> Calculated;

    //- Runtime type information
    TypeName("faPatchField");

    //- Debug switch to disallow the use of
    static int disallowGenericFaPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            patch,
            (
                const faPatch& p,
                const DimensionedField<Type, areaMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            patchMapper,
            (
                const faPatchField<Type>& ptf,
                const faPatch& p,
                const DimensionedField<Type, areaMesh>& iF,
                const faPatchFieldMapper& m
            ),
            (dynamic_cast<const faPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            dictionary,
            (
                const faPatch& p,
                const DimensionedField<Type, areaMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch and internal field and patch type
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const word& patchType
        );

        //- Construct from patch and internal field and patch field
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&,
            const bool valueRequired=false
        );

        //- Construct by mapping the given faPatchField onto a new patch
        faPatchField
        (
            const faPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        faPatchField(const faPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return tmp<faPatchField<Type>>(new faPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        faPatchField
        (
            const faPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type>>(new faPatchField<Type>(*this, iF));
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values).
        //  Allows override of constraint type
        static tmp<faPatchField<Type>> New
        (
            const word&,
            const word& actualPatchType,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<faPatchField<Type>> New
        (
            const word&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given faPatchField mapped onto a new patch
        static tmp<faPatchField<Type>> New
        (
            const faPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<faPatchField<Type>> New
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFaPatchField created on
        //  freestore without setting patchField values
//        template<class Type2>
        static tmp<faPatchField<Type>> NewCalculatedType
        (
            const faPatch&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<faPatchField<Type>> NewCalculatedType
        (
            const faPatchField<Type2>&
        );



    // Destructor

        virtual ~faPatchField() = default;


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
            const DimensionedField<Type, areaMesh>&
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

            //- Return the type of the calculated for of faPatchField
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

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
            }

            //- Return true if the matrix has already been manipulated
            bool manipulatedMatrix() const
            {
                return manipulatedMatrix_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const faPatchFieldMapper&
            );

            //- Reverse map the given faPatchField onto this faPatchField
            virtual void rmap
            (
                const faPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<Type>> snGrad() const;

            //- Return patch-normal gradient for coupled-patches
            //  using the deltaCoeffs provided
            virtual tmp<Field<Type>> snGrad
            (
                const scalarField& deltaCoeffs
            ) const
            {
                NotImplemented;
            }


            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();

            //- Update the coefficients associated with the patch field
            //  and apply weight field
            //  Sets Updated to true
            virtual void updateCoeffs(const scalarField& weights);

            //- Return internal field next to patch as patch field
            virtual tmp<Field<Type>> patchInternalField() const;

            //- Return internal field next to patch as patch field
            virtual void patchInternalField(Field<Type>&) const;

            //- Return patchField on the opposite patch of a coupled patch
            virtual tmp<Field<Type>> patchNeighbourField() const
            {
                NotImplemented;
            }

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            )
            {}

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );


            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<Field<scalar>>&
            ) const
            {
                NotImplemented;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<Field<scalar>>&
            ) const
            {
                NotImplemented;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const
            {
                NotImplemented;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this coupled patchField
            //  using the deltaCoeffs provided
            virtual tmp<Field<Type>> gradientInternalCoeffs
            (
                const scalarField& deltaCoeffs
            ) const
            {
                NotImplemented;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const
            {
                NotImplemented;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this coupled patchField
            //  using the deltaCoeffs provided
            virtual tmp<Field<Type>> gradientBoundaryCoeffs
            (
                const scalarField& deltaCoeffs
            ) const
            {
                NotImplemented;
            }

            //- Manipulate matrix
            virtual void manipulateMatrix(faMatrix<Type>& matrix);

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                faMatrix<Type>& matrix,
                const scalarField& weights
            );


        //- Write
        virtual void write(Ostream&) const;

        //- Helper function to write the keyword and entry only if the
        //  values are not equal. The value is then output as value2
        template<class EntryType>
        void writeEntryIfDifferent
        (
            Ostream& os,
            const word& entryName,
            const EntryType& value1,
            const EntryType& value2
        ) const;

        // Check

            //- Check faPatchField<Type> against given faPatchField<Type>
            void check(const faPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const faPatchField<Type>&);
        virtual void operator+=(const faPatchField<Type>&);
        virtual void operator-=(const faPatchField<Type>&);
        virtual void operator*=(const faPatchField<scalar>&);
        virtual void operator/=(const faPatchField<scalar>&);

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

        virtual void forceAssign(const faPatchField<Type>&);
        virtual void forceAssign(const Field<Type>&);
        virtual void forceAssign(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const faPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "fields/faPatchFields/faPatchField/faPatchField.C"
#include "fields/faPatchFields/basic/calculated/calculatedFaPatchField.H"
#endif

#define makeFaPatchField(faPatchTypeField)                                    \
                                                                              \
defineNamedTemplateTypeNameAndDebug(faPatchTypeField, 0);                     \
template<>                                                                    \
int faPatchTypeField::disallowGenericFaPatchField                             \
(                                                                             \
    debug::debugSwitch("disallowGenericFaPatchField", 0)                      \
);                                                                            \
defineTemplateRunTimeSelectionTable(faPatchTypeField, patch);                 \
defineTemplateRunTimeSelectionTable(faPatchTypeField, patchMapper);           \
defineTemplateRunTimeSelectionTable(faPatchTypeField, dictionary);


#define addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)   \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        PatchTypeField,                                                       \
        typePatchTypeField,                                                   \
        patch                                                                 \
    );                                                                        \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        PatchTypeField,                                                       \
        typePatchTypeField,                                                   \
        patchMapper                                                           \
    );                                                                        \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        PatchTypeField,                                                       \
        typePatchTypeField,                                                   \
        dictionary                                                            \
    );


// for non-templated patch fields
#define makeFaPatchTypeField(PatchTypeField, typePatchTypeField)              \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)

// for non-templated patch fields - use with caution
#define makeFaRemovablePatchTypeField(PatchTypeField, typePatchTypeField)     \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addRemovableToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField

// for templated patch fields
#define makeFaTemplatePatchTypeField(PatchTypeField, typePatchTypeField)      \
    defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);               \
    addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


#define makeFaPatchFields(type)                                               \
    makeFaTemplatePatchTypeField                                              \
    (                                                                         \
        faPatchScalarField,                                                   \
        type##FaPatchScalarField                                              \
    );                                                                        \
    makeFaTemplatePatchTypeField                                              \
    (                                                                         \
        faPatchVectorField,                                                   \
        type##FaPatchVectorField                                              \
    );                                                                        \
    makeFaTemplatePatchTypeField                                              \
    (                                                                         \
        faPatchSphericalTensorField,                                          \
        type##FaPatchSphericalTensorField                                     \
    );                                                                        \
    makeFaTemplatePatchTypeField                                              \
    (                                                                         \
        faPatchSymmTensorField,                                               \
        type##FaPatchSymmTensorField                                          \
    );                                                                        \
    makeFaTemplatePatchTypeField                                              \
    (                                                                         \
        faPatchTensorField,                                                   \
        type##FaPatchTensorField                                              \
    );


#define makeFaPatchTypeFieldTypeName(type)                                    \
defineNamedTemplateTypeNameAndDebug(type, 0);

#define makeFaPatchFieldsTypeName(type)                                       \
    defineNamedTemplateTypeNameAndDebug(type##FaPatchScalarField, 0);         \
    defineNamedTemplateTypeNameAndDebug(type##FaPatchVectorField, 0);         \
    defineNamedTemplateTypeNameAndDebug(type##FaPatchSphericalTensorField, 0);\
    defineNamedTemplateTypeNameAndDebug(type##FaPatchSymmTensorField, 0);     \
    defineNamedTemplateTypeNameAndDebug(type##FaPatchTensorField, 0)

#define makeFaPatchTypeFieldTypedefs(type)                                    \
    typedef type##FaPatchField<scalar> type##FaPatchScalarField;              \
    typedef type##FaPatchField<vector> type##FaPatchVectorField;              \
    typedef type##FaPatchField<sphericalTensor>                               \
        type##FaPatchSphericalTensorField;                                    \
    typedef type##FaPatchField<symmTensor> type##FaPatchSymmTensorField;      \
    typedef type##FaPatchField<tensor> type##FaPatchTensorField;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

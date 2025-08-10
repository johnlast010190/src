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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2013 H. Jasak
    (c) 2010-2023 Esi Ltd.

Class
    Foam::fvPatchField

Description
    Abstract base class with a fat-interface to all derived classes
    covering all possible ways in which they might be used.

    The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.

    The next level of derivation covers all the specialised types with
    specific evaluation procedures, particularly with respect to specific
    fields.

SourceFiles
    fvPatchField.C
    fvPatchFieldNew.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatchField_H
#define fvPatchField_H

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "volMesh/volMesh.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "fields/GeometricFields/GeometricField/GeometricFields.H"
#include "fields/Fields/fieldTypes.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/CoeffField/coeffFields.H"
#include "blockLduSystem/BlockLduSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class objectRegistry;
class dictionary;
class fvPatchFieldMapper;
class volMesh;
class gibFvPatchFieldMapper;


// Forward declaration of friend functions and operators

template<class Type>
class fvPatchField;

template<class Type>
class fvsPatchField;

template<class Type>
class calculatedFvPatchField;

template<class Type>
class fvMatrix;

template<class Type>
class fvBlockMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const fvPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                         Class fvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvPatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, volMesh>& internalField_;

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

        //- copy of the old field (needed by GIB)
        autoPtr<Field<Type>> oldTimeFieldPtr_;

        //- default value for the boundaries
        Type defaultGIBValue_;

public:

    typedef fvPatch Patch;
    typedef calculatedFvPatchField<Type> Calculated;


    //- Runtime type information
    TypeName("fvPatchField");

    //- Debug switch to disallow the use of genericFvPatchField
    static int disallowGenericFvPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patch,
            (
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patchMapper,
            (
                const fvPatchField<Type>& ptf,
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF,
                const fvPatchFieldMapper& m
            ),
            (dynamic_cast<const fvPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            dictionary,
            (
                const fvPatch& p,
                const DimensionedField<Type, volMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and value
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Type& value
        );

        //- Construct from patch and internal field and patch type
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const word& patchType
        );

        //- Construct from patch and internal field and patch field
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        fvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping the given fvPatchField onto a new patch
        fvPatchField
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fvPatchField(const fvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>(new fvPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        fvPatchField
        (
            const fvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>(new fvPatchField<Type>(*this, iF));
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<fvPatchField<Type>> New
        (
            const word&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values).
        //  Allows override of constraint type
        static tmp<fvPatchField<Type>> New
        (
            const word&,
            const word& actualPatchType,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given fvPatchField mapped onto a new patch
        static tmp<fvPatchField<Type>> New
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<fvPatchField<Type>> New
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        static tmp<fvPatchField<Type>> NewCalculatedType
        (
            const fvPatch&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<fvPatchField<Type>> NewCalculatedType
        (
            const fvPatchField<Type2>&
        );


    //- Destructor
    virtual ~fvPatchField() = default;


    // Member functions

        // Attributes

            //- Return the type of the calculated for of fvPatchField
            static const word& calculatedType();

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- for boundaries of blended type only
            //  provide a dummy implementation
            virtual List<bool> fixesValues() const
            {
                NotImplemented;
            }

            //- Return true if the value of the patch field
            //  is altered by assignment (the default)
            virtual bool assignable() const
            {
                return true;
            }

            //- for boundaries of blended type only
            //  provide a dummy implementation
            virtual List<bool> assignables() const
            {
                NotImplemented;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

            //- Return true if this patch field is region-coupled (coupled
            //  should be false)
            virtual bool regionCoupled() const
            {
                return false;
            }

            //- Zero 2nd derivative
            virtual bool extrapolated() const
            {
                return false;
            }

            //- Return true if this overrides the underlying constraint type
            bool overridesConstraint() const
            {
                if (type() == patch_.type())
                {
                    return false;
                }

                typename patchConstructorTable::iterator patchTypeCstrIter =
                    patchConstructorTable_().find(patch_.type());

                return patchTypeCstrIter != patchConstructorTable_().end();
            }


        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, volMesh>& internalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& primitiveField() const
            {
                return internalField_;
            }

            //- Return old field reference
            const Field<Type>& oldTimeField() const
            {
                return oldTimeFieldPtr_();
            }

            const Type& defaultGIBValue() const
            {
                return defaultGIBValue_;
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

            //- Lookup and return the patchField of the named field from the
            //  local objectRegistry.
            template<class GeometricField, class Type2>
            const typename GeometricField::Patch& lookupPatchField
            (
                const word& name
            ) const
            {
                return
                    patch().template lookupPatchFieldInDb
                    <
                        GeometricField,
                        Type2
                    >
                    (
                        db(),
                        name
                    );
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );

            //- Store old time and reset size
            virtual void updateGIB();

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);

            //- store the previous time of GIB
            virtual void storeGIB();

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

            //  Sets Updated to false
            virtual void resetUpdate();

            //- Return internal field next to patch as patch field
            virtual tmp<Field<Type>> patchInternalField() const;

            //- Return internal field next to patch as patch field
            virtual void patchInternalField(Field<Type>&) const;

            //- Initialise return of the patchField on the opposite patch of a
            //  coupled patch
            virtual void initPatchNeighbourField
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            ) const
            {}

            //- Return patchField on the opposite patch of a coupled patch
            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            ) const
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
            //  divergence operator in continuity - fluxes
            //  The operators from the above are used in momentum
            //  We need to use different coeffs contributions in some BCs
            virtual tmp<Field<Type>> valueDivInternalCoeffs
            (
                const tmp<Field<scalar>>& w
            ) const
            {
                return this->valueInternalCoeffs(w);
            }

            //- Return the matrix source coefficients corresponding to the
            //  divergence operator in continuity - fluxes
            //  The operators from the above are used in momentum
            //  We need to use different coeffs contributions in some BCs
            virtual tmp<Field<Type>> valueDivBoundaryCoeffs
            (
                const tmp<Field<scalar>>& w
            ) const
            {
                return this->valueBoundaryCoeffs(w);
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

            virtual tmp<CoeffField<Type>> gradientInternalBCoeffs() const;

            virtual tmp<Field<Type>> gradientBoundaryBCoeffs() const
            {
                return gradientBoundaryCoeffs();
            }

            virtual tmp<CoeffField<Type>> gradientInternalBCoeffs
            (
                const scalarField&
            ) const;

            virtual tmp<Field<Type>> gradientBoundaryBCoeffs
            (
                const scalarField& tDiff
            ) const
            {
                //- return RHS (source) segregated Coeffs
                return gradientBoundaryCoeffs(tDiff);
            }

            //- Return patch value for explicit gradient evaluation
            virtual tmp<Field<Type>> gradientBoundaryValue() const
            {
                return *this;
            }

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);

            //- Relax block matrix at the patch cells
            virtual void boundaryRelaxMatrix(fvBlockMatrix<Type>&) const;

            //- Operators specific to NS-eq.
                //- Momentum - gradp BC couplings (like total-pressure BC)
                virtual void addMomentumGradPCoupledBC
                (
                    fvBlockMatrix<Type>&,
                    const volScalarField&
                );

                //- Momentum - gradp BC couplings (like total-pressure BC)
                virtual void addContinuityCoupledBC
                (
                    BlockLduSystem<vector, scalar>&,
                    const volScalarField&,
                    const surfaceTensorField&
                ) const;

            //- Return fvOption-supplied boundary flux source
            virtual tmp<Field<Type>> boundarySources() const;

            //- Return fvOption-supplied boundary flux source using supplied
            // patch field values, and componentwise derivatives with respect
            // to patch field values
            virtual tmp<Field<Type>> boundarySources
            (
                const Field<Type>& f,
                Field<Type>& df
            ) const;

            //- Static function containing the actual implementation to avoid
            // circular include issues with fvOptions
            // This is a function template to allow partial specialisation
            template<class Type2>
            static tmp<Field<Type2>>
            boundarySources
            (
                const fvPatchField<Type2>& pf,
                const Field<Type2>& f,
                Field<Type2>& df
            );

            // Partial specialisations of above

            template<int N>
            static tmp<Field<VectorN<scalar,N>>>
            boundarySources
            (
                const fvPatchField<VectorN<scalar,N>>& pf,
                const Field<VectorN<scalar,N>>& f,
                Field<VectorN<scalar,N>>& df
            );

            template<int N>
            static tmp<Field<TensorN<scalar,N>>>
            boundarySources
            (
                const fvPatchField<TensorN<scalar,N>>& pf,
                const Field<TensorN<scalar,N>>& f,
                Field<TensorN<scalar,N>>& df
            );

            template<int N>
            static tmp<Field<DiagTensorN<scalar,N>>>
            boundarySources
            (
                const fvPatchField<DiagTensorN<scalar,N>>& pf,
                const Field<DiagTensorN<scalar,N>>& f,
                Field<DiagTensorN<scalar,N>>& df
            );

            template<int N>
            static tmp<Field<SphericalTensorN<scalar,N>>>
            boundarySources
            (
                const fvPatchField<SphericalTensorN<scalar,N>>& pf,
                const Field<SphericalTensorN<scalar,N>>& f,
                Field<SphericalTensorN<scalar,N>>& df
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;

            //- Helper function to write the keyword and entry only if the
            //  values are not equal. The value is then output as value2
            template<class EntryType>
            static void writeEntryIfDifferent
            (
                Ostream& os,
                const word& entryName,
                const EntryType& value1,
                const EntryType& value2
            );


        // Check

            //- Check fvPatchField<Type> against given fvPatchField<Type>
            void check(const fvPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvPatchField<Type>&);
        virtual void operator+=(const fvPatchField<Type>&);
        virtual void operator-=(const fvPatchField<Type>&);
        virtual void operator*=(const fvPatchField<scalar>&);
        virtual void operator/=(const fvPatchField<scalar>&);

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

        virtual void forceAssign(const fvPatchField<Type>&);
        virtual void forceAssign(const Field<Type>&);
        virtual void forceAssign(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const fvPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvPatchFields/fvPatchField/fvPatchFieldTemplates.C"
#endif


#define makeFvPatchField(fvPatchTypeField)                                     \
    defineNamedTemplateTypeNameAndDebug(fvPatchTypeField, 0);                  \
    template<>                                                                 \
    int fvPatchTypeField::disallowGenericFvPatchField                          \
    (                                                                          \
        debug::debugSwitch("disallowGenericFvPatchField", 0)                   \
    );                                                                         \
    defineTemplateRunTimeSelectionTable(fvPatchTypeField, patch);              \
    defineTemplateRunTimeSelectionTable(fvPatchTypeField, patchMapper);        \
    defineTemplateRunTimeSelectionTable(fvPatchTypeField, dictionary);


#define addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)    \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch                                                                  \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


#define addNamedToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField, name)\
    addNamedToRunTimeSelectionTable                                            \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch,                                                                 \
        name                                                                   \
    );                                                                         \
    addNamedToRunTimeSelectionTable                                            \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper,                                                           \
        name                                                                   \
    );                                                                         \
    addNamedToRunTimeSelectionTable                                            \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary,                                                            \
        name                                                                   \
    );


#define addSpecialNamedToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField, name, nameSanitised)\
    addSpecialNamedToRunTimeSelectionTable                                     \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch,                                                                 \
        name,                                                                  \
        nameSanitised                                                          \
    );                                                                         \
    addSpecialNamedToRunTimeSelectionTable                                     \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper,                                                           \
        name,                                                                  \
        nameSanitised                                                          \
    );                                                                         \
    addSpecialNamedToRunTimeSelectionTable                                     \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary,                                                            \
        name,                                                                  \
        nameSanitised                                                          \
    );


// Use with caution
#define addRemovableToPatchFieldRunTimeSelection(                              \
    PatchTypeField, typePatchTypeField)                                        \
                                                                               \
    addRemovableToRunTimeSelectionTable                                        \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch                                                                  \
    );                                                                         \
    addRemovableToRunTimeSelectionTable                                        \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addRemovableToRunTimeSelectionTable                                        \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


// For non-templated patch fields
#define makePatchTypeField(PatchTypeField, typePatchTypeField)                 \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


// For non-templated patch fields - use with caution
#define makeRemovablePatchTypeField(PatchTypeField, typePatchTypeField)        \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addRemovableToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


// For templated patch fields - explicitly instantiate
#define makeTemplatePatchTypeField(fieldType, type)                            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field),                     \
        0                                                                      \
    );                                                                         \
    addToPatchFieldRunTimeSelection                                            \
    (                                                                          \
        CAT3(fvPatch, CAPITALIZE(fieldType), Field),                           \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field)                      \
    )                                                                          \
    template class CAT3(type, FvPatch, Field)<fieldType>;


#define makePatchFields(type)                                                  \
    FOR_ALL_FIELD_TYPES(makeTemplatePatchTypeField, type)


#define makePatchFieldTypeName(fieldType, type)                                \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field),                     \
        0                                                                      \
    );                                                                         \
    template class CAT3(type, FvPatch, Field)<fieldType>;


#define makePatchFieldTypeNames(type)                                          \
    FOR_ALL_FIELD_TYPES(makePatchFieldTypeName, type)


#define makePatchTypeFieldTypedef(fieldType, type)                             \
    typedef type##FvPatchField<fieldType>                                      \
        CAT4(type, FvPatch, CAPITALIZE(fieldType), Field);


#define makePatchTypeFieldTypedefs(type)                                       \
    FOR_ALL_FIELD_TYPES(makePatchTypeFieldTypedef, type)


#define makePatchTypeFieldTypeName(type)                                       \
    defineNamedTemplateTypeNameAndDebug(type, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

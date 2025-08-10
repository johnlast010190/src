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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::GeometricBoundaryField

Description
    Generic GeometricBoundaryField class.

SourceFiles
    GeometricBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef GeometricBoundaryField_H
#define GeometricBoundaryField_H

#include "dimensionedTypes/dimensionedTypes.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "fields/FieldFields/FieldField/FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

// Forward declaration of friend functions and operators

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

/*---------------------------------------------------------------------------*\
                    Class GeometricBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricBoundaryField
:
    public FieldField<PatchField, Type>
{
public:

    // Public typedefs

        //- Type of boundary mesh on which this boundary is instantiated
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        //- Type of the internal field from which this GeometricField
        //  is derived
        typedef DimensionedField<Type, GeoMesh> Internal;


private:

    // Private data

        //- Reference to BoundaryMesh for which this field is defined
        const BoundaryMesh& bmesh_;


public:

    // Constructors

        //- Construct from a BoundaryMesh
        GeometricBoundaryField(const BoundaryMesh&);

        //- Construct from a BoundaryMesh, reference to the internal field
        //  and a patch field type
        GeometricBoundaryField
        (
            const BoundaryMesh&,
            const Internal&,
            const word&,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct from a BoundaryMesh, reference to the internal field,
        //  a wordList of patch field types and optionally
        //  the actual patch types (to override constraint patches)
        GeometricBoundaryField
        (
            const BoundaryMesh&,
            const Internal&,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct from a BoundaryMesh, reference to the internal field
        //  and a PtrList<PatchField<Type>>
        GeometricBoundaryField
        (
            const BoundaryMesh&,
            const Internal&,
            const PtrList<PatchField<Type>>&
        );

        //- Construct as copy setting the reference to the internal field
        GeometricBoundaryField(const Internal&, const GeometricBoundaryField&);

        //- Construct as copy
        //  Dangerous because Field may be set to a field which gets deleted
        //  Need new type of BoundaryField, one which is part of a geometric
        //  field for which snGrad etc. may be called and a free standing
        //  BoundaryField for which such operations are unavailable.
        GeometricBoundaryField(const GeometricBoundaryField&);

        //- Construct from dictionary
        GeometricBoundaryField
        (
            const BoundaryMesh&,
            const Internal&,
            const dictionary&
        );


    // Member Functions

        //- Read the boundary field
        void readField(const Internal& field, const dictionary& dict);

        //- Update the boundary condition coefficients
        void updateCoeffs();

        //- Update the values stored for patches by re-copying internal data
        void updatePatchFields();

        //- Evaluate boundary conditions
        void evaluate();

        //- Evaluate boundary conditions - static helper which can operate
        //  on a FieldField
        static void evaluatePatchFields
        (
            const typename GeoMesh::Mesh& mesh,
            PtrList<PatchField<Type>>& ff
        );

        //- Reset boundary field using provided boundary field,
        //  maintaining consistent reference to the internal field
        void reset(const GeometricBoundaryField&);

        //- Return a list of the patch field types
        wordList types() const;

        //- Return a list of the actual patch types
        wordList patchTypes() const;

        //- Combine two lists of patch types
        wordList combinePatchTypes(const wordList&, const wordList&) const;

        //- Return BoundaryField of the cell values neighbouring the boundary
        tmp<GeometricBoundaryField> boundaryInternalField() const;

        //- Return BoundaryField of the values on the other side of couples
        tmp<GeometricBoundaryField> boundaryNeighbourField() const;

        //- Return a list of pointers for each patch field with only those
        //  pointing to interfaces being set
        LduInterfaceFieldPtrsList<Type> interfaces() const;

        //- Return a list of pointers for each patch field with only those
        //  pointing to interfaces being set
        lduInterfaceFieldPtrsList scalarInterfaces() const;

        //- Return a list of pointers for each patch field with only those
        //  pointing to block-coupled interfaces being set
        typename BlockLduInterfaceFieldPtrsList<Type>::Type
        blockInterfaces() const;

        //- Return a copy of the old time field
        List<Field<Type>> oldTimeField() const;

        //- Write boundary field as dictionary entry
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Write dictionary entries of the individual boundary fields
        void writeEntries(Ostream& os) const;


    // Member Operators

        //- Assignment operator
        void operator=(const GeometricBoundaryField&);

        //- Assignment to FieldField<PatchField, Type>
        void operator=(const FieldField<PatchField, Type>&);

        //- Assignment to FieldField<OtherPatchField, Type>
        template<template<class> class OtherPatchField>
        void operator=(const FieldField<OtherPatchField, Type>&);

        //- Assignment to Type
        void operator=(const Type&);


        //- Forced assignment to
        //  BoundaryField<Type, PatchField, BoundaryMesh>
        void forceAssign(const GeometricBoundaryField&);

        //- Forced assignment to FieldField<PatchField, Type>
        void forceAssign(const FieldField<PatchField, Type>&);

        //- Forced assignment to FieldField<OtherPatchField, Type>
        template<template<class> class OtherPatchField>
        void forceAssign(const FieldField<OtherPatchField, Type>&);

        //- Forced assignment to Type
        void forceAssign(const Type&);
};


template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/GeometricFields/GeometricField/GeometricBoundaryField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

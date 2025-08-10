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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd.

Class
    Foam::fvPatch

Description
    A finiteVolume patch using a polyPatch and a fvBoundaryMesh

SourceFiles
    fvPatch.C
    fvPatchNew.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatch_H
#define fvPatch_H

#include "meshes/polyMesh/polyPatches/directPolyPatch/directPolyPatch.H"
#include "primitives/ints/lists/labelList.H"
#include "containers/Lists/SubList/SubList.H"
#include "db/typeInfo/typeInfo.H"
#include "memory/tmp/tmp.H"
#include "fields/Fields/primitiveFields.H"
#include "fields/Fields/Field/SubField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldsFwd.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "db/objectRegistry/objectRegistry.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvBoundaryMesh;
class surfaceInterpolation;

/*---------------------------------------------------------------------------*\
                           Class fvPatch Declaration
\*---------------------------------------------------------------------------*/

class fvPatch
{
    // Private data

        //- Reference to the underlying polyPatch
        const polyPatch& polyPatch_;

        //- Reference to boundary mesh
        const fvBoundaryMesh& boundaryMesh_;

        bool active_;


    // Private Member Functions

        //- Disallow construct as copy
        fvPatch(const fvPatch&);

        //- Disallow assignment
        void operator=(const fvPatch&);


protected:

    // Protected Member Functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const;

        //- Initialise the patches for moving points
        virtual void initMovePoints();

        //- Correct patches after moving points
        virtual void movePoints();


public:

    typedef fvBoundaryMesh BoundaryMesh;

    friend class fvBoundaryMesh;
    friend class surfaceInterpolation;

    //- Runtime type information
    TypeName(directPolyPatch::typeName_());


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvPatch,
            polyPatch,
            (const polyPatch& patch, const fvBoundaryMesh& bm),
            (patch, bm)
        );


    // Constructors

        //- Construct from polyPatch and fvBoundaryMesh
        fvPatch(const polyPatch&, const fvBoundaryMesh&);


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<fvPatch> New
        (
            const polyPatch&,
            const fvBoundaryMesh&
        );


    //- Destructor
    virtual ~fvPatch();


    // Member Functions

        // Access

            //- Return the polyPatch
            const polyPatch& patch() const
            {
                return polyPatch_;
            }

            //- Return name
            const word& name() const
            {
                return polyPatch_.name();
            }

            //- Return start label of this patch in the polyMesh face list
            virtual label start() const
            {
                return polyPatch_.start();
            }

            //- Return size
            virtual label size() const
            {
                return polyPatch_.size();
            }

            //- Return true if this patch is constrained
            virtual bool constrained() const
            {
                return false;
            }

            //- Return true if this patch is enables
            //- Used inside the GIB
            virtual bool isActive() const
            {
                return active_;
            }

            //- Function to change the status of the patch
            //- Used inside the GIB
            virtual bool& active()
            {
                return active_;
            }

            //- Return true if this patch is coupled
            virtual bool coupled() const
            {
                return polyPatch_.coupled();
            }

            //- Return true if the given type is a constraint type
            static bool constraintType(const word& pt);

            //- Return a list of all the constraint patch types
            static wordList constraintTypes();

            //- Return the index of this patch in the fvBoundaryMesh
            label index() const
            {
                return polyPatch_.index();
            }

            //- Return boundaryMesh reference
            const fvBoundaryMesh& boundaryMesh() const
            {
                return boundaryMesh_;
            }

            //- Slice list to patch
            template<class T>
            const typename List<T>::subList patchSlice(const List<T>& l) const
            {
                return typename List<T>::subList(l, size(), start());
            }

            //- Return faceCells
            virtual const labelUList& faceCells() const
            {
                return polyPatch_.faceCells();
            }


        // Access functions for geometrical data

            //- Return face centres
            const vectorField& Cf() const;

            //- Return neighbour cell centres
            tmp<vectorField> Cn() const;

            //- Return face area vectors
            const vectorField& Sf() const;

            //- Return face area magnitudes
            const scalarField& magSf() const;

            //- Return face normals
            tmp<vectorField> nf() const;

            //- Return cell-centre to face-centre vector
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre vector is returned
            virtual tmp<vectorField> delta() const;


        // Access functions for demand driven data

            //- Return patch weighting factors
            const scalarField& weights() const;

            //- Return the face - cell distance coeffient
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre distance coeffient is returned
            const scalarField& deltaCoeffs() const;


        // Evaluation functions

            //- Return given internal field next to patch as patch field
            template<class Type>
            tmp<Field<Type>> patchInternalField(const UList<Type>&) const;

            //- Return given internal field next to patch as patch field
            template<class Type>
            void patchInternalField(const UList<Type>&, Field<Type>&) const;

            //- Return the corresponding patchField of the named field
            template<class GeometricField, class Type>
            const typename GeometricField::Patch& patchField
            (
                const GeometricField&
            ) const;

            //- Lookup and return the patchField of the named field from the
            //  local objectRegistry.
            //  N.B.  The dummy pointer arguments are used if this function is
            //  instantiated within a templated function to avoid a bug in gcc.
            //  See inletOutletFvPatchField.C and outletInletFvPatchField.C
            template<class GeometricField, class Type>
            const typename GeometricField::Patch& lookupPatchField
            (
                const word& name,
                const GeometricField* = nullptr,
                const Type* = nullptr
            ) const;

            template<class GeometricField, class Type>
            const typename GeometricField::Patch& lookupPatchFieldInDb
            (
            const objectRegistry& obr,
            const word& name,
            const GeometricField* = nullptr,
            const Type* = nullptr
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh/fvPatches/fvPatch/fvPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

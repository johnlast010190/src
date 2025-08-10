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
    (c) 2016-2023 Esi Ltd.

Class
    Foam::regionState

Description

    Intermediate abstract base class for regioned state functions. Provides
    mesh modification and field instantiation and initialisation capability.
    and fields

SourceFiles
    regionState.C
    regionStateTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef regionState_H
#define regionState_H

#include "stateFunction/stateFunction.H"
#include "fieldInitialisation/fieldInit/fieldInit.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                         Class regionState Declaration
\*---------------------------------------------------------------------------*/

class regionState
:
    public stateFunction
{
    // Private data

        //- Region mesh
        fvMesh* meshPtr_;

        autoPtr<fvSolutionRegistry> solutionRegPtr_;

        //- Boundary type modification components

            //- stored field names and internal values
            //- needed for field maintenance during boundary type modification
            PtrList<Tuple2<word, Field<scalar>>> scalarFields_;
            PtrList<Tuple2<word, Field<vector>>> vectorFields_;
            PtrList<Tuple2<word, Field<tensor>>> tensorFields_;
            PtrList<Tuple2<word, Field<symmTensor>>> symmTensorFields_;
            PtrList<Tuple2<word, Field<sphericalTensor>>>
                sphericalTensorFields_;

            //- stored field boundary dictionaries
            //- needed for field maintenance during boundary type modification
            autoPtr<dictionary> initialFieldBoundaryDictsPtr_;



    // Private Member Functions

        //- Create volume fields with different types
        template<class Type>
        void createVolField
        (
            const fvMesh& mesh,
            const dictionary& fieldDict,
            const word fieldName
        );

        //- Create surface fields with different types
        template<class Type>
        void createSurfaceField
        (
            const surfaceMesh& mesh,
            const dictionary& fieldDict,
            const word fieldName
        );

        //- Create point fields with different types
        template<class Type>
        void createPointField
        (
            const pointMesh& mesh,
            const dictionary& fieldDict,
            const word fieldName
        );

        //- Create fields with different variable types
        template<class Type>
        void createTypeField
        (
            const dictionary& fd,
            const word fieldName
        );


        //- Read and store existing fields
        template
        <
            class Type,
            template<class> class PatchField,
            class Mesh,
            class Mesh2
        >
        void readInitialGeoField
        (
            const objectRegistry& runTime,
            const Mesh2& mesh,
            const dictionary& fd,
            const word& fieldName
        );


        //- read all types of existing fields
        template<class Type>
        void readInitialTypeField
        (
            const dictionary& fd,
            const word fieldName
        );

        //- Re-do boundary initialisation after internal field was initialised
        void reInitBoundary(const fieldInit& f);

        template<template<class> class PatchType, class MeshType>
        void reInitBoundaryTypes(const fieldInit& f);

        template<class Type, template<class> class PatchType, class MeshType>
        void reInitBoundaryType(const fieldInit& f);

        //- Clear residual files from a previous setup with non-conformal meshes
        void clearNonConformalSetup(const bool initFields = false);

        //- Generate initial default patch dictionary
        Xfer<dictionary> initialPatchDict(const label& patchI) const;

        //- Correct handling of inlet/outlet sub-type patches
        word correctPatchPhysicalType(dictionary& patchDict);

        //- Check validity of multi-patch modification operations
        //- no renaming or restricted patch types allowed
        bool validMultiPatchModify
        (
            const dictionary& defaults,
            const dictionary& patchDict
        );

        //- Modify patch definition
        void modifyBoundaryPatch
        (
            const fvMesh& mesh,
            const dictionary& dict,
            label patchI
        );

        //- Modify boundary mesh names and types
        void modifyBoundaryMesh
        (
            const dictionary& defaults,
            List<bool>& modifiedPatch
        );

        //- store internal field and field name on local variables
        template<class Type>
        void storeInternalField
        (
            const word& name,
            const Field<Type>&
        );

        //- search storage for named field and reset if found
        template<class Type>
        bool setInternalValue(const word& name);

        //- read old fields for possible reuse
        void readInitialFields();

        //- instantiate fields
        void instantiateFields();

        //- update initialisation depending on switches
        void updateInitialisation(fieldInit&);

        //- Initialise fields
        void initialiseFields();

        //- clear fields not initialised by caseSetup
        void clearFieldFiles();

        //- clear fields from current instance (multi-instance cases)
        void clearInstanceFields();

        //- synchronise parallel operations
        void sync() const;

        //- Disallow default bitwise copy construct
        regionState(const regionState&);

        //- Disallow default bitwise assignment
        void operator=(const regionState&);


public:

    //- Runtime type information
        TypeName("region");

    // Static data members

    // Constructors

        //- Construct meshed region from components
        regionState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    // Selectors


    //- Destructor
    ~regionState() override;


    // Member Functions

        // Access

            //- const access to fvMesh
            inline const fvMesh& mesh() const
            {
                return *meshPtr_;
            }

            //- Return const reference to the local database
            inline const fvSolutionRegistry& localDb() const
            {
                return *solutionRegPtr_;
            }
        // Check

        // Edit

            //- Generate and output fields
            void createObjects(const Time& time) override;

            //- Modify mesh
            void modifyMesh(scalar writePause) override;

            //- Create new fields
            void createFields() override;


};


// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

template<>
void regionState::storeInternalField<scalar>
(
    const word& name,
    const Field<scalar>& v
);

template<>
void regionState::storeInternalField<vector>
(
    const word& name,
    const Field<vector>& v
);

template<>
void regionState::storeInternalField<tensor>
(
    const word& name,
    const Field<tensor>& v
);

template<>
void regionState::storeInternalField<symmTensor>
(
    const word& name,
    const Field<symmTensor>& v
);

template<>
void regionState::storeInternalField<sphericalTensor>
(
    const word& name,
    const Field<sphericalTensor>& v
);


template<>
bool regionState::setInternalValue<scalar>(const word& name);
template<>
bool regionState::setInternalValue<vector>(const word& name);
template<>
bool regionState::setInternalValue<tensor>(const word& name);
template<>
bool regionState::setInternalValue<symmTensor>(const word& name);
template<>
bool regionState::setInternalValue<sphericalTensor>(const word& name);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "regionStateTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

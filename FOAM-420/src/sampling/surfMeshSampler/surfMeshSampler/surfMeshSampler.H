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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::surfMeshSampler

Group
    grpUtilitiesFunctionObjects

Description
    An abstract class for surfMesh with sampling.

SourceFiles
    surfMeshSampler.C
    surfMeshSamplerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef surfMeshSampler_H
#define surfMeshSampler_H

#include "surfMesh/surfMesh.H"
#include "surfFields/surfFields/surfFields.H"

#include "meshes/primitiveShapes/point/pointField.H"
#include "primitives/strings/word/word.H"
#include "primitives/ints/lists/labelList.H"
#include "meshes/meshShapes/face/faceList.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "memory/autoPtr/autoPtr.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "meshes/polyMesh/polyMesh.H"
#include "interpolation/interpolation/interpolation/interpolation.H"
#include "db/error/error.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class surfMeshSampler;

/*---------------------------------------------------------------------------*\
                       Class surfMeshSampler Declaration
\*---------------------------------------------------------------------------*/

class surfMeshSampler
{
    // Private data

        //- Name for surfMesh lookup
        word name_;

        //- Mesh reference
        const polyMesh& mesh_;


    // Private Member Functions

    // Service methods


protected:

    // Protected Member Functions

        //- Get existing or create new surfMesh
        surfMesh& getOrCreateSurfMesh() const;

        //- Get existing or create new surfField.
        //  Create with same name and dimensions as the 'parent' volField.
        template<class Type>
        DimensionedField<Type, surfGeoMesh>&
        getOrCreateSurfField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;


        // //- Transfer field from object registry to surfField
        // template<class Type>
        // bool transferField
        // (
        //     const objectRegistry& store,
        //     const word& fieldName
        // ) const;


        //- Write the given fields
        template<class Type>
        label writeFields(const wordReList& select) const;

public:

    //- Runtime type information
    TypeName("surfMeshSampler");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        surfMeshSampler,
        word,
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        ),
        (name, mesh, dict)
    );


   // iNew helper class

        //- Class used for the PtrLists read-construction
        class iNew
        {
            //- Reference to the volume mesh
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<surfMeshSampler> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);

                return surfMeshSampler::New(name, mesh_, dict);
            }
        };


    // Constructors

        //- Construct from name, mesh
        surfMeshSampler
        (
            const word& name,
            const polyMesh& mesh
        );

        //- Construct from dictionary
        surfMeshSampler
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


        //- Clone
        autoPtr<surfMeshSampler> clone() const
        {
            NotImplemented;
        }



    // Selectors

        //- Return a reference to the selected surface
        static autoPtr<surfMeshSampler> New
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfMeshSampler();


    // Member Functions

      // Access

        //- Access to the underlying mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Name of surface
        const word& name() const
        {
            return name_;
        }


        //- Create surfMesh if required.
        void create() const;

        //- Return existing surfMesh.
        virtual const surfMesh& surface() const;


        //- Does the surface need an update?
        virtual bool needsUpdate() const = 0;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire() = 0;

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update() = 0;


      // Edit

        //- Rename
        virtual void rename(const word& newName)
        {
            name_ = newName;
        }


      // Execute

        //- Sample from specified volume field to surface
        virtual bool sample(const word& fieldName) const = 0;

        //- Sample from volume fields to specified surface fields.
        virtual label sample(const UList<word>& fields) const;


      // Write

        //- Write specified fields
        virtual label write(const wordReList& fieldSelection) const;

        //- Write
        virtual void print(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "surfMeshSampler/surfMeshSampler/surfMeshSamplerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

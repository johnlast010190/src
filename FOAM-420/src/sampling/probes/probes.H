/*---------------------------------------------------------------------------*\
e       o        |
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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2017-2022 Esi Ltd.

Class
    Foam::probes

Group
    grpUtilitiesFunctionObjects

Description
    Set of locations to sample.

    Call write() to sample and write files.

    Example of function object specification:
    \verbatim
    probes
    {
        type            probes;
        libs            ("libsampling.so");

        // Name of the directory for probe data
        name            probes;

        // Write at same frequency as fields
        writeControl    outputTime;
        writeInterval   1;

        // Fields to be probed
        fields          (p U);

        // Optional: do not recalculate cells if mesh moves
        fixedLocations  false;

        // Optional: interpolation scheme to use (default is cell)
        interpolationScheme cellPoint;

        probeLocations
        (
            ( 1e-06 0 0.01 )      // at inlet
            (0.21 -0.20999 0.01)  // at outlet1
            (0.21 0.20999 0.01)   // at outlet2
            (0.21 0 0.01)         // at central block
        );

        // Optional: filter out points that haven't been found. Default
        //           is to include them (with value -VGREAT)
        includeOutOfBounds  true;

        // Optional: Putting probes to reference frame
        // Reference frame name
        referenceFrame rotatingFrame;

        // Are the points defined in the frame?
        definedInFrame false;

        // Local coordinate system for output
        localOnOutput true;
    }
    \endverbatim

SourceFiles
    probes.C

\*---------------------------------------------------------------------------*/

#ifndef probes_H
#define probes_H

#include "db/functionObjects/regionFunctionObject/regionFunctionObject.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "surfaceMesh/surfaceMesh.H"
#include "primitives/strings/lists/wordReList.H"
#include "global/columnatedFileWriter/columnatedFileWriter.H"
#include "referenceFrames/coordinateFrame.H"

using namespace Foam::functionObjects;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class objectRegistry;
class dictionary;
class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class probes Declaration
\*---------------------------------------------------------------------------*/

class probes
:
    public regionFunctionObject,
    public pointField
{
protected:

    // Protected classes

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:
            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0)
            {}
        };


    // Protected data

        //- Copy of the dictionary
        // This is inelegant, but at least it doesn't go out of scope when read
        // by columnatedFileWriter.
        const dictionary dict_;

        //- Switch to determine whether or not to use the columnatedFileWriter
        const bool newFileFormat_ = false;

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        //- Flag if mesh updated when fixedLocations active
        bool meshUpdated_;

        //- Coordinate frame
        coordinateFrame* coorFramePtr_;

        //- Probe location reported in local coordinate frame?
        Switch localOnOutput_;

        //- Store time 0 points
        tmp<pointField> time0Points_;

      // Read from dictonary

        //- Names of fields to probe
        wordReList fieldSelection_;

        //- Fixed locations, default = yes
        //  Note: set to false for moving mesh calculations where locations
        //        should move with the mesh
        bool fixedLocations_;

        //- Interpolation scheme name
        //  Note: only possible when fixedLocations_ is true
        word interpolationScheme_;

        //- Include probes that were not found
        bool includeOutOfBounds_;

        //- if point is outside mesh, use nearest face/cell,
        bool nearestElement_;


      // Calculated

        //- Categorized scalar/vector/tensor vol fields
        fieldGroup<scalar> scalarFields_;
        fieldGroup<vector> vectorFields_;
        fieldGroup<sphericalTensor> sphericalTensorFields_;
        fieldGroup<symmTensor> symmTensorFields_;
        fieldGroup<tensor> tensorFields_;

        //- Categorized scalar/vector/tensor surf fields
        fieldGroup<scalar> surfaceScalarFields_;
        fieldGroup<vector> surfaceVectorFields_;
        fieldGroup<sphericalTensor> surfaceSphericalTensorFields_;
        fieldGroup<symmTensor> surfaceSymmTensorFields_;
        fieldGroup<tensor> surfaceTensorFields_;

        //- Cells to be probed (obtained from the locations)
        labelList elementList_;

        //- Faces to be probed
        labelList faceList_;

        //- Processor holding the cell or face (-1 if point not found
        //  on any processor)
        labelList processor_;

        //- Current open files. This is for the legacy output, which is
        //  strangely complicated.  I'm going to leave the legacy output and
        //  non-legacy outputs completely separate for now, to ensure I don't
        //  break backwards compatability.
        HashPtrTable<OFstream> probeFilePtrs_;

        //- columnatedFileWriter wraps an OFStream, which inherits from
        //  OFstreamAllocator, which has an old-fashioned raw pointer in it.  It
        //  therefore can't be copied, and the OpenFOAM containers require copy
        //  constructors.  Therefore, we have to use this HashPtrTable thing
        //  (which again, uses raw pointers).  This should be fixed properly,
        //  but that would involve making all the OpenFOAM containers work with
        //  move constructors and I don't have time right now.
        HashPtrTable<columnatedFileWriter, word> fileWriters_;


    // Protected Member Functions

        //- Clear old field groups
        void clearFieldGroups();

        //- Classify field types, returns the number of fields
        label classifyFields();

        //- Find cells and faces containing probes
        virtual void findElements(const fvMesh& mesh);

        //- Classify field type and open/close file streams,
        //  returns number of fields to sample
        label prepare();

        wordHashSet getAllFieldNames() const;
        const fileName getOutputDir() const;

        void initialiseLegacyOutputFiles
        (
            wordHashSet& currentFields,
            const fileName& probeDir
        );

        void initialiseOutputFiles
        (
            wordHashSet& currentFields,
            const fileName& probeDir
        );

private:

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );


        //- Sample and write a particular surface field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );

        //- Sample and write all the fields of the given type
        template<class Type>
        void sampleAndWrite(const fieldGroup<Type>&);

        //- Sample and write all the surface fields of the given type
        template<class Type>
        void sampleAndWriteSurfaceFields(const fieldGroup<Type>&);

        //- Disallow default bitwise copy construct
        probes(const probes&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const probes&) = delete;


public:

    //- Runtime type information
    TypeName("probes");


    // Constructors

        //- Construct from Time and dictionary
        probes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );

    //- Destructor
    virtual ~probes();


    // Member Functions

        //- Return names of fields to probe
        virtual const wordReList& fieldNames() const
        {
            return fieldSelection_;
        }

        //- Return locations to probe
        virtual const pointField& probeLocations() const
        {
            return *this;
        }

        //- Return location for probe i
        virtual const point& probe(const label i) const
        {
            return operator[](i);
        }

        //- Cells to be probed (obtained from the locations)
        const labelList& elements() const
        {
            return elementList_;
        }

        //- Read the probes
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state)
        {}

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Sample a single vol field on all sample locations
        template<class Type>
        tmp<Field<Type>> sample(const word& fieldName) const;

        //- Sample a single scalar field on all sample locations
        template<class Type>
        tmp<Field<Type>> sampleSurfaceFields(const word& fieldName) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "probes/probesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

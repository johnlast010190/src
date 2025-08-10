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
    (c) 2019 Esi Ltd.
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::functionObjects::fieldCalcs::surfaceFieldCalc

Group
    grpFieldFunctionObjects

Description
    Provides a 'face regionType' variant of the fieldCalcs function object.

    Perform a list user-specified operations on selected faces. Output
    calculated field values and a generate report of field on selected faces

Usage
    Examples of function object specification:
    \verbatim
    surfaceFieldCalc1
    {
        type            surfaceFieldCalc;
        libs            ("libfieldFunctionObjects.so");

        log             true;
        writeControl    writeTime;
        writeFields     true;

        surfaceFormat   vtk;
        regionType      sampledSurface;

        probeLocations ((0 0 0));

        sampledSurfaceDict
        {
            type sampledTriSurfaceMesh;
            surface downstream.stl;
            source cells;//boundaryFaces;//cells;
            interpolate false;
        }

        name  downstream;

        operations
        (
            {
                operation swirlAngle;
                fieldName swirlAngle;
                writeFields true;
                writeReport true;
                subsetDict
                (
                    {
                        type sphere;
                        centre (0 0 0);
                        radius 1;
                        inside true;
                    }

                    {
                        type box;
                        boundBox (0 0 0) (1 1 1);
                        inside true;
                    }
                );
            }
            {
                operation field;
                fieldName U;
            }
            {
                operation field;
                fieldName p;
            }
        );

    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: surfaceFieldCalc   | yes         |
        log          | write data to standard output | no    | no
        surfaceFormat | output value format    | no          |
        regionType   | face regionType: see below  | yes         |
        name         | name of face regionType if required  | no |
        operations    | operation to perform    | yes         |
        faceZoneFlip  | geometric face zone flip based on zero face  | no |
        probeLocations | output field values at given probeLocations | no
    \endtable

    Where \c regionType is defined by
    \plaintable
        faceZone     | requires a 'name' entry to specify the faceZone
        patch        | requires a 'name' entry to specify the patch
        surface      | requires a 'name' entry to specify the surfMesh
        sampledSurface | requires a 'sampledSurfaceDict' sub-dictionary
    \endplaintable

    The \c operations is one of:
    \plaintable
       field            | field
       swirlAngle       | swirl angle
       axialDeviation   | field axial deviation from mean (%)
    \endplaintable

Note
    - faces on empty patches get ignored
    - if the field is a volField the \c faceZone can only consist of boundary
      faces
    - Using \c surface:
      - The keyword %subRegion should not be used to select surfaces.
        Specify instead the regionType 'surface' and provide the surface name.
    - using \c sampledSurface:
        - not available for surface fields
        - if interpolate=true they use \c interpolationCellPoint
          otherwise they use cell values
        - each triangle in \c sampledSurface is logically only in one cell
          so interpolation will be wrong when triangles are larger than
          cells.  This can only happen for sampling on a \c triSurfaceMesh
        - take care when using isoSurfaces - these might have duplicate
          triangles and so integration might be wrong

See also
    Foam::fieldCalcs
    Foam::functionObject

SourceFiles
    surfaceFieldCalc.C
    surfaceFieldCalcTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceFieldCalc_H
#define functionObjects_surfaceFieldCalc_H

#include "fieldCalcs/fieldCalc/fieldCalc.H"
#include "primitives/enums/NamedEnum.H"
#include "meshedSurf/meshedSurf.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"
#include "fields/volFields/volFieldsFwd.H"
#include "surfFields/surfFields/surfFieldsFwd.H"
#include "meshedSurf/meshedSurfRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sampledSurface;
class surfaceWriter;

namespace functionObjects
{
namespace fieldCalcs
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFieldCalc Declaration
\*---------------------------------------------------------------------------*/

class surfaceFieldCalc
:
    public fieldCalc
{
public:

    // Public data types

        //- Region type enumeration
        enum regionTypes
        {
            stFaceZone,
            stPatch,
            stSurface,
            stSampledSurface
        };

        //- Region type names
        static const NamedEnum<regionTypes, 4> regionTypeNames_;

        //- Operation type enumeration
        enum calcTypes
        {
            ctField,                //!< field
            ctSwirlAngle,           //!< swirl angle
            ctBulkSwirlAngle,       //!< the bulk swirl angle
            ctAxialDeviation        //!< field deviation
        };

        //- Operation type names
        static const NamedEnum<calcTypes, 4> calcTypeNames_;

private:

    // Private Member Functions

        //- Set faces to evaluate based on a face zone
        void setProbeIndex(const pointField& Cf);

        //- Set faces to evaluate based on a face zone
        void setFaceZoneFaces();

        //- Set faces to evaluate based on a patch
        void setPatchFaces();

        //- Combine mesh faces and points from multiple processors
        void combineMeshGeometry
        (
            faceList& faces,
            pointField& points
        ) const;

        //- Combine surface faces and points from multiple processors
        void combineSurfaceGeometry
        (
            faceList& faces,
            pointField& points
        ) const;

        //- Calculate and return total area of the surfaceFieldCalc: sum(magSf)
        scalar totalArea() const;


protected:

    // Protected data

        //- Region type
        regionTypes regionType_;

        //List of opeartions to perform on surface
        PtrList<dictionary> operations_;

        //List of probe locations to interrogate surface over
        pointField probeLocations_;

        //List of probe surface indexes
        labelList probeIndex_;

        //List of writers for operations
        PtrList<writeFile> operationsWriter_;

        //List of writers for probes
        List<autoPtr<writeFile>> probesWriter_;

        //- Total area of the surfaceFieldCalc
        scalar totalArea_;

        //- Global number of faces
        label nFaces_;

        //- Per face whether active
        DynamicList<label> filterFaces_;

        // If operating on mesh faces (faceZone, patch)

            //- Local list of face IDs
            labelList faceId_;

            //- Local list of patch ID per face
            labelList facePatchId_;

            //- List representing the face flip map
            //  (false: use as-is, true: negate)
            boolList faceFlip_;

        //- Underlying sampledSurface (if operating on sampledSurface)
        autoPtr<sampledSurface> surfacePtr_;

        //- Surface writer
        autoPtr<surfaceWriter> surfaceWriterPtr_;

    // Protected Member Functions

        //- The volume mesh or surface registry being used
        const objectRegistry& obr() const;

        //- Return the local list of face IDs
        inline const labelList& faceId() const;

        //- Return the local list of patch ID per face
        inline const labelList& facePatch() const;

        //- Return the local true/false list representing the face flip map
        inline const boolList& faceFlip() const;

        //- True if the specified operation needs a surface Sf or Cf
        bool needsGeom() const;

        //- Initialise, e.g. face addressing
        void initialise(const dictionary& dict);

        //- Return true if the field name is known and a valid type
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Return field values by looking up field name
        template<class Type>
        tmp<Field<Type>> getFieldCalcs
        (
            const word& fieldName,
            const bool mustGet = false
        ) const;

        //- Filter a surface field according to faceIds
        template<class Type>
        tmp<Field<Type>> filterField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field
        ) const;

        //- Filter a volume field according to faceIds
        template<class Type>
        tmp<Field<Type>> filterField
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        ) const;

        //- Output file header information
        virtual void writeFileHeader(writeFile& writer) const;

        //- Output file header information
        virtual void writeProbesFileHeader(writeFile& writer) const;

        void calcGeomFilterFaces
        (
            const vectorField& Cf,
            const PtrList<dictionary>& geomFilters
        );

public:

    //- Run-time type information
    TypeName("surfaceFieldCalc");


    // Constructors

        //- Construct from name, Time and dictionary
        surfaceFieldCalc
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        surfaceFieldCalc
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceFieldCalc();


    // Public Member Functions

        //- Return the region type
        inline const regionTypes& regionType() const;

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate and write
        virtual bool write();

        //- Calculate surface expressions
        void calculate
        (
            const meshedSurfRef& surfToWrite
        );

    // surface field evaluation functions

        scalar integrate
        (
            const Field<scalar> values,
            const vectorField& Sf
        );

        scalar integrate
        (
            const Field<vector> values,
            const vectorField& Sf
        );

        template<class Type>
        void writeReport
        (
            const Field<Type> values,
            const word& fieldName,
            const vectorField& Sf,
            writeFile& writer,
            autoPtr<writeFile>& probesWriter
        );

        template<class Type>
        bool writeField
        (
            const Field<Type> values,
            const word& fieldName,
            const meshedSurf& surfToWrite,
            const writeFile& writer
        );

        //- Calculation on specified field
        void field
        (
            const dictionary& dict,
            const vectorField& Sf,
            const meshedSurf& surfToWrite,
            writeFile& writer,
            autoPtr<writeFile>& probesWriter
        );

        //- Calculation of deviation of specified field
        void axialDeviation
        (
            const dictionary& dict,
            const vectorField& Sf,
            const meshedSurf& surfToWrite,
            writeFile& writer,
            autoPtr<writeFile>& probesWriter
        );

        //- Calculate the swirl angle
        void swirlAngle
        (
            const dictionary& dict,
            const vectorField& Cf,
            const vectorField& Sf,
            const meshedSurf& surfToWrite,
            writeFile& writer,
            autoPtr<writeFile>& probesWriter
        );

        //- Calculate the swirl angle based on bulk angular velocity
        void bulkSwirlAngle
        (
            const dictionary& dict,
            const vectorField& Cf,
            const vectorField& Sf,
            const meshedSurf& surfToWrite,
            writeFile& writer,
            autoPtr<writeFile>& probesWriter
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldCalcs
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldCalcs/surfaceFieldCalc/surfaceFieldCalcI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldCalcs/surfaceFieldCalc/surfaceFieldCalcTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

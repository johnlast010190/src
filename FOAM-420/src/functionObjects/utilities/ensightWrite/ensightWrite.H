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
    Foam::functionObjects::ensightWrite

Group
    grpUtilitiesFunctionObjects

Description
    Writes fields in ensight format.

    Example of function object specification:
    \verbatim
    ensight
    {
        type            ensightWrite;
        libs            ("libutilityFunctionObjects.so");
        writeControl    writeTime;
        writeInterval   1;
        format          binary;

        overwrite       true;
        width           12;
        directory       "EnSight";

        fields
        (
            U
            p
        );
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description                  | Required  | Default value
        type         | Type name: ensightWrite      | yes       |
        fields       | Fields to output             | yes       |
        writeControl | Output control               | recommended | timeStep
        directory    | The output directory name    | no        | "ensightWrite"
        overwrite    | Remove existing directory    | no        | false
        format       | ASCII or binary format       | no        | same as simulation
        width        | Mask width for \c data/XXXX  | no        | 8
        noPatches    | Suppress writing patches     | no        | false
        patches      | Select patches to write      | no        |
        faceZones    | Select faceZones to write    | no        |
        consecutive  | Consecutive output numbering | no        | false
    \endtable

    Note that if the \c patches entry is an empty list, this will select all
    patches and suppress writing the internalMesh.
    Consecutive output numbering can be used in conjunction with \c overwrite.

See also
    Foam::functionObjects::vtkWrite
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    ensightWrite.C
    ensightWriteTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ensightWrite_H
#define functionObjects_ensightWrite_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "ensight/file/ensightCase.H"
#include "ensight/mesh/ensightMesh.H"

#include "primitives/strings/lists/wordReList.H"
#include "interpolation/interpolation/interpolation/interpolation.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMeshSubsetProxy/fvMeshSubsetProxy.H"
#include "primitives/enums/Enum.H"
#include "sets/topoSetSource/topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class ensightWrite Declaration
\*---------------------------------------------------------------------------*/

class ensightWrite
:
    public fvMeshFunctionObject
{
public:

    // Public static data
        static const Enum<topoSetSource::setAction> actionNames;

private:

    // Private data

        //- Mesh subset handler
        fvMeshSubset meshSubset_;

        //- Writer options
        ensightMesh::options writeOpts_;
        ensightCase::options caseOpts_;

        //- Name of fields to process
        wordReList selectFields_;

        //- Dictionary of volume selections
        dictionary selection_;

        //- Output directory name
        fileName dirName_;

        //- Consecutive output numbering
        bool consecutive_;

        //- Track changes in mesh geometry
        enum polyMesh::readUpdateState meshState_;

        //- Ensight case handler
        autoPtr<ensightCase> ensCase_;

        //- Ensight mesh handler
        autoPtr<ensightMesh> ensMesh_;


    // Private Member Functions

        //- Ensight case handler
        ensightCase& ensCase()
        {
            return ensCase_();
        }

        //- Ensight mesh handler
        ensightMesh& ensMesh()
        {
            return *ensMesh_;
        }


        //- Apply for the volume field type
        template<class Type>
        int writeVolField
        (
            const fvMeshSubset& proxy,
            const word& inputName,
            int& state
        );

        //- Create Ensight case and return if performing a static restart
        bool createCase();

        //- Update mesh subset according to zones, geometry, bounds
        bool updateSubset(fvMeshSubset& subsetter);

        bool readSelection(const dictionary& dict);

        //- Update meshes, subMeshes etc.
        bool update();

        //- Process by trying to apply for various volume field types.
        int process(const word& inputName);

        //- Disallow default bitwise copy construct
        ensightWrite(const ensightWrite&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ensightWrite&) = delete;

public:

    //- Runtime type information
    TypeName("ensightWrite");


    // Constructors

        //- Construct from runTime and dictionary.
        ensightWrite
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ensightWrite();


    // Member Functions

        //- Read the ensightWrite specification
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write fields
        virtual bool write();

        //- Execute at the final time-loop, flush case file
        virtual bool end();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightWrite/ensightWriteTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

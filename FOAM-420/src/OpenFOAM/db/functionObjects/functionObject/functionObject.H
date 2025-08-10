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
    (c) 2017-2018 OpenCFD Ltd.

Namespace
    Foam::functionObjects

Description
    Namespace for functionObjects.

    OpenFOAM includes a collection of functionObjects selected by the user at
    run-time to manipulate the simulation and provide mechanisms to extract
    field and derived quantities.  Alternatively, the same actions can be
    executed after the simulation using the \c -postProcess command-line option.

    \subsection secFunctionObjects Using function objects

    FunctionObjects are selected by additional entries in the
    $FOAM_CASE/system/controlDict dictionary.  Each object is listed in the \c
    functions sub-dictionary, e.g. to select the \c functionObjectType
    functionObject the following entry would be specified:

    \verbatim
    functions
    {
        <functionObjectName>
        {
            type                functionObjectType;
            libs                ("libMyFunctionObjectlib.so");
            region              defaultRegion;
            enabled             yes;
            timeStart           0;
            timeEnd             10;
            writeControl        writeTime;
            writeInterval       1;
            ...
        }
    }
    \endverbatim

    Where:
    \table
        Property | Description                           | Required | Default
        type     | Type of function object               | yes      |
        libs     | Libraries containing implementation   | yes      |
        region   | Name of region for multi-region cases | no       |
        enabled  | On/off switch                         | no       | yes
        log      | Log information to standard output    | no       | yes
        timeStart| Start time                            | no       |
        timeEnd  | End time                              | no       |
        executeControl  | See time controls below        | no       | timeStep
        executeInterval | Steps between each execute phase | no     |
        writeControl    | See time controls below        | no       | timeStep
        writeInterval   | Steps between each write phase | no       |
    \endtable

    Time controls:
    \table
        Option            | Description
        timeStep          | Execute/write every 'Interval' time-steps
        writeTime         | Execute/write every 'Interval' output times
        adjustableRunTime | Execute/write every 'Interval' run time period
        runTime           | Execute/write every 'Interval' run time period
        clockTime         | Execute/write every 'Interval' clock time period
        cpuTime           | Execute/write every 'Interval' CPU time period
        none              | Execute/write disabled
    \endtable

    The sub-dictionary name \c \<functionObjectName\> is chosen by the user, and
    is typically used as the name of the output directory for any data written
    by the functionObject.  The \c type entry defines the type of function
    object properties that follow.  FunctionObjects are packaged into separate
    libraries and the \c libs entry is used to specify which library should be
    loaded.

    Each function object has two separate run phases:

      - The \c execute phase is meant to be used for updating calculations
        or for management tasks.
      - The \c write phase is meant for writing the calculated data to disk.

    For each phase the respective time controls are provided, as listed above.

Class
    Foam::functionObject

Description
    Abstract base-class for Time/database function objects.

See also
    Foam::functionObjectList
    Foam::functionObjects::timeControl

SourceFiles
    functionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObject_H
#define functionObject_H

#include "db/typeInfo/typeInfo.H"
#include "memory/autoPtr/autoPtr.H"
#include "primitives/bools/Switch/Switch.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class functionObject Declaration
\*---------------------------------------------------------------------------*/

class functionObject
{
    // Private data

        //- Name
        const word name_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObject(const functionObject&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const functionObject&) = delete;


protected:

    // Protected Member Functions

        //- Return a scoped name, e.g. used to construct local field names
        word scopedName(const word& name) const;


public:

    // Forward declarations
    class unavailableFunctionObject;

    //- Runtime type information
    virtual const word& type() const = 0;

    static int debug;

    //- Global post-processing mode switch
    static bool postProcess;

    //- Directory prefix
    static word outputPrefix;

    //- Switch write log to Info
    Switch log;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            functionObject,
            dictionary,
            (const word& name, const Time& runTime, const dictionary& dict),
            (name, runTime, dict)
        );


    // Constructors

        //- Construct from components
        functionObject(const word& name);

        //- Return clone
        autoPtr<functionObject> clone() const
        {
            NotImplemented;
        }


    // Selectors

        //- Select from dictionary, based on its "type" entry
        static autoPtr<functionObject> New
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~functionObject() = default;


    // Member Functions

        //- Return the name of this functionObject
        const word& name() const;

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary& dict);

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        virtual bool execute() = 0;

        //- Execute using the specified subIndex.
        //  The base implementation is a no-op.
        //  \param subIndex an execution sub-index corresponding to a
        //      sub-cycle or something similar.
        virtual bool execute(const label subIndex);

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write() = 0;

        //- Called when Time::run() determines that the time-loop exits.
        //  By default it simply calls execute().
        virtual bool end();

        //- Called at the end of Time::adjustDeltaT() if adjustTime is true
        virtual bool adjustTimeStep();

        //- Did any file get changed during execution?
        virtual bool filesModified() const;

        //- Update for changes of mesh
        //  The base implementation is a no-op.
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        //  The base implementation is a no-op.
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

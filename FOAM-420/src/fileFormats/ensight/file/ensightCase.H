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
    Foam::ensightCase

Description
    Supports writing of ensight cases as well as providing common factory
    methods to open new files.

SourceFiles
    ensightCase.C

\*---------------------------------------------------------------------------*/

#ifndef ensightCase_H
#define ensightCase_H

#include "memory/autoPtr/autoPtr.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "db/IOstreams/IOstreams/InfoProxy.H"
#include "containers/HashTables/Map/Map.H"
#include "include/OSspecific.H"
#include "db/IOstreams/Pstreams/Pstream.H"

#include "ensight/file/ensightGeoFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class ensightCase;
class instant;
class Time;

/*---------------------------------------------------------------------------*\
                         Class ensightCase Declaration
\*---------------------------------------------------------------------------*/

class ensightCase
{
public:

    // Forward declarations
    class options;

    // Public Data

        //- The name for "data" subdirectory
        static const char* dataDirName;

        //- The name for geometry files
        static const char* geometryName;

private:

    // Private data

        //- Case writing options
        const options* options_;

        //- Output path (absolute)
        fileName ensightDir_;

        //- Case name (with ".case" ending)
        word caseName_;

        //- Output stream
        mutable OFstream* os_;

        //- Track state changes since last write
        mutable bool changed_;

        //- Time index (timeset 1)
        label timeIndex_;

        //- Time value (timeset 1)
        scalar timeValue_;

        //- Record of time index/value used (eg, field values).
        //  These values will be used for timeset 1.
        Map<scalar> timesUsed_;

        //- Record time indices when geometry is written.
        //  These values will be used to decide if timeset 1
        //  or a separate timeset are used.
        //  The special index '-1' is used static geometry.
        mutable labelHashSet geomTimes_;

        //- Record time indices when clouds are written.
        //  These values will be used to decide if timeset 1
        //  or a separate timeset are used.
        mutable labelHashSet cloudTimes_;

        //- Fields/Variables with the ensight type
        mutable HashTable<string> variables_;

        //- Cloud names and variables
        mutable HashTable<HashTable<string>> cloudVars_;


    // Private Member Functions

        //- The data directory
        fileName dataDir() const;

        //- Initial file management (master only)
        void initialize();

        //- Check if timeset uses different times than from time-set 1
        label checkTimeset(const labelHashSet& lookup) const;

        //- Write the header into the case file.
        void writeHeader() const;

        //- Write the timeset 1 into the case file.
        //  Return the time correction in effect
        scalar writeTimeset() const;

        //- Write the timeset into the case file.
        void writeTimeset
        (
            const label ts,
            const labelHashSet& lookup,
            const scalar timeCorrection = 0
        ) const;


        //- Note geometry being used
        void noteGeometry(const bool moving) const;

        //- Note cloud being used
        void noteCloud(const word& cloudName) const;

        //- Note cloud/variable being used
        void noteCloud
        (
            const word& cloudName,
            const word& varName,
            const char* ensightType
        ) const;

        //- Note field variable being used
        void noteVariable
        (
            const word& varName,
            const char* ensightType
        ) const;


        //- Open stream for new data file (on master), using the current index.
        //  File is without initial description lines.
        autoPtr<ensightFile> createDataFile(const word&) const;

        //- Open stream for new cloud file (on master).
        //  File is without initial description lines.
        autoPtr<ensightFile> createCloudFile
        (
            const word& cloudName,
            const word& name
        ) const;


        //- Disallow default bitwise copy construct
        ensightCase(const ensightCase&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ensightCase&) = delete;


public:


    // Constructors

        //- Construct from components
        ensightCase
        (
            const fileName& ensightDir,
            const word& caseName,
            const options& opts
        );

        //- Construct from components with all default options
        ensightCase
        (
            const fileName& ensightDir,
            const word& caseName,
            const IOstream::streamFormat format = IOstream::BINARY
        );



    //- Destructor
    ~ensightCase();


    // Member Functions

    // Access

        //- Reference to the case options
        inline const ensightCase::options& option() const;

        //- Ascii/Binary file output
        inline IOstream::streamFormat format() const;

        //- The nominal path to the case file
        inline const fileName& path() const;

        //- The output '*' mask
        inline const word& mask() const;

        //- Consistent zero-padded integer value
        inline word padded(const label i) const;

        //- Use values per nodes instead of per element
        inline bool nodeValues() const;

        //- Write clouds into their own directory instead in "data" directory
        inline bool separateCloud() const;


    // Edit

        //- Set time for time-set 1, using next available index.
        //  Create corresponding sub-directory.
        //  Do not mix between nextTime and setTime in an application.
        void nextTime(const scalar t);

        //- Set time for time-set 1, using next available index.
        //  Create corresponding sub-directory.
        //  Do not mix between nextTime and setTime in an application.
        void nextTime(const instant& t);

        //- Set current index and time for time-set 1.
        //  Create corresponding sub-directory
        //  Do not mix between nextTime and setTime in an application.
        void setTime(const scalar t, const label index);

        //- Set current index and time for time-set 1.
        //  Create corresponding sub-directory
        //  Do not mix between nextTime and setTime in an application.
        void setTime(const instant& t, const label index);

        // set geometry update times
        void setGeomTimes(const labelList& ts)
        {
            geomTimes_ = labelHashSet(ts);
        };

    // Addition of entries to case file

        //- Open stream for new geometry file (on master).
        autoPtr<ensightGeoFile> newGeometry(const bool moving = false) const;


        //- Open stream for new cloud positions (on master).
        //  Note the use of ensightFile, not ensightGeoFile.
        autoPtr<ensightFile> newCloud
        (
            const word& cloudName
        ) const;


        //- Open stream for new data file (on master), using the current index.
        template<class Type>
        autoPtr<ensightFile> newData(const word& varName) const;


        //- Open stream for new cloud data file (on master), using the current index.
        template<class Type>
        autoPtr<ensightFile> newCloudData
        (
            const word& cloudName,
            const word& varName
        ) const;


    // Output

        //- Rewind the output stream (master only).
        void rewind() const;

        //- Write the case file
        void write() const;

        //- Output stream (master only).
        inline Ostream& operator()() const;

        //- Print some general information.
        Ostream& printInfo(Ostream&) const;
};


//- Configuration options for the ensightCase
class ensightCase::options
{
private:

    //- Ascii/Binary file output
    IOstream::streamFormat format_;

    //- Width of mask for subdirectories
    label width_;

    //- The '*' mask appropriate for subdirectories
    word mask_;

    //- The printf format for zero-padded subdirectory numbers
    string printf_;

    //- Remove existing directory and sub-directories on creation
    bool overwrite_;

    //- Write values at nodes
    bool nodeValues_;

    //- Write clouds into their own directory
    bool separateCloud_;

public:

    // Constructors

        //- Construct with the specified format (default is binary)
        options(IOstream::streamFormat format = IOstream::BINARY);


    // Member Functions

    // Access

        //- Ascii/Binary file output
        IOstream::streamFormat format() const;

        //- The '*' mask appropriate for sub-directories
        const word& mask() const;

        //- Consistent zero-padded integer value
        word padded(const label i) const;

        //- Return current width of mask and padded.
        label width() const;

        //- Remove existing directory and sub-directories on creation
        bool overwrite() const;

        //- Use values per nodes instead of per element
        bool nodeValues() const;

        //- Write clouds into their own directory instead in "data" directory
        bool separateCloud() const;


    // Edit

        //- Set width of mask and padded.
        //  Default width is 8 digits, max width is 31 digits.
        void width(const label i);

        //- Remove existing directory and sub-directories on creation
        void overwrite(bool);

        //- Use values per nodes instead of per element
        void nodeValues(bool);

        //- Write clouds into their own directory instead in "data" directory
        void separateCloud(bool);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensight/file/ensightCaseI.H"

#ifdef NoRepository
    #include "ensight/file/ensightCaseTemplates.C"
#endif

#endif

// ************************************************************************* //

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
    (c) 2011-17 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd

Class
    Foam::OFstream

Description
    Output to file stream, using a OSstream

SourceFiles
    OFstream.C

\*---------------------------------------------------------------------------*/

#ifndef OFstream_H
#define OFstream_H

#include "db/IOstreams/Sstreams/OSstream.H"
#include "primitives/strings/fileName/fileName.H"
#include "db/typeInfo/className.H"

#include <fstream>
using std::ofstream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class OFstreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- A std::ostream with the ability to handle compressed files
class OFstreamAllocator
{
protected:

    // Member Data

        //- The allocated stream pointer (ofstream or ogzstream).
        // Actually allocated on the heap, should use smart pointers instead
        std::ostream* allocatedPtr_;

    // Constructors

        //- Construct from pathname
        OFstreamAllocator
        (
            const fileName& pathname,
            IOstream::compressionType compression=IOstream::UNCOMPRESSED,
            const bool append = false
        );

        // Explicitly delete copy constructors, because the destructor
        // de-allocates allocatedPtr_
        OFstreamAllocator& operator=( const OFstreamAllocator & ) = delete;
        OFstreamAllocator ( const OFstreamAllocator & ) = delete;

        // At some stage, we should probably do the following:
        // 1. Use std::unique_ptr instead of raw pointers for allocatedPtr_,
        // 2. Remove the custom destructor (no longer needed, because
        //    std::unique_ptr gives you RAII for free),
        // 3. Probably still explicitly delete the copy constructor for
        //    OFStreamAllocator (is there a sensible reason to allow copying of
        //    a file handle? I can't think of one...),
        // 4. Define default move constructors for OFStreamAllocator (safe
        //    because we're no longer using new/delete).
        //
        // However, this isn't terribly helpful for me at the moment, because
        // the OpenFOAM containers require copy construction, so we'd need to
        // update all of them to allow their use with move.

        // OFstreamAllocator& operator=( OFstreamAllocator && ) = default;
        // OFstreamAllocator ( OFstreamAllocator && ) = default;

    //- Destructor
    ~OFstreamAllocator();


    // Protected Member Functions

        //- Delete the stream pointer
        void deallocate();

};


/*---------------------------------------------------------------------------*\
                          Class OFstream Declaration
\*---------------------------------------------------------------------------*/

class OFstream
:
    public OFstreamAllocator,
    public OSstream
{
public:

    // Declare name of the class and its debug switch
    ClassName("OFstream");


    // Constructors

        //- Construct from pathname
        OFstream
        (
            const fileName& pathname,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED,
            const bool append = false
        );

        // Explicitly delete the copy and move constructors, because
        // OFstreamAllocator plays with pointers


    //- Destructor
    ~OFstream();


    // Member functions

      // Access

        //- Read/write access to the name of the stream
        using OSstream::name;


      // STL stream

        //- Access to underlying std::ostream
        virtual std::ostream& stdStream();

        //- Const access to underlying std::ostream
        virtual const std::ostream& stdStream() const;


      // Print

        //- Print description of IOstream to Ostream
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global predefined null output stream "/dev/null"
extern OFstream Snull;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2017 OpenCFD Ltd.

InClass
    Foam::StringStream

Description
    Input/output from string buffers.

SourceFiles
    StringStream.C

\*---------------------------------------------------------------------------*/

#ifndef StringStream_H
#define StringStream_H

#include "db/IOstreams/Sstreams/ISstream.H"
#include "db/IOstreams/Sstreams/OSstream.H"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class StringStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- Allocator for variants of a std stringstream
template<class StreamType>
class StringStreamAllocator
{
protected:

    // Member Data

        //- The allocated stream pointer.
        StreamType theStream_;


    // Constructors

        //- Construct null
        StringStreamAllocator()
            : theStream_()
        {}

        //- Construct from pointer, taking ownership
        StringStreamAllocator(StreamType* ptr)
            : theStream_(ptr)
        {}

        //- Construct from string
        StringStreamAllocator(const std::string& buffer)
            :  theStream_(buffer)
        {}


public:

    // Public Member Functions

        //- Get the string
        Foam::string str() const
        {
            return theStream_.str();
        }

        //- Set the string
        void str(const std::string& s)
        {
            theStream_.str(s);
        }
};


/*---------------------------------------------------------------------------*\
                        Class IStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Input from string buffer, using a ISstream
class IStringStream
:
    public StringStreamAllocator<std::istringstream>,
    public ISstream
{
public:

    // Constructors

        //- Construct from string
        IStringStream
        (
            const Foam::string& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            StringStreamAllocator<std::istringstream>(buffer),
            ISstream(theStream_, name, format, version)
        {}

        //- Construct from string and name
        IStringStream
        (
            const string& name,
            const string& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            StringStreamAllocator<std::istringstream>(buffer),
            ISstream
            (
                theStream_,
                name,
                format,
                version
            )
        {}

        //- Construct from char*
        IStringStream
        (
            const char* buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            StringStreamAllocator<std::istringstream>(buffer),
            ISstream(theStream_, name, format, version)
        {}


    //- Destructor
    ~IStringStream()
    {}


    // Member Functions

        //- Print description to Ostream
        void print(Ostream& os) const;


    // Member operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary:
        //  e.g. thing thisThing(IFstream("thingFileName")());
        Istream& operator()() const
        {
            return const_cast<IStringStream&>(*this);
        }

};


/*---------------------------------------------------------------------------*\
                        Class OStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Output to string buffer, using a OSstream
class OStringStream
:
    public StringStreamAllocator<std::ostringstream>,
    public OSstream
{
public:

    // Constructors

        //- Construct and set stream status
        OStringStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            StringStreamAllocator<std::ostringstream>(),
            OSstream(theStream_, "output", format, version)
        {}

        //- Construct as copy
        OStringStream(const OStringStream& oss)
        :
            StringStreamAllocator<std::ostringstream>(oss.str()),
            OSstream(theStream_, oss.name(), oss.format(), oss.version())
        {}


    //- Destructor
    ~OStringStream()
    {}


    // Member Functions

        //- Reset the output buffer and rewind the stream
        void reset()
        {
            this->str("");  // No other way to reset the end
            this->rewind();
        }

        //- Rewind the output stream
        void rewind()
        {
            // pubseekpos() instead of seekp() for symmetry with other classes
            theStream_.rdbuf()->pubseekpos(0, std::ios_base::out);
        }

        //- Print description to Ostream
        void print(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

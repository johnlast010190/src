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
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::UOPstream

Description
    Output inter-processor communications stream operating on external
    buffer.

SourceFiles
    UOPstream.C

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/Pstreams/Pstream.H"

#ifndef UOPstream_H
#define UOPstream_H

#include "db/IOstreams/Pstreams/UPstream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "db/IOstreams/Pstreams/PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class UOPstream Declaration
\*---------------------------------------------------------------------------*/

class UOPstream
:
    public UPstream,
    public Ostream
{
    // Private data

        int toProcNo_;

        DynamicList<char>& sendBuf_;

        const int tag_;

        const label comm_;

        const bool sendAtDestruct_;


    // Private Member Functions

        //- Write a T to the transfer buffer
        template<class T>
        inline void writeToBuffer(const T&);

        //- Write a char to the transfer buffer
        inline void writeToBuffer(const char&);

        //- Write data to the transfer buffer
        inline void writeToBuffer(const void* data, size_t count, size_t align);


public:

    // Constructors

        //- Construct given process index to send to and optional buffer size,
        //  write format and IO version
        UOPstream
        (
            const commsTypes commsType,
            const int toProcNo,
            DynamicList<char>& sendBuf,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool sendAtDestruct = true,
            streamFormat format=BINARY,
            versionNumber version=currentVersion
        );

        //- Construct given buffers
        UOPstream(const int toProcNo, PstreamBuffers&);


    //- Destructor
    ~UOPstream();


    // Member functions

        // Inquiry

            //- Return flags of output stream
            ios_base::fmtflags flags() const override
            {
                return ios_base::fmtflags(0);
            }


        // Write functions

            //- Write given buffer to given processor
            static bool write
            (
                const commsTypes commsType,
                const int toProcNo,
                const char* buf,
                const std::streamsize bufSize,
                const int tag = UPstream::msgType(),
                const label communicator = 0
            );

            //- Write next token to stream
            Ostream& write(const token&) override;

            //- Write character
            Ostream& write(const char) override;

            //- Write character string
            Ostream& write(const char*) override;

            //- Write word
            Ostream& write(const word&) override;

            //- Write string
            Ostream& write(const string&) override;

            //- Write std::string surrounded by quotes.
            //  Optional write without quotes.
            Ostream& writeQuoted
            (
                const std::string&,
                const bool quoted=true
            ) override;

            //- Write int32_t
            Ostream& write(const int32_t) override;

            //- Write int64_t
            Ostream& write(const int64_t) override;

            //- Write floatScalar
            Ostream& write(const floatScalar) override;

            //- Write doubleScalar
            Ostream& write(const doubleScalar) override;

            //- Write binary block
            Ostream& write(const char*, std::streamsize) override;

            //- Add indentation characters
            void indent() override
            {}


        // Stream state functions

            //- Flush stream
            void flush() override
            {}

            //- Add newline and flush stream
            void endl() override
            {}

            //- Get width of output field
            int width() const override
            {
                return 0;
            }

            //- Set width of output field (and return old width)
            int width(const int) override
            {
                 return 0;
            }

            //- Get precision of output field
            int precision() const override
            {
                 return 0;
            }

            //- Set precision of output field (and return old precision)
            int precision(const int) override
            {
                 return 0;
            }


        // Edit

            //- Set flags of stream
            ios_base::fmtflags flags(const ios_base::fmtflags) override
            {
                return ios_base::fmtflags(0);
            }


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

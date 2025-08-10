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
    (c) 2016-2017 OpenCFD Ltd.

Class
    Foam::fileFormats::STLCore

Description
    Core routines used when reading/writing STL files.

SourceFiles
    STLCore.C

\*---------------------------------------------------------------------------*/

#ifndef STLCore_H
#define STLCore_H

#include "stl/STLpoint.H"
#include "stl/STLtriangle.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class fileFormats::STLCore Declaration
\*---------------------------------------------------------------------------*/

class STLCore
{
public:

    // Public data types

        //- Enumeration for the format of data in the stream
        enum STLFormat
        {
            ASCII,  //!< ASCII
            BINARY, //!< BINARY
            BLOCKBINARY, //!< MULTIBLOCK BINARY
            UNKNOWN //!< Detect based on (input) content or (output) extension
        };


protected:

    // Protected Member Functions

        //- Detect 'stlb' extension as binary
        static bool isBinaryName
        (
            const fileName& filename,
            const STLFormat& format
        );


        //- Check contents to detect if the file is a binary STL.
        //  Return the estimated number of triangles or 0 on error.
        //  Additionally set  if a multi-block binary
        static int detectBinaryHeader
        (
            const fileName& filename,
            bool& blockBinary
        );


        //- Read STL binary file header.
        //  Return the opened file stream and estimated number of triangles.
        //  The stream is invalid and number of triangles is 0 on error.
        static autoPtr<std::istream> readBinaryHeader
        (
            const fileName& filename,
            label& nTrisEstimated
        );

        //- Write STL binary file and number of triangles to stream
        static void writeBinaryHeader(ostream& os, uint32_t nTris);


    // Constructors

        //- Construct null
        STLCore();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

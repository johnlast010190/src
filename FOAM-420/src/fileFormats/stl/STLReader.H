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
    (c) 2011 OpenFOAM Foundation
    (c) 2016-2017 OpenCFD Ltd.

Class
    Foam::fileFormats::STLReader

Description
    Internal class used by the STLsurfaceFormat and triSurface.

SourceFiles
    STLReader.C
    STLsurfaceFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef STLReader_H
#define STLReader_H

#include "stl/STLCore.H"
#include "meshes/meshShapes/labelledTri/labelledTri.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class fileFormats::STLReader Declaration
\*---------------------------------------------------------------------------*/

class STLReader
:
    public STLCore
{
    // Private Data

        bool sorted_;

        //- The points supporting the facets
        List<STLpoint> points_;

        //- The zones associated with the faces
        List<label> zoneIds_;

        //- The solid names, in the order of their first appearance
        List<word> names_;

        //- The solid count, in the order of their first appearance
        List<label> sizes_;

        //- The STL format used
        STLFormat format_;


    // Private Member Functions

        //- Read ASCII
        bool readASCII(const fileName& filename);

        //- Read BINARY
        bool readBINARY(const fileName& filename);

        //- Read Multi-Block BINARY
        bool readBlockBINARY(const fileName& filename);

        //- Read ASCII or BINARY
        bool readFile(const fileName& filename, const STLFormat& format);


        //- Disallow default bitwise copy construct
        STLReader(const STLReader&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const STLReader&) = delete;


public:

    // Constructors

        //- Read from file, filling in the information.
        //  Auto-detect ASCII/BINARY format.
        STLReader(const fileName& filename);

        //- Read from file, filling in the information.
        //  Manually selected choice of ASCII/BINARY/UNKNOWN(detect) formats.
        STLReader(const fileName& filename, const STLFormat& format);


    //- Destructor
    ~STLReader();


    // Member Functions

        //- Flush all values
        void clear();

        //- Calculate merge points mapping, return old to new pointMap.
        //  The merge tolerance based on ASCII or BINARY input format.
        //  \return number of unique points
        label mergePointsMap(labelList& pointMap) const;

        //- Calculate merge points mapping, return old to new pointMap.
        //  \return number of unique points
        label mergePointsMap(const scalar mergeTol, labelList& pointMap) const;

        //- File read was already sorted?
        inline bool sorted() const
        {
            return sorted_;
        }

        //- Return full access to the points
        inline List<STLpoint>& points()
        {
            return points_;
        }

        //- Return full access to the zoneIds
        inline List<label>& zoneIds()
        {
            return zoneIds_;
        }

        //- The list of solid names in the order of their first appearance
        inline List<word>& names()
        {
            return names_;
        }

        //- The list of solid sizes in the order of their first appearance
        inline List<label>& sizes()
        {
            return sizes_;
        }

        //- The STL format used (ASCII or BINARY)
        inline enum STLFormat stlFormat() const
        {
            return format_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

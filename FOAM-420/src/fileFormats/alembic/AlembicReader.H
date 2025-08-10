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
    (c) 2023 Esi Ltd.

Class
    Foam::fileFormats::AlembicReader

Description
    Internal class used by triSurface.

SourceFiles
    AlembicReader.C

\*---------------------------------------------------------------------------*/

#ifndef AlembicReader_H
#define AlembicReader_H

#include "meshes/meshShapes/labelledTri/labelledTri.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "include/Alembic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class fileFormats:: AlembicReader Declaration
\*---------------------------------------------------------------------------*/

class AlembicReader
{
    // Private Data

        //- Points supporting the facets
        DynamicList<point> points_;

        //- Groups associated with the faces
        HashTable<label> groupToPatch_;

        //- Tri faces
        DynamicList<labelledTri> faces_;

        //- Store xform and mesh data found in archive
        DynamicList<Tuple2<AbcG::IXform, DynamicList<AbcG::IPolyMesh>>> meshdata_;

        //- for interpolation, store non-transformed points
        bool storeBasePoints_;
        DynamicList<point> basePoints_;
        DynamicList<Abc::M44d> mats_;


    // Private Member Functions

        //- Read file
        bool readFile(const fileName& filename);

        //- Check format (Ogawa or HDF5(not compiled))
        void checkCoreType(AbcF::IFactory::CoreType& coreType);

        //- Read archive metadata
        void getMetaData(Abc::IArchive& archive);

        //- Print time series info
        void printTimeSeriesData
        (
            Abc::IArchive& archive
        );

        //- Print properties info
        void printProperties
        (
            Abc::ICompoundProperty iCProp,
            Abc::PropertyHeader header
        );

        //- Print object info
        void printObjectInfo(AbcG::IObject iObj);

        //- Recursively traverse compound properties
        void traverseCompoundProps
        (
            Abc::ICompoundProperty iProp,
            bool first = false
        );

        //- Recursively traverse data hierarchy and print info
        void traverseDataTree
        (
            AbcG::IObject iObj,
            bool first = false
        );

        //- Disallow default bitwise copy construct
        AlembicReader(const  AlembicReader&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const  AlembicReader&) = delete;


public:

    //- Runtime type information
    ClassName("AlembicReader");


    // Constructors

        //- Read from file
        AlembicReader(const fileName& filename, bool storeBase=false);


    //- Destructor
    ~AlembicReader();


    // Member Functions

        //- Flush mesh data
        void clear();


        //- Return full access to the points
        inline DynamicList<point>& points()
        {
            return points_;
        }

        //- Return full access to the zoneIds
        inline HashTable<label>& groupToPatch()
        {
            return groupToPatch_;
        }

        //- The list of solid names in the order of their first appearance
        inline DynamicList<labelledTri>& faces()
        {
            return faces_;
        }

        inline DynamicList<point>& basePoints()
        {
            return basePoints_;
        }

        inline DynamicList<Abc::M44d>& trafo()
        {
            return mats_;
        }

        //- Populate the mesh data lists given the time instance
        void setTimeData(index_t timeIdx);
        void setTimeData
        (
            Abc::chrono_t time,
            Abc::ISampleSelector::TimeIndexType idxType
                = Abc::ISampleSelector::TimeIndexType::kNearIndex
        );
        void setTimeData(Abc::ISampleSelector& iss);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

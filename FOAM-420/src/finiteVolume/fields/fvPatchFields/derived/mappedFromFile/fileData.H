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
    (c) 2009-2009 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

Class
    Foam::fileData

Description
    Class for readfing and mapping fluent profile input files to patches
    Currently only for scalars

SourceFiles
    fileData.C

\*---------------------------------------------------------------------------*/

#ifndef fileData_H
#define fileData_H

#include "db/dictionary/dictionary.H"
#include "db/IOstreams/IOstreams/Istream.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "triSurface/triSurface.H"

#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class fileData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fileData
{
    // Private data

        //- write back side
        bool mapBack_;

    // Private Member Functions

        ////- Disallow default bitwise copy construct
        //fileData(const fileData&);

        //- Disallow default bitwise assignment
        void operator=(const fileData&);

        //- Read fields
        void readFields();
        void readFieldsNAS();
        void readFieldsNTL();

        // Extract scalar from NAS/NTL file
        scalar parseScalar(const string& s);

        // Read a column of a given width from a fixed-format NTL file
        std::string readToken
        (
            const string& line,
            const size_t& width,
            size_t& index
        );

protected:

        //- fields
        static PtrList<Field<Type>> fields_;

        //- geometry triSurface
        static triSurface surface_;

        //- input file name
        static fileName filename_;

        //- input field name
        static word fieldName_;

        //- convert to meters
        static bool convertToMillimeters_;

        //- Octree
        static autoPtr<indexedOctree<treeDataPoint>> pointTreePtr_;

        //- geometry is loaded
        static bool geometryLoaded_;

public:

    // Constructors

        //- Construct null
        fileData();

        //- Construct from dictionary
        fileData(const dictionary&);

        //- Construct from copy
        fileData(const fileData&);

    //- Destructor
    ~fileData();


    // Member Functions

        //- map requested field to provided points
        Field<Type>* mapData
        (
            const pointField&,
            const vectorField&
        );


        const fileName& filename() const
        {
            return filename_;
        }

        const word& fieldName() const
        {
            return fieldName_;
        }

        const bool& mapBack() const
        {
            return mapBack_;
        }

    //  // Read
    //      void read(ISstream&);

    //  // unexpected character read error
    //  static void checkIO(const token&, const token::punctuationToken&);

    //  // go to next valid token
    //  static void nextValidToken(ISstream&);

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvPatchFields/derived/mappedFromFile/fileData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

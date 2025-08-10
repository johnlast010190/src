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
    Foam::fileFormats::FIREMeshReader

Description
    Read AVL/FIRE fpma, fpmb files.

Note
    Does not handle compressed versions (fpmaz, fpmbz) of these files.

SourceFiles
    FIREMeshReader.C

\*---------------------------------------------------------------------------*/

#ifndef FIREMeshReader_H
#define FIREMeshReader_H

#include "common/reader/meshReader.H"
#include "fire/FIRECore.H"
#include "primitives/ints/lists/labelList.H"
#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// forward declarations
class polyMesh;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                 Class fileFormats::FIREMeshReader Declaration
\*---------------------------------------------------------------------------*/

class FIREMeshReader
:
    public meshReader,
    public FIRECore
{

protected:

    // Protected Data

        labelList  owner_;
        labelList  neigh_;

        labelList  faceZoneId_;
        wordList   faceNames_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        FIREMeshReader(const FIREMeshReader&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const FIREMeshReader&) = delete;


        //- Read the mesh from the file(s)
        virtual bool readGeometry(const scalar scaleFactor = 1.0);

        //- Read points from file
        void readPoints(ISstream&, const scalar scaleFactor = 1.0);

        //- Read points from file
        void readFaces(ISstream&);

        //- Read cell connectivities from file
        void readCells(ISstream&);

        //- Read cell/face selections from file
        void readSelections(ISstream&);

        //-
        void reorganize();

        void addPatches(polyMesh&) const;


public:

    // Constructors

        //- Construct by reading file, optionally with scaling
        FIREMeshReader(const fileName&, const scalar scaleFactor = 1.0);


    //- Destructor
    virtual ~FIREMeshReader();


    // Member Functions

        //- Create and return polyMesh
        virtual autoPtr<polyMesh> mesh(const objectRegistry&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

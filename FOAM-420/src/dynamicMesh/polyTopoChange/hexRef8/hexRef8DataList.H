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
    (c) 2021 Esi Ltd.

Class
    Foam::hexRef8DataList

Description
    Various for reading/decomposing/reconstructing/distributing refinement
    data for a list of meshes.

SourceFiles
    hexRef8DataList.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef8DataList_H
#define hexRef8DataList_H

#include "primitives/ints/lists/labelIOList.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "containers/Lists/UPtrList/UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapPolyMesh;
class mapDistributePolyMesh;
class refinementHistory;
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class hexRef8DataList Declaration
\*---------------------------------------------------------------------------*/

class hexRef8DataList
{

private:

    // Private data

        //- mesh addressing (which proc has which mesh)
        labelList meshToProc_;

        labelListList meshNoInProc_;

        //- hexRef8 data
        PtrList<labelIOList> cellLevelPtrs_;

        PtrList<labelIOList> pointLevelPtrs_;

        PtrList<uniformDimensionedScalarField> level0EdgePtrs_;

        PtrList<refinementHistory> refHistoryPtrs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        hexRef8DataList(const hexRef8DataList&);

        //- Disallow default bitwise assignment
        void operator=(const hexRef8DataList&);


public:

    // Constructors

        //- Construct read. Assumes only first mesh in list
        //  having the files
        hexRef8DataList
        (
            PtrList<fvMesh>& myProcMeshes,
            const labelList& meshToProc,
            const labelListList& meshNoInProc
        );


    //- Destructor
    ~hexRef8DataList();


    // Member Functions

        //- In-place distribute
        void distribute
        (
            PtrList<mapDistributePolyMesh>& myProcDist
        );

        //- In-place reorder
        void reorder
        (
            PtrList<mapPolyMesh>& myProcMap
        );

        //- Write
        bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

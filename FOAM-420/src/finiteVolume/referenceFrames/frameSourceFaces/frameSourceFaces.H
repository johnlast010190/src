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
    (c) 2020-2023 Esi Ltd.

 Description
    Data structure that holds information regarding the faces inside an MRF

SourceFiles
    frameSourceFaces.C
\*---------------------------------------------------------------------------*/

#ifndef frameSourceFaces_H
#define frameSourceFaces_H

#include "primitives/strings/lists/wordList.H"
#include "primitives/ints/lists/labelList.H"
#include "fvMesh/fvMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class frameSourceFaces Declaration
\*---------------------------------------------------------------------------*/

class frameSourceFaces
:
    public regIOobject
{

// Private member data

    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Reference to the cells in the source cell set
    const labelList& cells_;

    //- Dictionary with the coefficients
    dictionary coeffsDict_;

    //- Number of faces in the zone
    label nZoneFaces_;

    //- This one is being updated by updateSourceFaces function
    labelHashSet includedPatchSet_;

    //- Included internal faces
    labelList internalFaces_;

    //- Included moving patch faces
    labelListList includedPatchFaces_;

    //- Default including moving patch faces
    labelListList defaultMovingPatchFaces_;

    //- Type of the face notMoving/moving/defaultMoving
    labelList faceType_;

    //- Marking cells that are in the frame
    boolList cellInFrame_;

    //- Types of the faces
    enum {notMoving, moving, defaultMoving};


protected:


    // Protected Member Functions

        //- Marking faces that are included and excluded in the frame source
        void categorizeFaces();

        //- Create included/excluded face lists
        void createFrameFaceLists();

        //- Return if the patch is excluded or not
        void initializeDataStructures();

        //- Printing information for debugging purposes
        void printDebugInformation();

public:

    //- Type name of the class
    TypeName("frameSourceFaces");

    //- Constructor from name, mesh, cells, active, dict
    frameSourceFaces
    (
        const word& name,
        const fvMesh& mesh,
        const labelList& cells,
        const bool& active,
        const dictionary& dict
    );

    //- Destructor
    virtual ~frameSourceFaces()
    {}

    //- Faces that are included in the frame motion
    void updateSourceFaces(const labelList& patches);

    //- Internal faces included in the source
    const labelListList& includedFaces() const;

    //- Faces that are excluded from the pick up (added contributions as well)
    const labelListList& excludedFaces() const;

    //- Internal faces included in the source
    const labelList& internalFaces() const;

    //- Write
    virtual bool writeData(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

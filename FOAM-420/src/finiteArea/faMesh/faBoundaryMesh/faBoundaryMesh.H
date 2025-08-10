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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Class
    faBoundaryMesh

Description
    Finite area boundary mesh

SourceFiles
    faBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef faBoundaryMesh_H
#define faBoundaryMesh_H

#include "faMesh/faPatches/faPatch/faPatchList.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterfacePtrsList.H"
#include "primitives/strings/lists/wordList.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class faMesh;

/*---------------------------------------------------------------------------*\
                      Class faBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class faBoundaryMesh
:
    public faPatchList,
    public regIOobject
{
    // private data

        //- Reference to mesh
        const faMesh& mesh_;

        //- Disallow construct as copy
        faBoundaryMesh(const faBoundaryMesh&);

        //- Disallow assignment
        void operator=(const faBoundaryMesh&);


public:

    //- Runtime type information
    TypeName("faBoundaryMesh");


    // Constructors

        //- Construct from dictionary
        faBoundaryMesh
        (
            const IOobject& io,
            const faMesh& fam
        );

        //- Construct given size
        faBoundaryMesh
        (
            const IOobject& io,
            const faMesh& fam,
            const label size
        );


    // Destructor - default


    // Member functions

        // Access

            //- Calculate the geometry for the patches
            //  (transformation tensors etc.)
            void calcGeometry();

            //- Return the mesh reference
            const faMesh& mesh() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            lduInterfacePtrsList interfaces() const;

            //- Return a list of patch types
            wordList types() const;

            //- Return a list of patch names
            wordList names() const;

            //- Find patch index given a name
            label findPatchID(const word& patchName) const;

            //- Return patch index for a given edge label
            label whichPatch(const label edgeIndex) const;

            //- Find patch indices given a name
            labelList findIndices(const keyType&, const bool useGroups) const;

            //- Check boundary definition
            bool checkDefinition(const bool report = false) const;

        // Edit

            //- Correct faBoundaryMesh after moving points
            void movePoints(const pointField&);

            //- Correct faBoundaryMesh after topology update
            void updateMesh();

            //- writeData member function required by regIOobject
            bool writeData(Ostream&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const faBoundaryMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

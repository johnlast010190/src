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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::extrude2DMesh

Description
    Given a 2D mesh insert all the topology changes to extrude. Does not work
    in parallel.

SourceFiles
    extrude2DMesh.C

\*---------------------------------------------------------------------------*/

#ifndef extrude2DMesh_H
#define extrude2DMesh_H

#include "db/typeInfo/typeInfo.H"
#include "primitives/ints/label/label.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/ints/lists/labelList.H"
#include "db/dictionary/dictionary.H"
#include "extrudeModel/extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class mapPolyMesh;
class mapDistributePolyMesh;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class extrude2DMesh Declaration
\*---------------------------------------------------------------------------*/

class extrude2DMesh
{
    // Private data

        //- Reference to 2D mesh
        polyMesh& mesh_;

        const dictionary dict_;

        //const dictionary patchDict_;

        const extrudeModel& model_;

        const word modelType_;

        const word patchType_;

        label frontPatchi_;
        label backPatchi_;

    // Private Member Functions

        //- Check the mesh is 2D
        void check2D() const;

        //- Find extrusion direction
        //void findExtrudeDirection();

        //- Disallow default bitwise copy construct
        extrude2DMesh(const extrude2DMesh&);

        //- Disallow default bitwise assignment
        void operator=(const extrude2DMesh&);

public:

    //- Runtime type information
    ClassName("extrude2DMesh");


    // Constructors
        extrude2DMesh
        (
            polyMesh&,
            const dictionary& dict,
            const extrudeModel& model
        );


    //- Destructor
    ~extrude2DMesh();


    // Member Functions

        //- Add front and back patches
        void addFrontBackPatches();

        //- Play commands into polyTopoChange to extrude mesh.
        void setRefinement(polyTopoChange&);

        //- Force recalculation of locally stored data on topological change
        void updateMesh(const mapPolyMesh&)
        {}

        //- Force recalculation of locally stored data for mesh distribution
        void distribute(const mapDistributePolyMesh&)
        {}

        label frontPatchi() const
        {
            return frontPatchi_;
        }

        label backPatchi() const
        {
            return backPatchi_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

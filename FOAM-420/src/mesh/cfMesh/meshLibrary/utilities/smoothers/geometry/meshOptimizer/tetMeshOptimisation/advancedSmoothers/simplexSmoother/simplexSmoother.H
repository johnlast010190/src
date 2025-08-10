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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    simplexSmoother

Description
    Mesh smoothing without any topological changes. The vertex is put into
    location which minimises the objective function.

SourceFiles
    simplexSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef simplexSmoother_H
#define simplexSmoother_H

#include "meshes/primitiveShapes/point/point.H"
#include "utilities/containers/DynList/DynList.H"
#include "utilities/meshes/primitives/partTet/partTet.H"
#include "meshes/boundBox/boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMeshSimplex;

/*---------------------------------------------------------------------------*\
                    Class simplexSmoother Declaration
\*---------------------------------------------------------------------------*/

//- class for volume optimizer
class simplexSmoother
{
    protected:

    // Protected data
        //- mesh points
        DynList<point, 128>& points_;

        //- list tets around the given vertex
        const DynList<partTet, 128>& tets_;

        //- label of the point
        const label pointI_;

        //- bound box
        boundBox bb_;

private:

    // Private member functions
        //- Disallow default bitwise copy construct
        simplexSmoother(const simplexSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const simplexSmoother&);

    public:

    // Constructor
        //- construct from partTetMeshSimplex
        simplexSmoother(partTetMeshSimplex& simplex);

    // Destructor
        virtual ~simplexSmoother();

    // Member functions
        //- improve the position of the centre node
        virtual void optimizeNodePosition(const scalar tol = 0.001) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

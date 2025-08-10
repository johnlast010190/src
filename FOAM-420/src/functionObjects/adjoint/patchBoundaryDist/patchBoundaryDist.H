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
    (c) 1991-2009 OpenCFD Ltd.

Class
    Foam::patchBoundaryDist

Description
    Distance calculation boundary faces from edge of patch.
    Compiles list of bounding edges for patch using syncTools and then uses
    octree to find closest point.

Note: Not using wallDist since this would result in significant errors under
      some circumstances where adjacent and near boundaries are parallel to
      current patch.

      A patch based edge-face wave would be a more cost effective solution,
      but requires significant investment in patch wave parallelisation
      architecture

SourceFiles
    patchBoundaryDist.C

\*---------------------------------------------------------------------------*/

#ifndef patchBoundaryDist_H
#define patchBoundaryDist_H

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class patchBoundaryDist Declaration
\*---------------------------------------------------------------------------*/

class patchBoundaryDist
:
    public scalarField
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- pacth index
        label index_;


    // Private Member Functions

        //- Do all calculations.
        void calculate();

        //- Disallow default bitwise copy construct
        patchBoundaryDist(const patchBoundaryDist&);

        //- Disallow default bitwise assignment
        void operator=(const patchBoundaryDist&);


public:

    // Constructors

        //- Construct from components
        patchBoundaryDist(const fvMesh& mesh, label patchIndex);


    // Destructor

        virtual ~patchBoundaryDist();


    // Member Functions

        const scalarField& y() const
        {
            return *this;
        }

        //- Correct for mesh geom/topo changes
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

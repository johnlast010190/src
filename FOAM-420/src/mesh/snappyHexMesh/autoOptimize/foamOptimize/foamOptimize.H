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
    (c) 2017 Esi Ltd.

Class
    Foam::foamOptimize

Description
    Perform optimisation of mesh usign sphericity metric

SourceFiles
    foamOptimize.C

\*---------------------------------------------------------------------------*/

#ifndef foamOptimize_H
#define foamOptimize_H

#include "fvMesh/fvMesh.H"
#include "hessianMeshOptimization/meshOptimization/hessianMeshOptimization.H"
#include "autoOptimize/autoOptimize.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class foamOptimize Declaration
\*---------------------------------------------------------------------------*/

class foamOptimize
:
    public autoOptimize
{

private:

      // Private data
        fvMesh& mesh_;

        dictionary coeffsDict_;

        autoPtr<hessianMeshOptimization> optimPtr_;

public:

    //- Runtime type information
    TypeName("foamOptimize");


        //- Update input points field
        void movePoints(pointField& newPoints);

        //- optimise mesh points
        void optimize();

        // Constructors
        foamOptimize
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~foamOptimize();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2014-2015 OpenFOAM Foundation

Class
    Foam::wallDependentModel

Description
    A class which provides on-demand creation and caching of wall distance and
    wall normal fields for use by multiple models.

SourceFiles
    wallDependentModel.C

\*---------------------------------------------------------------------------*/

#ifndef wallDependentModel_H
#define wallDependentModel_H

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class wallDependentModel Declaration
\*---------------------------------------------------------------------------*/

class wallDependentModel
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        wallDependentModel(const wallDependentModel&);

        //- Disallow default bitwise assignment
        void operator=(const wallDependentModel&);


public:

    // Constructors

        //- Construct from a mesh
        wallDependentModel(const fvMesh& mesh);


    //- Destructor
    virtual ~wallDependentModel();

    // Member Functions

        // Return the wall distance, creating and storing it if necessary
        const volScalarField& yWall() const;

        // Return the wall normal, creating and storing it if necessary
        const volVectorField& nWall() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

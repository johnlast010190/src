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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::patchSwitchFvMesh

Description
    Foam::patchSwitchFvMesh

SourceFiles
    patchSwitchFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef patchSwitchFvMesh_H
#define patchSwitchFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"
#include "patchSwitchFvMesh/patchSwitches/patchSwitch/patchSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class patchSwitchFvMesh Declaration
\*---------------------------------------------------------------------------*/

class patchSwitchFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- DynamicMesh dictionary
        dictionary dynamicMeshCoeffs_;

        //- Switches to activate patches
        PtrList<patchSwitch> pSwitchPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        patchSwitchFvMesh(const patchSwitchFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const patchSwitchFvMesh&);


public:

    //- Runtime type information
    TypeName("patchSwitchFvMesh");


    // Constructors

        //- Construct from IOobject
        patchSwitchFvMesh(const IOobject& io);


    //- Destructor
    ~patchSwitchFvMesh(){};


    // Member Functions

        //- Dummy update function which does not change the mesh
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

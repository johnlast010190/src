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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#ifndef fvSolutionRegistry_H
#define fvSolutionRegistry_H

#include "db/solutionRegistry/solutionRegistry.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvSolutionRegistry
:
    public solutionRegistry
{
private:
    //- Disallow default bitwise copy construct
    fvSolutionRegistry(const fvSolutionRegistry&);

    //- Disallow default bitwise assignment
    void operator=(const fvSolutionRegistry&);

protected:

    fvMesh& mesh_;

public:

    fvSolutionRegistry
    (
        const Time& runTime,
        fvMesh&,
        const word& regionName,
        bool defaultBehavior = false,
        bool adjointMultipoint = false
    );

    // Member Functions

    // Access

    //- access to mesh
    inline fvMesh& mesh() const
    {
        return mesh_;
    }
};

}

#endif

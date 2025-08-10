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

\*---------------------------------------------------------------------------*/

#ifndef solutionRegistry_H
#define solutionRegistry_H

#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/objectRegistry/objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionRegistry
:
    public objectRegistry
{
private:
    //- Disallow default bitwise copy construct
    solutionRegistry(const solutionRegistry&);

    //- Disallow default bitwise assignment
    void operator=(const solutionRegistry&);

protected:

    word regionName_;

    word meshName_;

    polyMesh& mesh_;

    bool default_;

public:

    solutionRegistry
    (
        const Time& runTime,
        polyMesh&,
        const word& regionName,
        bool defaultBehavior = false,
        bool adjointMultipoint = false
    );

    ~solutionRegistry() override;

    // Member Functions

    // Access

    inline const word& meshName() const
    {
        return meshName_;
    }

    inline const word& regionName() const
    {
        return regionName_;
    }

    const objectRegistry& registry() const;

};

}

#endif

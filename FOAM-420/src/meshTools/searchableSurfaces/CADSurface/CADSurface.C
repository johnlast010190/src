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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "CADSurface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(CADSurface, 0);
addToRunTimeSelectionTable(searchableSurface, CADSurface, dict);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CADSurface::CADSurface
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh
    (
        io,
        triSurface
        (
            io.globalFilePath(typeName),
            dict.lookupOrDefault<scalar>("linDeflection", 0.01),
            dict.lookupOrDefault<scalar>("angDeflection", 0.5),
            dict.lookupOrDefault<bool>("forceTriangulation", false),
            dict.lookupOrDefault<word>("nameRegionsBy", "solid"),
            wordHashSet(dict.lookupOrDefault<wordList>("ignoreRegions", wordList()))
        )
    )
{
    writeStats(Info);
    //checkGeometry(); //TODO
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CADSurface::~CADSurface()
{}


// ************************************************************************* //

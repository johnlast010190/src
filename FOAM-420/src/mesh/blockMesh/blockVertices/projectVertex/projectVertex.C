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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "blockVertices/projectVertex/projectVertex.H"
#include "global/unitConversion/unitConversion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "searchableSurfaces/searchableSurfacesQueries/searchableSurfacesQueries.H"
#include "fields/pointPatchFields/pointPatchField/pointConstraint/pointConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockVertices
{
    defineTypeNameAndDebug(projectVertex, 0);
    addToRunTimeSelectionTable(blockVertex, projectVertex, Istream);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockVertices::projectVertex::projectVertex
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    Istream& is
)
:
    pointVertex(dict, index, geometry, is),
    geometry_(geometry)
{
    wordList names(is);
    surfaces_.setSize(names.size());
    forAll(names, i)
    {
        surfaces_[i] = geometry_.findSurfaceID(names[i]);

        if (surfaces_[i] == -1)
        {
            FatalIOErrorInFunction(is)
                << "Cannot find surface " << names[i] << " in geometry"
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::blockVertices::projectVertex::operator point() const
{
    pointField start(1, pointVertex::operator point());

    pointField boundaryNear(start);
    List<pointConstraint> boundaryConstraint;


    // Note: how far do we need to search? Probably not further than
    //       span of surfaces themselves. Make sure to limit in case
    //       of e.g. searchablePlane which has infinite bb.
    boundBox bb(searchableSurfacesQueries::bounds(geometry_, surfaces_));
    bb.min() = max(bb.min(), point(-GREAT, -GREAT, -GREAT));
    bb.max() = min(bb.max(), point(GREAT, GREAT, GREAT));

    searchableSurfacesQueries::findNearest
    (
        geometry_,
        surfaces_,
        start,
        scalarField(start.size(), magSqr(bb.span())),
        boundaryNear,
        boundaryConstraint
    );

    return boundaryNear[0];
}


// ************************************************************************* //

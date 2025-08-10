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
    (c) 2014-2016 OpenFOAM Foundation
    (c) 2019-2021 OpenCFD Ltd

\*---------------------------------------------------------------------------*/

#include "blockEdges/BSplineEdge/BSplineEdge.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{
    defineTypeNameAndDebug(BSplineEdge, 0);

    addToRunTimeSelectionTable
    (
        blockEdge,
        BSplineEdge,
        Istream
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdges::BSplineEdge::BSplineEdge
(
    const pointField& points,
    const edge& fromTo,
    const pointField& internalPoints
)
:
    blockEdge(points, fromTo),
    BSpline
    (
        polyLine::concat(firstPoint(), internalPoints, lastPoint())
    )
{}


Foam::blockEdges::BSplineEdge::BSplineEdge
(
    const pointField& points,
    const label from,
    const label to,
    const pointField& internalPoints
)
:
    BSplineEdge(points, edge(from,to), internalPoints)
{}


Foam::blockEdges::BSplineEdge::BSplineEdge
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces&,
    const pointField& points,
    Istream& is
)
:
    blockEdge(dict, index, points, is),
    BSpline
    (
        polyLine::concat(firstPoint(), pointField(is), lastPoint())
    )
{
    token tok(is);
    is.putBack(tok);

    // Discard unused start/end tangents
    if (tok == token::BEGIN_LIST)
    {
        vector tangent0Ignored(is);
        vector tangent1Ignored(is);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::blockEdges::BSplineEdge::position(const scalar mu) const
{
    return BSpline::position(mu);
}


Foam::scalar Foam::blockEdges::BSplineEdge::length() const
{
    return BSpline::length();
}


// ************************************************************************* //

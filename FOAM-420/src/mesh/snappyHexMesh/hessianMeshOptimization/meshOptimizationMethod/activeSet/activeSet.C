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
    (c) 2014 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/meshOptimizationMethod/activeSet/activeSet.H"

namespace Foam
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
activeSet::activeSet
(
    const labelList& pointSet,
    const labelList& faceSet,
    const labelList& cellSet,
    const labelList& faceCellSet,
    const labelList& edgeSet
)
:
    pointSet_(pointSet),
    faceSet_(faceSet),
    cellSet_(cellSet),
    faceCellSet_(faceCellSet),
    edgeSet_(edgeSet),
    updated_(true)
{
}

Foam::activeSet::activeSet
(
    const labelList& pointSet
)
:
    pointSet_(pointSet),
    faceSet_(0),
    cellSet_(0),
    faceCellSet_(0),
    edgeSet_(0),
    updated_(true)
{
}

activeSet::activeSet()
:
    pointSet_(0),
    faceSet_(0),
    cellSet_(0),
    faceCellSet_(0),
    edgeSet_(0),
    updated_(false)
{}

activeSet::~activeSet()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void activeSet::update
(
    const labelList& pointSet,
    const labelList& faceSet,
    const labelList& cellSet,
    const labelList& faceCellSet,
    const labelList& edgeSet
)
{
    pointSet_.setSize(pointSet.size());
    pointSet_ = pointSet;

    faceSet_.setSize(faceSet.size());
    faceSet_ = faceSet;

    edgeSet_.setSize(edgeSet.size());
    edgeSet_ = edgeSet;

    cellSet_.setSize(cellSet.size());
    cellSet_ = cellSet;

    faceCellSet_.setSize(faceCellSet.size());
    faceCellSet_ = faceCellSet;

    updated_ = true;
}

void activeSet::updatePointSet(const labelList& pointSet)
{
    pointSet_.setSize(pointSet.size());
    pointSet_ = pointSet;
}
bool activeSet::updated() const
{
    return updated_;
}

} /* namespace Foam */

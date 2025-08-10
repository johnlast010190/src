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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

\*---------------------------------------------------------------------------*/

#include "utilities/anisotropicMeshing/edgeMeshGeometryModification/edgeMeshGeometryModification.H"
#include "include/demandDrivenData.H"
#include "db/dictionary/dictionary.H"
#include "edgeMesh/edgeMesh.H"

namespace Foam
{

// * * * * * * * * * * * * * * Private member functions* * * * * * * * * * * //

void edgeMeshGeometryModification::checkModification()
{
    if (meshDict_.found("anisotropicSources"))
    {
        modificationActive_ = true;

        const dictionary& anisotropicDict =
            meshDict_.subDict("anisotropicSources");

        coordinateModifierPtr_ = new coordinateModifier(anisotropicDict);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

edgeMeshGeometryModification::edgeMeshGeometryModification
(
    const edgeMesh& em,
    const dictionary& meshDict
)
:
    edgeMesh_(em),
    meshDict_(meshDict),
    coordinateModifierPtr_(nullptr),
    modificationActive_(false)
{
    checkModification();
}

edgeMeshGeometryModification::~edgeMeshGeometryModification()
{
    deleteDemandDrivenData(coordinateModifierPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool edgeMeshGeometryModification::activeModification() const
{
    return modificationActive_;
}

const edgeMesh* edgeMeshGeometryModification::modifyGeometry() const
{
    if (!modificationActive_)
    {
        WarningInFunction
            << "Modification is not active" << endl;

        return nullptr;
    }

    const pointField& pts = edgeMesh_.points();

    pointField newPts(pts.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(pts, pointI)
        newPts[pointI] = coordinateModifierPtr_->modifiedPoint(pts[pointI]);

    const edgeMesh* newEdgeMesh = new edgeMesh(newPts, edgeMesh_.edges());

    return newEdgeMesh;
}

const edgeMesh* edgeMeshGeometryModification::
revertGeometryModification() const
{
    if (!modificationActive_)
    {
        WarningInFunction
            << "Modification is not active" << endl;

        return nullptr;
    }

    const pointField& pts = edgeMesh_.points();

    pointField newPts(pts.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(pts, pointI)
    {
        newPts[pointI] =
            coordinateModifierPtr_->backwardModifiedPoint(pts[pointI]);
    }

    const edgeMesh* newEdgeMeshPtr = new edgeMesh(newPts, edgeMesh_.edges());

    return newEdgeMeshPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

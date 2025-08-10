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

Description

\*---------------------------------------------------------------------------*/

#include "utilities/boundaryLayers/triangulateNonPlanarBaseFaces/triangulateNonPlanarBaseFaces.H"
#include "db/dictionary/dictionary.H"
#include "include/demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triangulateNonPlanarBaseFaces::triangulateNonPlanarBaseFaces
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    invertedCell_(mesh_.cells().size(), false),
    decomposeFace_(mesh_.faces().size(), false),
    tol_(0.5)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

triangulateNonPlanarBaseFaces::~triangulateNonPlanarBaseFaces()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triangulateNonPlanarBaseFaces::setRelativeTolerance(const scalar tol)
{
    tol_ = tol;
}

void triangulateNonPlanarBaseFaces::triangulateLayers()
{
    if (findNonPlanarBoundaryFaces())
    {
        Info<< "Decomposing twisted boundary faces" << endl;

        decomposeBoundaryFaces();

        decomposeCellsIntoPyramids();
    }
    else
    {
        Info<< "All boundary faces are flat" << endl;
    }
}

void triangulateNonPlanarBaseFaces::readSettings
(
    const dictionary& meshDict,
    triangulateNonPlanarBaseFaces& triangulator
)
{
    if (meshDict.found("boundaryLayers"))
    {
        const dictionary& layersDict = meshDict.subDict("boundaryLayers");

        if (layersDict.found("optimisationParameters"))
        {
            const dictionary& optLayerDict =
                layersDict.subDict("optimisationParameters");

            if (optLayerDict.found("relFlatnessTol"))
            {
                const scalar relTol =
                    readScalar(optLayerDict.lookup("relFlatnessTol"));

                triangulator.setRelativeTolerance(relTol);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

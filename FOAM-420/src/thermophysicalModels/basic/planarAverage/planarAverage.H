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
    (c) 2023 Esi Ltd

SourceFiles
    planarAverage.C

\*---------------------------------------------------------------------------*/

#include "fields/Fields/primitiveFieldsFwd.H"
#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "memory/tmp/tmp.H"
#include "fields/volFields/volFieldsFwd.H"


namespace Foam
{

class faMesh;

namespace planarAverage
{

void addToAverage
(
    const vector& ownerPoint,
    const vector& nbrPoint,
    const scalar& hMin,
    const scalar& hMax,
    const label nPlanes,
    const vector& n,
    const scalar& ownerf,
    const scalar& nbrf,
    scalarField& avg,
    labelField& nFound
);

//- Average the field to a profile with nodes laid out as a function of height
template <class GeoScalarField, class GeoVectorField>
tmp<scalarField> createHeightProfile
(
    const GeoScalarField& f,
    const label nPlanes,
    const scalar& hMin,
    const scalar& hMax,
    const GeoVectorField& points,
    const labelList& ownerList,
    const labelList& neighbourList,
    const vector& n
);

//- Assign a supplied interpolation to field at given points
void interpolateHeightProfileToPoints
(
    const scalarField& profileNodes,
    const scalar& hMin,
    const scalar& hMax,
    const vectorField& points,
    const vector& n,
    scalarField& f
);

//- Average the supplied patch field over all planes with normal n
tmp<scalarField> planarAverage
(
    const faMesh& aMesh,
    const scalarField& f,
    const vector& n,
    label nPlanes
);

//- Average the supplied field over all planes with normal n
tmp<volScalarField>
planarAverage(const volScalarField& f, const vector& n, label nPlanes);

}
}

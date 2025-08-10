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
    (c) 2017 OpenCFD Ltd.

InNamespace
    Foam::fvc

Description
    Correct flux-U difference in the internal loop  using relaxation factor

SourceFiles
    fvcCorrectAlpha.C

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcCorrectAlpha.H"
#include "fvMesh/fvMesh.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<GeometricField<scalar, fvsPatchField, surfaceMesh>> alphaCorr
(
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& phiU,
    const bool finalIter
)
{
    const fvMesh& mesh = U.mesh();
    const word fieldName = U.select(finalIter);

    scalar alpha = 1;
    if (mesh.solution().relaxEquation(fieldName))
    {
        alpha = mesh.solution().equationRelaxationFactor(fieldName);
    }

    return
        (1 - alpha)
       *(phiU.prevIter() - (fvc::interpolate(U.prevIter()) & mesh.Sf()));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

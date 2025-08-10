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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Description

\*---------------------------------------------------------------------------*/

#include "interpolation/edgeInterpolation/edgeInterpolationScheme/edgeInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<scalar>, Mesh);
defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<vector>, Mesh);
defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<tensor>, Mesh);
defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<symmTensor>, Mesh);
defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<sphericalTensor>, Mesh);

defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<scalar>, MeshFlux);
defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<vector>, MeshFlux);
defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<tensor>, MeshFlux);
defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<symmTensor>, MeshFlux);
defineTemplateRunTimeSelectionTable(edgeInterpolationScheme<sphericalTensor>, MeshFlux);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

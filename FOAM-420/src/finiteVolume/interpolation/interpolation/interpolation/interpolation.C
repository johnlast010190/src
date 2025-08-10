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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "interpolation/interpolation/interpolation/interpolation.H"
#include "fields/volFields/volFields.H"
#include "meshes/polyMesh/polyMesh.H"
#include "fields/pointPatchFields/basic/calculated/calculatedPointPatchFields.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolation<Type>::interpolation
(
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    psi_(psi),
    pMesh_(psi.mesh()),
    pMeshPoints_(pMesh_.points()),
    pMeshFaces_(pMesh_.faces()),
    pMeshFaceCentres_(pMesh_.faceCentres()),
    pMeshFaceAreas_(pMesh_.faceAreas())
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolation/interpolation/interpolation/interpolationNew.C"

// ************************************************************************* //

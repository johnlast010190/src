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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcCurl.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
curl
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    word nameCurlVf = "curl(" + vf.name() + ')';

    // Gausses theorem curl
    // tmp<GeometricField<Type, fvPatchField, volMesh>> tcurlVf =
    //     fvc::surfaceIntegrate(vf.mesh().Sf() ^ fvc::interpolate(vf));

    // Calculate curl as the Hodge dual of the skew-symmetric part of grad
    tmp<GeometricField<Type, fvPatchField, volMesh>> tcurlVf =
        2.0*(*skew(fvc::grad(vf, nameCurlVf)));

    tcurlVf.ref().rename(nameCurlVf);

    return tcurlVf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
curl
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> Curl(fvc::curl(tvf()));
    tvf.clear();
    return Curl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

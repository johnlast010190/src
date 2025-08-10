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
    (c) 2023 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcSnGrad.H"
#include "fvMesh/fvMesh.H"
#include "finiteVolume/snGradSchemes/snGradScheme/snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
snGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name,
    const word& forceGradSchemeName
)
{
    return fv::snGradScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().snGradScheme(name),
        forceGradSchemeName
    )().snGrad(vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
snGrad
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const word& name,
    const word& forceGradSchemeName
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> SnGrad
    (
        fvc::snGrad(tvf(), name, forceGradSchemeName)
    );
    tvf.clear();
    return SnGrad;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
snGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::snGrad(vf, "snGrad(" + vf.name() + ')');
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
snGrad
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> SnGrad
    (
        fvc::snGrad(tvf())
    );
    tvf.clear();
    return SnGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

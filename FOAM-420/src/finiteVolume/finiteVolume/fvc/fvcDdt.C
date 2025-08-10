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
    (c) 2022 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcDdt.H"
#include "fvMesh/fvMesh.H"
#include "finiteVolume/ddtSchemes/ddtScheme/ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const dimensioned<Type> dt,
    const fvMesh& mesh
)
{
    return fv::ddtScheme<Type>::New
    (
        mesh,
        mesh.schemes().ddtScheme("ddt(" + dt.name() + ')')
    ).ref().fvcDdt(dt);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::getDdtScheme(vf).ref().fvcDdt(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme(name)
    ).ref().fvcDdt(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::getDdtScheme(rho, vf).ref().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme(name)
    ).ref().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::getDdtScheme(rho, vf).ref().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme(name)
    ).ref().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::getDdtScheme(alpha, rho, vf).ref().fvcDdt(alpha, rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme(name)
    ).ref().fvcDdt(alpha, rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
ddt
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    return fv::ddtScheme<Type>::New
    (
        sf.mesh(),
        sf.mesh().schemes().ddtScheme("ddt(" + sf.name() + ')')
    ).ref().fvcDdt(sf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const one&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return ddt(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ddt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const one&
)
{
    return ddt(vf);
}


template<class Type>
tmp<GeometricField<typename Foam::flux<Type>::type, fvsPatchField, surfaceMesh>>
ddtCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf,
    const word interpolationName
)
{
    return fv::getDdtScheme(U).ref().fvcDdtUfCorr(U, Uf, interpolationName);
}


template<class Type>
tmp<GeometricField<typename Foam::flux<Type>::type, fvsPatchField, surfaceMesh>>
ddtCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField
    <
        typename Foam::flux<Type>::type,
        fvsPatchField,
        surfaceMesh
    >& phi,
    const word interpolationName
)
{
    return fv::getDdtScheme(U).ref().fvcDdtPhiCorr(U, phi, interpolationName);
}


template<class Type>
tmp<GeometricField<typename Foam::flux<Type>::type, fvsPatchField, surfaceMesh>>
ddtCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf,
    const word interpolationName
)
{
    return fv::getDdtScheme(rho, U).ref().fvcDdtUfCorr
    (
        rho, U, Uf, interpolationName
    );
}


template<class Type>
tmp<GeometricField<typename Foam::flux<Type>::type, fvsPatchField, surfaceMesh>>
ddtCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField
    <
        typename Foam::flux<Type>::type,
        fvsPatchField,
        surfaceMesh
    >& phi,
    const word interpolationName
)
{
    return fv::getDdtScheme(rho, U).ref().fvcDdtPhiCorr
    (
        rho, U, phi, interpolationName
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

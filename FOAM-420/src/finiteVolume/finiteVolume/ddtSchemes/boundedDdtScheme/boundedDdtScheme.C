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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "finiteVolume/ddtSchemes/boundedDdtScheme/boundedDdtScheme.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcDdt.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
boundedDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    return scheme_.ref().fvcDdt(dt);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
boundedDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_.ref().fvcDdt(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
boundedDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_.ref().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
boundedDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_.ref().fvcDdt(rho, vf) - fvc::ddt(rho)*vf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
boundedDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_.ref().fvcDdt(alpha, rho, vf) - fvc::ddt(alpha, rho)*vf;
}


template<class Type>
tmp<fvMatrix<Type>>
boundedDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_.ref().fvmDdt(vf);
}


template<class Type>
tmp<fvMatrix<Type>>
boundedDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_.ref().fvmDdt(rho, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
boundedDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_.ref().fvmDdt(rho, vf) - fvm::Sp(fvc::ddt(rho), vf);
}


template<class Type>
tmp<fvMatrix<Type>>
boundedDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return
        scheme_.ref().fvmDdt(alpha, rho, vf)
      - fvm::Sp(fvc::ddt(alpha, rho), vf);
}


template<class Type>
tmp<typename boundedDdtScheme<Type>::fluxFieldType>
boundedDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf,
    const word interpolationName
)
{
    return scheme_.ref().fvcDdtUfCorr(U, Uf, interpolationName);
}


template<class Type>
tmp<typename boundedDdtScheme<Type>::fluxFieldType>
boundedDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const word interpolationName
)
{
    return scheme_.ref().fvcDdtPhiCorr(U, phi, interpolationName);
}


template<class Type>
tmp<typename boundedDdtScheme<Type>::fluxFieldType>
boundedDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf,
    const word interpolationName
)
{
    return scheme_.ref().fvcDdtUfCorr(rho, U, Uf, interpolationName);
}


template<class Type>
tmp<typename boundedDdtScheme<Type>::fluxFieldType>
boundedDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const word interpolationName
)
{
    return scheme_.ref().fvcDdtPhiCorr(rho, U, phi, interpolationName);
}


template<class Type>
tmp<surfaceScalarField> boundedDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_.ref().meshPhi(vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

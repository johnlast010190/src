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
    (c) 2020-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/fvc/fvcDDt.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "fvMesh/fvMesh.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
DDt
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> ddtDivPhiPsi
        = fvc::ddt(psi) + fvc::div(phi, psi);

    if (phi.mesh().moving())
    {
        return ddtDivPhiPsi - fvc::div(fvc::absolute(phi, psi))*psi;
    }
    else
    {
        return ddtDivPhiPsi - fvc::div(phi)*psi;
    }
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
DDt
(
    const geometricOneField& rho,
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    return fvc::DDt(phi, psi);
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
DDt
(
    const volScalarField& rho,
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    if (phi.dimensions() == dimMass/dimTime)
    {
        return
            1/rho
           *(
               fvc::ddt(rho, psi) + fvc::div(phi, psi)
             - (fvc::ddt(rho) + fvc::div(phi))*psi
            );
    }
    else
    {
        return fvc::DDt(phi, psi);
    }
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
DDt
(
    const tmp<surfaceScalarField>& tphi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> DDtPsi
    (
        fvc::DDt(tphi(), psi)
    );
    tphi.clear();
    return DDtPsi;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> DDt
(
    const geometricOneField& rho,
    const tmp<surfaceScalarField>& tphi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> DDtPsi
    (
        fvc::DDt(rho, tphi(), psi)
    );
    tphi.clear();
    return DDtPsi;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> DDt
(
    const volScalarField& rho,
    const tmp<surfaceScalarField>& tphi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> DDtPsi
    (
        fvc::DDt(rho, tphi(), psi)
    );
    tphi.clear();
    return DDtPsi;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> DDt
(
    const tmp<geometricOneField>& trho,
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> DDtPsi
    (
        fvc::DDt(trho(), phi, psi)
    );
    trho.clear();
    return DDtPsi;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> DDt
(
    const tmp<volScalarField>& trho,
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> DDtPsi
    (
        fvc::DDt(trho(), phi, psi)
    );
    trho.clear();
    return DDtPsi;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> DDt
(
    const tmp<geometricOneField>& trho,
    const tmp<surfaceScalarField>& tphi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> DDtPsi
    (
        fvc::DDt(trho(), tphi(), psi)
    );
    tphi.clear();
    trho.clear();
    return DDtPsi;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> DDt
(
    const tmp<volScalarField>& trho,
    const tmp<surfaceScalarField>& tphi,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> DDtPsi
    (
        fvc::DDt(trho(), tphi(), psi)
    );
    tphi.clear();
    trho.clear();
    return DDtPsi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

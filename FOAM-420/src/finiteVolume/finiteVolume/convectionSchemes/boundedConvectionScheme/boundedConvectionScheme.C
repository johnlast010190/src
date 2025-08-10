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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/convectionSchemes/boundedConvectionScheme/boundedConvectionScheme.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const convectionScheme<Type>& boundedConvectionScheme<Type>::scheme() const
{
    return scheme_();
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
boundedConvectionScheme<Type>::interpolate
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return scheme_().interpolate(phi, vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
boundedConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return scheme_().flux(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
boundedConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> absFlux = absoluteFlux(faceFlux, vf);
    return
        scheme_().fvmDiv(faceFlux, vf)
      + fvm::SuSp(-fvc::surfaceIntegrate(absFlux), vf);
}


template<class Type>
tmp<fvMatrix<Type>>
boundedConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const volScalarField& vfMult,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> absFlux = absoluteFlux(faceFlux, vf);
    return
        scheme_().fvmDiv(faceFlux, vfMult, vf)
      + fvm::SuSp(-fvc::surfaceIntegrate(absFlux)*vfMult, vf);
}


template<class Type>
tmp<fvBlockMatrix<Type>>
boundedConvectionScheme<Type>::fvmBDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<fvBlockMatrix<Type>> tbfvm
    (
        new fvBlockMatrix<Type>
        (
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>(vf)
        )
    );

    tmp<fvMatrix<Type>> segConv = fvmDiv(faceFlux, vf);

    tbfvm->insertEquation(0, segConv.ref());

    return tbfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
boundedConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> absFlux = absoluteFlux(faceFlux, vf);
    return
        scheme_().fvcDiv(faceFlux, vf)
      - fvc::surfaceIntegrate(absFlux)*vf;
}


template<class Type>
tmp<surfaceScalarField> boundedConvectionScheme<Type>::absoluteFlux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (faceFlux.mesh().moving() && faceFlux.dimensions() == dimVolume/dimTime)
    {
        return fvc::absolute(faceFlux, vf);
    }
    else
    {
        return tmp<surfaceScalarField>(faceFlux);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

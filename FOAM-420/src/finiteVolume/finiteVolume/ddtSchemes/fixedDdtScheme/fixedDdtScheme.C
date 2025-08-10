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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/ddtSchemes/fixedDdtScheme/fixedDdtScheme.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcDdt.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "primitives/functions/Function1/Table/TableBase.H"
#include "primitives/functions/Function1/Table/Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fixedDdtScheme<Type>::fixedDdtScheme
(
    const fvMesh& mesh, Istream& is
)
:
    ddtScheme<Type>(mesh, is),
    scheme_
    (
        fv::ddtScheme<Type>::New(mesh, is)
    )
{
    token firstToken(is);

    if (firstToken.isNumber())
    {
        const scalar scaleC = firstToken.number();
        scale_ = new Function1Types::Constant<scalar>
        (
            "scale",
            scaleC
        );
    }
    else if (firstToken.isWord())
    {
        const word functionType = firstToken.wordToken();
        if (functionType == "constant")
        {
            dictionary functionDict=dictionary();
            scalar scale(readScalar(is));
            scale_ = new Function1Types::Constant<scalar>
            (
                "scale",
                scale
            );
        }
        else if (functionType == "table")
        {
            dictionary functionDict=dictionary();
            List<Tuple2<scalar, scalar>> tup(is);
            scale_ = new Function1Types::Table<scalar>
            (
                "scale",
                tup
            );
        }
        else
        {
            FatalError << "functionType = " <<  functionType
                       << " is not supported."
                       << "Only constant and table types are supported."
                       << endl;
        }
    }
    else
    {
        FatalError << "Specify either a scalar or a Function1: constant/table"
                   << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
fixedDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    return scale()*scheme_.ref().fvcDdt(dt);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
fixedDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scale()*scheme_.ref().fvcDdt(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
fixedDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scale()*scheme_.ref().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
fixedDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scale()*scheme_.ref().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
fixedDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scale()*scheme_.ref().fvcDdt(alpha, rho, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
fixedDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scale()*scheme_.ref().fvmDdt(vf);
}


template<class Type>
tmp<fvMatrix<Type>>
fixedDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scale()*scheme_.ref().fvmDdt(rho, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
fixedDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scale()*scheme_.ref().fvmDdt(rho, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
fixedDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scale()*scheme_.ref().fvmDdt(alpha, rho, vf);
}


template<class Type>
tmp<typename fixedDdtScheme<Type>::fluxFieldType>
fixedDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf,
    const word interpolationName
)
{
    return scale()*scheme_.ref().fvcDdtUfCorr(U, Uf, interpolationName);
}


template<class Type>
tmp<typename fixedDdtScheme<Type>::fluxFieldType>
fixedDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const word interpolationName
)
{
    return scale()*scheme_.ref().fvcDdtPhiCorr(U, phi, interpolationName);
}


template<class Type>
tmp<typename fixedDdtScheme<Type>::fluxFieldType>
fixedDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf,
    const word interpolationName
)
{
    return scale()*scheme_.ref().fvcDdtUfCorr(rho, U, Uf, interpolationName);
}


template<class Type>
tmp<typename fixedDdtScheme<Type>::fluxFieldType>
fixedDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const word interpolationName
)
{
    return scale()*scheme_.ref().fvcDdtPhiCorr(rho, U, phi, interpolationName);
}


template<class Type>
tmp<surfaceScalarField> fixedDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scale()*scheme_.ref().meshPhi(vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

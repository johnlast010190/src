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

#include "finiteVolume/ddtSchemes/smartBoundedDdtScheme/smartBoundedDdtScheme.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcDdt.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "cfdTools/general/fvOptions/fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<volScalarField>
smartBoundedDdtScheme<Type>::fvOptSrc(const volScalarField& rho) const
{
    if (!fvOptionsPtr_)
    {
        // Looks up or creates the fvOptions
        fvOptionsPtr_ = &fv::options::New(rho.mesh(), rho.db());
    }
    return fvOptionsPtr_->operator()(const_cast<volScalarField&>(rho))&rho;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
smartBoundedDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (thermo(vf.db()).isochoric() && !rho.mesh().moving())
    {
        return rho*this->scheme().ref().fvcDdt(vf) + fvOptSrc(rho)*vf;
    }
    else
    {
        tmp<fv::ddtScheme<scalar>> rhoScheme
        (
            fv::ddtScheme<scalar>::New
            (
                vf.mesh(), vf.mesh().schemes().ddtScheme
                (
                    "ddt(" + rho.name() + ',' + vf.name() + ')'
                )
            )
        );
        return this->scheme().ref().fvcDdt(rho, vf)
            - rhoScheme.ref().fvcDdt(rho)*vf + fvOptSrc(rho)*vf;
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
smartBoundedDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (thermo(vf.db()).isochoric() && !rho.mesh().moving())
    {
        return
            this->scheme().ref().fvcDdt(alpha, rho, vf)
          + alpha*fvOptSrc(rho)*vf;
    }
    else
    {
        tmp<fv::ddtScheme<scalar>> rhoScheme
        (
            fv::ddtScheme<scalar>::New
            (
                vf.mesh(), vf.mesh().schemes().ddtScheme
                (
                    "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')'
                )
            )
        );
        return
            this->scheme().ref().fvcDdt(alpha, rho, vf)
          - rhoScheme.ref().fvcDdt(alpha, rho)*vf
          + alpha*fvOptSrc(rho)*vf;
    }
}


template<class Type>
tmp<fvMatrix<Type>>
smartBoundedDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (thermo(vf.db()).isochoric() && !rho.mesh().moving())
    {
        return
            rho*this->scheme().ref().fvmDdt(vf)
          + fvm::SuSp(fvOptSrc(rho), vf);
    }
    else
    {
        tmp<fv::ddtScheme<scalar>> rhoScheme
        (
            fv::ddtScheme<scalar>::New
            (
                vf.mesh(), vf.mesh().schemes().ddtScheme
                (
                    "ddt(" + rho.name() + ',' + vf.name() + ')'
                )
            )
        );
        return
            this->scheme().ref().fvmDdt(rho, vf)
          + fvm::SuSp
            (
                -rhoScheme.ref().fvcDdt(rho) + fvOptSrc(rho), vf
            );
    }
}


template<class Type>
tmp<fvMatrix<Type>>
smartBoundedDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (thermo(vf.db()).isochoric() && !rho.mesh().moving())
    {
        return
            this->scheme().ref().fvmDdt(alpha, rho, vf)
          + fvm::SuSp(fvOptSrc(rho), vf);
    }
    else
    {
        tmp<fv::ddtScheme<scalar>> rhoScheme
        (
            fv::ddtScheme<scalar>::New
            (
                vf.mesh(), vf.mesh().schemes().ddtScheme
                (
                    "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')'
                )
            )
        );
        return
            this->scheme().ref().fvmDdt(alpha, rho, vf)
          + fvm::SuSp(-rhoScheme.ref().fvcDdt(alpha, rho) + fvOptSrc(rho), vf);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

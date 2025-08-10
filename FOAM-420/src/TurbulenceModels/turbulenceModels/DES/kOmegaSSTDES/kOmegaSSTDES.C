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
    (c) 2016 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "kOmegaSSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTDES<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2,
    bool createFvOptions
)
{
    // Correct the turbulence viscosity
    kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::correctNut(S2, createFvOptions);

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTDES<BasicTurbulenceModel>::correctNut()
{
    correctNut(true);
}


template<class BasicTurbulenceModel>
void kOmegaSSTDES<BasicTurbulenceModel>::correctNut(bool createFvOptions)
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), createFvOptions);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& magGradU,
    const volScalarField& CDES
) const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;

    return min(CDES*this->delta(), sqrt(k)/(this->betaStar_*omega));
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTDES<BasicTurbulenceModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU,
    const bool deferred
) const
{
    volScalarField CDES(this->CDES(F1));
    return sqrt(this->k_())/dTilda(mag(gradU), CDES)()();
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTDES<BasicTurbulenceModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return GbyNu0; // Unlimited
}

template<class BasicTurbulenceModel>
void kOmegaSSTDES<BasicTurbulenceModel>::correct()
{
    kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::correct();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTDES<BasicTurbulenceModel>::kOmegaSSTDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSSTBase<DESModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    CDESkom_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDESkom",
            this->coeffDict_,
            0.82
        )
    ),
    CDESkeps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDESkeps",
            this->coeffDict_,
            0.60
        )
    )
{
    correctNut(false);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::read())
    {
        kappa_.readIfPresent(this->coeffDict());
        CDESkom_.readIfPresent(this->coeffDict());
        CDESkeps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicTurbulenceModel>::LESRegion() const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const volVectorField& U = this->U_;

    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k) & fvc::grad(omega))/omega
    );

    const volScalarField F1(this->F1(CDkOmega));

    tmp<volScalarField> tLESRegion
    (
        new volScalarField
        (
            IOobject
            (
                "DES::LESRegion",
                this->mesh_.time().timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neg
            (
                dTilda
                (
                    mag(fvc::grad(U)),
                    F1*CDESkom_ + (1 - F1)*CDESkeps_
                )
              - sqrt(k)/(this->betaStar_*omega)
            )
        )
    );

    return tLESRegion;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //

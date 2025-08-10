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
    (c) 2015 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "EddyDiffusivity/EddyDiffusivity.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::EddyDiffusivity<BasicTurbulenceModel>::correctNut()
{
    alphat_ = this->rho_*this->nut()/Prt_;
    alphat_.correctBoundaryConditions();

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::EddyDiffusivity<BasicTurbulenceModel>::EddyDiffusivity
(
    const word& type,
    const alphaField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
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

    // coeffsDict not available yet
    Prt_("Prt", dimless, 0.85),
    Sct_("Sct", dimless, 0.7),

    alphat_
    (
        IOobject
        (
            IOobject::groupName("alphat", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Foam::EddyDiffusivity<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        Prt_.readIfPresent(this->coeffDict());
        Sct_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::EddyDiffusivity<BasicTurbulenceModel>::Dmt() const
{
    tmp<volScalarField> Dmt = this->rho_*this->nut()/Sct_;
    Dmt->correctBoundaryConditions();

    return Dmt;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::EddyDiffusivity<BasicTurbulenceModel>::Dmt(const label patchi) const
{
    tmp<scalarField> Dmt =
        this->rho_.boundaryField()[patchi]*this->nut(patchi)/Sct_.value();
    return Dmt;
}



template<class BasicTurbulenceModel>
void Foam::EddyDiffusivity<BasicTurbulenceModel>::correctEnergyTransport()
{
    EddyDiffusivity<BasicTurbulenceModel>::correctNut();
}


// ************************************************************************* //

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
    (c) 2010-2016 Esi Ltd.
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "IncompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::IncompressibleTurbulenceModel<TransportModel>::
IncompressibleTurbulenceModel
(
    const word& type,
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transport,
    const word& propertiesName
)
:
    TurbulenceModel
    <
        geometricOneField,
        geometricOneField,
        incompressibleTurbulenceModel,
        TransportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::IncompressibleTurbulenceModel<TransportModel>>
Foam::IncompressibleTurbulenceModel<TransportModel>::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const TransportModel& transport,
    const word& propertiesName
)
{
    autoPtr<IncompressibleTurbulenceModel> modelPtr
    (
        static_cast<IncompressibleTurbulenceModel*>(
        TurbulenceModel
        <
            geometricOneField,
            geometricOneField,
            incompressibleTurbulenceModel,
            TransportModel
        >::New
        (
            geometricOneField(),
            geometricOneField(),
            U,
            phi,
            phi,
            transport,
            propertiesName
        ).ptr())
    );
    modelPtr->read();
    return modelPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleTurbulenceModel<TransportModel>::devReff() const
{
    return devRhoReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U
) const
{
    return divDevRhoReff(U);
}

template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U, const word& scheme
) const
{
    return divDevRhoReff(U, scheme);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleTurbulenceModel<TransportModel>::
devRhoReff() const
{
    NotImplemented;
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleTurbulenceModel<TransportModel>::
divDevRhoReff
(
    volVectorField& U
) const
{
    NotImplemented;
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleTurbulenceModel<TransportModel>::
divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;
}

template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleTurbulenceModel<TransportModel>::
divDevRhoReff
(
    volVectorField& U, const word& scheme
) const
{
    NotImplemented;
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleTurbulenceModel<TransportModel>::
divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U,
    const word& scheme
) const
{
    NotImplemented;
}


template<class TransportModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::IncompressibleTurbulenceModel<TransportModel>::bDivDevReff
(
    volVectorField& U
) const
{
    return bDivDevRhoReff(U);
}

template<class TransportModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::IncompressibleTurbulenceModel<TransportModel>::
bDivDevRhoReff
(
    volVectorField& U
) const
{
    NotImplemented;
}


template<class TransportModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::IncompressibleTurbulenceModel<TransportModel>::bDivDevReff
(
    volVectorField& U,
    const word& scheme
) const
{
    return bDivDevRhoReff(U, scheme);
}

template<class TransportModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::IncompressibleTurbulenceModel<TransportModel>::
bDivDevRhoReff
(
    volVectorField& U,
    const word& scheme
) const
{
    NotImplemented;
}


// ************************************************************************* //

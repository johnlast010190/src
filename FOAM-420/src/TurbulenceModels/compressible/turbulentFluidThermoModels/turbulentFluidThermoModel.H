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
    (c) 2015-2016 OpenFOAM Foundation

Typedef
    Foam::compressible::turbulenceModel

Typedef
    Foam::compressible::RASModel

Typedef
    Foam::compressible::LESModel

Description
    Typedefs for turbulence, RAS and LES models for compressible flow
    based on the standard laminar transport package.

SourceFiles
    turbulentFluidThermoModel.C
    turbulentFluidThermoModels.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentFluidThermoModel_H
#define turbulentFluidThermoModel_H

#include "CompressibleTurbulenceModel/CompressibleTurbulenceModel.H"
#include "ThermalDiffusivity/ThermalDiffusivity.H"
#include "EddyDiffusivity/EddyDiffusivity.H"
#include "laminar/laminarModel/laminarModel.H"
#include "RAS/RASModel/RASModel.H"
#include "LES/LESModel/LESModel.H"
#include "fluidThermo/fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace compressible
    {
        typedef ThermalDiffusivity<CompressibleTurbulenceModel<fluidThermo>>
            turbulenceModel;

        typedef laminarModel<turbulenceModel> laminarModel;
        typedef RASModel<EddyDiffusivity<turbulenceModel>> RASModel;
        typedef LESModel<EddyDiffusivity<turbulenceModel>> LESModel;

        template<class BasicCompressibleTurbulenceModel>
        autoPtr<BasicCompressibleTurbulenceModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const typename BasicCompressibleTurbulenceModel::transportModel&
                transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentFluidThermoModels/turbulentFluidThermoModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

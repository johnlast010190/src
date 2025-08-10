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
    (c) 2013-2017 OpenFOAM Foundation

Class
    Foam::CompressibleTurbulenceModel

Description
    Templated abstract base class for single-phase compressible
    turbulence models.

SourceFiles
    CompressibleTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef CompressibleTurbulenceModel_H
#define CompressibleTurbulenceModel_H

#include "TurbulenceModel/TurbulenceModel.H"
#include "compressibleTurbulenceModel.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class CompressibleTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

template<class TransportModel>
class CompressibleTurbulenceModel
:
    public TurbulenceModel
    <
        geometricOneField,
        volScalarField,
        compressibleTurbulenceModel,
        TransportModel
    >
{

public:

    typedef geometricOneField alphaField;
    typedef volScalarField rhoField;
    typedef TransportModel transportModel;


    // Constructors

        //- Construct
        CompressibleTurbulenceModel
        (
            const word& type,
            const geometricOneField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<CompressibleTurbulenceModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~CompressibleTurbulenceModel()
    {}


    // Member Functions

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const
        {
            return this->transport_.mu();
        }

        //- Return the laminar dynamic viscosity on patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
            return this->transport_.mu(patchi);
        }

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return this->transport_.mu()/this->rho_;
        }

        //- Return the laminar viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return
                this->transport_.mu(patchi)/this->rho_.boundaryField()[patchi];
        }

        //- Return the turbulence dynamic viscosity
        virtual tmp<volScalarField> mut() const
        {
            return this->rho_*this->nut();
        }

        //- Return the turbulence dynamic viscosity on patch
        virtual tmp<scalarField> mut(const label patchi) const
        {
            return this->rho_.boundaryField()[patchi]*this->nut(patchi);
        }

        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const
        {
            return mut() + mu();
        }

        //- Return the effective dynamic viscosity on patch
        virtual tmp<scalarField> muEff(const label patchi) const
        {
            return mut(patchi) + mu(patchi);
        }

        //- Return the source term for the momentum equation
        virtual tmp<fvBlockMatrix<vector>> bDivDevReff
        (
            volVectorField& U
        ) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvBlockMatrix<vector>> bDivDevRhoReff
        (
            volVectorField& U
        ) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvBlockMatrix<vector>> bDivDevReff
        (
            volVectorField& U,
            const word& scheme
        ) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvBlockMatrix<vector>> bDivDevRhoReff
        (
            volVectorField& U,
            const word& scheme
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompressibleTurbulenceModel/CompressibleTurbulenceModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

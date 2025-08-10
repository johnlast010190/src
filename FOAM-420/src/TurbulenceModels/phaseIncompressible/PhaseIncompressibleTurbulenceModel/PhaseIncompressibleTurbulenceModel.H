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
    (c) 2016 Esi Ltd.
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::PhaseIncompressibleTurbulenceModel

Description
    Templated abstract base class for multiphase incompressible
    turbulence models.

SourceFiles
    PhaseIncompressibleTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseIncompressibleTurbulenceModel_H
#define PhaseIncompressibleTurbulenceModel_H

#include "TurbulenceModel/TurbulenceModel.H"
#include "incompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class PhaseIncompressibleTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

template<class TransportModel>
class PhaseIncompressibleTurbulenceModel
:
    public TurbulenceModel
    <
        volScalarField,
        geometricOneField,
        incompressibleTurbulenceModel,
        TransportModel
    >
{

public:

    typedef volScalarField alphaField;
    typedef geometricOneField rhoField;
    typedef TransportModel transportModel;


    // Constructors

        //- Construct
        PhaseIncompressibleTurbulenceModel
        (
            const word& type,
            const alphaField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const TransportModel& transportModel,
            const word& propertiesName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<PhaseIncompressibleTurbulenceModel> New
        (
            const alphaField& alpha,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const TransportModel& transportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~PhaseIncompressibleTurbulenceModel()
    {}


    // Member Functions

        //- Return the density
        tmp<volScalarField> rho() const
        {
            return this->transport_.rho();
        }

        //- Return the heat capacity
        tmp<volScalarField> Cp() const
        {
            return this->transport_.Cp();
        }

        //- Return the conductivity
        tmp<volScalarField> lambda() const
        {
            return this->transport_.lambda();
        }

        //- Return the turbulent Prandtl number
        tmp<volScalarField> Prt() const
        {
            return this->transport_.Prt();
        }

        //- Return the laminar Prandtl number
        tmp<volScalarField> Pr() const
        {
            return this->transport_.Pr();
        }

        //- Return the laminar thermal diffusivity
        tmp<volScalarField> alphaLam() const
        {
            return this->transport_.alphaLam();
        }

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<volScalarField> pPrime() const;

        //- Return the face-phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<surfaceScalarField> pPrimef() const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff
        (
            volVectorField& U,
            const word& scheme
        ) const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            volVectorField& U,
            const word& scheme
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PhaseIncompressibleTurbulenceModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

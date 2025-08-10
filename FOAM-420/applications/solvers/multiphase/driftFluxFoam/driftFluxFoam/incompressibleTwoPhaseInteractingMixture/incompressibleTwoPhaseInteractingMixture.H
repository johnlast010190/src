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
    (c) 2014-2015 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

Class
    Foam::incompressibleTwoPhaseInteractingMixture

Description
    A two-phase incompressible transportModel for interacting phases
    requiring the direct evaluation of the mixture viscosity,
    e.g. activated sludge or slurry.

SourceFiles
    incompressibleTwoPhaseInteractingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhaseInteractingMixture_H
#define incompressibleTwoPhaseInteractingMixture_H

#include "compressibleTransportModel/compressibleTransportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "mixtureViscosityModel/mixtureViscosityModel.H"
#include "twoPhaseMixture/twoPhaseMixture.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class incompressibleTwoPhaseInteractingMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseInteractingMixture
:
    public IOdictionary,
    public compressibleTransportModel,
    public twoPhaseMixture
{
protected:

    // Protected data

        autoPtr<mixtureViscosityModel> muModel_;
        autoPtr<viscosityModel> nucModel_;

        dimensionedScalar rhod_;
        dimensionedScalar rhoc_;

        //- Optional diameter of the dispersed phase particles
        dimensionedScalar dd_;

        //- Optional maximum dispersed phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField mu_;


public:

    TypeName("incompressibleTwoPhaseInteractingMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoPhaseInteractingMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~incompressibleTwoPhaseInteractingMixture()
    {}


    // Member Functions

        //- Return const-access to the mixture viscosityModel
        const mixtureViscosityModel& muModel() const
        {
            return muModel_();
        }

        //- Return const-access to the continuous-phase viscosityModel
        const viscosityModel& nucModel() const
        {
            return nucModel_();
        }

        //- Return const-access to the dispersed-phase density
        const dimensionedScalar& rhod() const
        {
            return rhod_;
        }

        //- Return const-access to continuous-phase density
        const dimensionedScalar& rhoc() const
        {
            return rhoc_;
        };

        //- Return the diameter of the dispersed-phase particles
        const dimensionedScalar& dd() const
        {
            return dd_;
        }

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic mixture viscosity
        tmp<volScalarField> mu() const
        {
            return mu_;
        }

        //- Return the dynamic mixture viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
            return mu_.boundaryField()[patchi];
        }

        //- Return the mixture density
        virtual tmp<volScalarField> rho() const
        {
            return alpha1_*rhod_ + alpha2_*rhoc_;
        }

        //- Return the mixture density for patch
        virtual tmp<scalarField> rho(const label patchi) const
        {
            return
                alpha1_.boundaryField()[patchi]*rhod_.value()
              + alpha2_.boundaryField()[patchi]*rhoc_.value();
        }

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const
        {
            return mu_/rho();
        }

        //- Return the mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return mu_.boundaryField()[patchi]/rho(patchi);
        }

        //- Return the mixture heat capacity
        virtual tmp<volScalarField> Cp() const
        {
            NotImplemented;
        }

        //- Return the mixture conductivity
        virtual tmp<volScalarField> lambda() const
        {
            NotImplemented;
        }

        //- Return the mixture Prandtl number
        virtual tmp<volScalarField> Prt() const
        {
            NotImplemented;
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            mu_ = muModel_->mu(rhoc_*nucModel_->nu());
        }

        //- Read base transportProperties dictionary
        virtual bool read()
        {
            return readImpl();
        }

private:
    bool readImpl();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

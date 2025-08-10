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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::incompressible::RASModels::LamBremhorstKE

Group
    grpIcoRASTurbulence

Description
    Lam and Bremhorst low-Reynolds number k-epsilon turbulence model
    for incompressible flows

    This turbulence model is described in:
    \verbatim
        Lam, C. K. G., & Bremhorst, K. (1981).
        A modified form of the k-ε model for predicting wall turbulence.
        Journal of Fluids Engineering, 103(3), 456-460.
    \endverbatim

SourceFiles
    LamBremhorstKE.C

\*---------------------------------------------------------------------------*/

#ifndef LamBremhorstKE_H
#define LamBremhorstKE_H

#include "turbulentTransportModels/turbulentTransportModel.H"
#include "eddyViscosity/eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class LamBremhorstKE Declaration
\*---------------------------------------------------------------------------*/

class LamBremhorstKE
:
    public eddyViscosity<incompressible::RASModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        LamBremhorstKE(const LamBremhorstKE&);
        void operator=(const LamBremhorstKE&);

        tmp<volScalarField> Rt() const;
        tmp<volScalarField> fMu(const volScalarField& Rt) const;
        tmp<volScalarField> f1(const volScalarField& fMu) const;
        tmp<volScalarField> f2(const volScalarField& Rt) const;
        void correctNut(const volScalarField& fMu);


protected:

    // Protected data

        dimensionedScalar Cmu_;
        dimensionedScalar Ceps1_;
        dimensionedScalar Ceps2_;
        dimensionedScalar sigmaEps_;

        volScalarField k_;
        volScalarField epsilon_;

        //- Wall distance
        //  Note: different to wall distance in parent RASModel
        //  which is for near-wall cells only
        //const volScalarField& y_; // moved to each function call
        // this is done to auto re-create meshObject should it be deleted
        // i.e. during load balancing


    // Protected Member Functions

        virtual void correctNut();


public:

    //- Runtime type information
    TypeName("LamBremhorstKE");


    // Constructors

        //- Construct from components
        LamBremhorstKE
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~LamBremhorstKE()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

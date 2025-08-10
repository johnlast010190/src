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
    Foam::incompressible::RASModels::kkLOmega

Group
    grpIcoRASTurbulence

Description
    Low Reynolds-number k-kl-omega turbulence model for
    incompressible flows.

    This turbulence model is described in:
    \verbatim
        Walters, D. K., & Cokljat, D. (2008).
        A three-equation eddy-viscosity model for Reynolds-averaged
        Navier–Stokes simulations of transitional flow.
        Journal of Fluids Engineering, 130(12), 121401.
    \endverbatim

    however the paper contains several errors which must be corrected for the
    model to operation correctly as explained in

    \verbatim
        Furst, J. (2013).
        Numerical simulation of transitional flows with laminar kinetic energy.
        Engineering MECHANICS, 20(5), 379-388.
    \endverbatim

    All these corrections and updates are included in this implementation.

    The default model coefficients are
    \verbatim
        kkLOmegaCoeffs
        {
            A0             4.04
            As             2.12
            Av             6.75
            Abp            0.6
            Anat           200
            Ats            200
            CbpCrit        1.2
            Cnc            0.1
            CnatCrit       1250
            Cint           0.75
            CtsCrit        1000
            CrNat          0.02
            C11            3.4e-6
            C12            1.0e-10
            CR             0.12
            CalphaTheta    0.035
            Css            1.5
            CtauL          4360
            Cw1            0.44
            Cw2            0.92
            Cw3            0.3
            CwR            1.5
            Clambda        2.495
            CmuStd         0.09
            Prtheta        0.85
            Sigmak         1
            Sigmaw         1.17
        }
    \endverbatim

SourceFiles
    kkLOmega.C

\*---------------------------------------------------------------------------*/

#ifndef kkLOmega_H
#define kkLOmega_H

#include "RAS/RASModel/RASModel.H"
#include "eddyViscosity/eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kkLOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kkLOmega
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private member functions

        // Disallow default bitwise copy construct and assignment
        kkLOmega(const kkLOmega&);
        void operator=(const kkLOmega&);

        tmp<volScalarField> fv(const volScalarField& Ret) const;

        tmp<volScalarField> fINT() const;

        tmp<volScalarField> fSS(const volScalarField& omega) const;

        tmp<volScalarField> Cmu(const volScalarField& S) const;

        tmp<volScalarField> BetaTS(const volScalarField& Rew) const;

        tmp<volScalarField> fTaul
        (
            const volScalarField& lambdaEff,
            const volScalarField& ktL,
            const volScalarField& omega
        ) const;

        tmp<volScalarField> alphaT
        (
            const volScalarField& lambdaEff,
            const volScalarField& fv,
            const volScalarField& ktS
        ) const;

        tmp<volScalarField> fOmega
        (
            const volScalarField& lambdaEff,
            const volScalarField& lambdaT
        ) const;

        tmp<volScalarField> phiBP(const volScalarField& omega) const;

        tmp<volScalarField> phiNAT
        (
            const volScalarField& ReOmega,
            const volScalarField& fNatCrit
        ) const;

        tmp<volScalarField> D(const volScalarField& k) const;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar A0_;
            dimensionedScalar As_;
            dimensionedScalar Av_;
            dimensionedScalar Abp_;
            dimensionedScalar Anat_;
            dimensionedScalar Ats_;
            dimensionedScalar CbpCrit_;
            dimensionedScalar Cnc_;
            dimensionedScalar CnatCrit_;
            dimensionedScalar Cint_;
            dimensionedScalar CtsCrit_;
            dimensionedScalar CrNat_;
            dimensionedScalar C11_;
            dimensionedScalar C12_;
            dimensionedScalar CR_;
            dimensionedScalar CalphaTheta_;
            dimensionedScalar Css_;
            dimensionedScalar CtauL_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar CwR_;
            dimensionedScalar Clambda_;
            dimensionedScalar CmuStd_;
            dimensionedScalar Prtheta_;
            dimensionedScalar Sigmak_;
            dimensionedScalar Sigmaw_;


        // Fields

            volScalarField kt_;
            volScalarField kl_;
            volScalarField omega_;
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

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("kkLOmega");


    // Constructors

        //- Construct from components
        kkLOmega
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kkLOmega()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaT/Sigmak_ + this->nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaT/Sigmaw_ + this->nu())
            );
        }

        //- Return the laminar kinetic energy
        virtual tmp<volScalarField> kl() const
        {
            return kl_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> kt() const
        {
            return kt_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the total fluctuation kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "k",
                        this->mesh_.time().timeName(),
                        this->db_
                    ),
                    kt_ + kl_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the total fluctuation kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Validate the turbulence fields after construction
        //  Update turbulence viscosity and other derived fields as requires
        virtual void validate();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kkLOmega.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

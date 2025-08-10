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
    (c) 2019 Esi Ltd.

Class
    Foam::RASModels::ellipticBlendingLagKE

Group
    grpRASTurbulence

Description
    The model applies the Stress Strain Lag method for the k-epsilon model.
    Using the elliptic blending relaxation form of the Reynolds Stress Model
    (RSM) terms appearing in the derivation of the equations.

    The variant implemented employs N=6, such that f=0 on walls.

    Wall boundary conditions are:

        k       = kLowReWallFunction
        epsilon = epsilonWallFunction
        v2      = v2WallFunction
        f       = fWallFunction

    These are applicable to both low- and high-Reynolds number flows.

    Inlet values can be approximated by:

        v2      = 2/3 k
        f       = zero-gradient

    References:
    \verbatim
    Lardeau, Sylvain & Billard, Flavien. (2016).
    Development of an elliptic-blending lag model for industrial applications.
    54th AIAA Aerospace Sciences Meeting. 10.2514/6.2016-1600.

    \endverbatim

    The default model coefficients are
    \verbatim
        ellipticBlendingKECoeffs
        {
            Cmu         0.22;
            CmuKEps     0.09;
            C1          1.4;
            C2          0.3;
            CL          0.23;
            Ceta        70;
            Ceps2       1.9;
            Ceps3       -0.33;
            sigmaEps    1.3;
            sigmaK      1;
            sigmaPhi
        }
    \endverbatim

Note
    If the kLowReWallFunction is employed, a velocity variant of the turbulent
    viscosity wall function should be used, e.g. nutUWallFunction.  Turbulence
    k variants (nutk...) for this case will not behave correctly.

SourceFiles
    ellipticBlendingLagKE.C

\*---------------------------------------------------------------------------*/

#ifndef ellipticBlendingLagKE_H
#define ellipticBlendingLagKE_H

#include "Base/ellipticBlendingKE/ellipticBlendingKEBase.H"
#include "RAS/RASModel/RASModel.H"
#include "eddyViscosity/eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                             Class ellipticBlendingLagKE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class ellipticBlendingLagKE
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>,
    public ellipticBlendingKEBase
{

protected:

    // Protected data

        // Model coefficients
            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar C3star_;
            dimensionedScalar C4_;
            dimensionedScalar C5_;
            dimensionedScalar CL_;
            dimensionedScalar CT_;
            dimensionedScalar CtLim_;
            dimensionedScalar Ceta_;
            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Ceps3_;
            dimensionedScalar Ceps4_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;


        // Fields

            //- Turbulence kinetic energy
            volScalarField k_;

            //- Turbulence dissipation
            volScalarField epsilon_;

            //- Turbulence stress normal to streamlines
            volScalarField v2k_;

            //- Damping function
            volScalarField a_;

            //- Rotational source in v2k
            volScalarField rotSource_;

            //- Nut over k
            volScalarField nutOverk_;

            volScalarField fh2_;
        // Bounding values
            dimensionedScalar v2kMin_;
            dimensionedScalar v2kMax_;
            dimensionedScalar aMin_;

        // Number of fixed point iterations
            dimensioned<label> fixedPointIterations_;

        // v2k equations source bound for stabilization
            Switch boundV2kSource_;

    // Protected Member Functions

        virtual void correctNut(const volScalarField& S2);

        virtual void correctNut(const volTensorField& gradU);

        virtual void correctNut();

        //- Return time scale, Ts
        tmp<volScalarField> Ts() const;

        //- Return length scale, Ls
        tmp<volScalarField> Ls() const;

        //- Return time scale, TsLim
        tmp<volScalarField> TsLim(const volTensorField& gradU) const;

        //- Solve turbulence equations
        void solveTurb(const volTensorField& gradU);
public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("ellipticBlendingLagKE");


    // Constructors

        //- Construct from components
        ellipticBlendingLagKE
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
    virtual ~ellipticBlendingLagKE()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nut_/sigmaK_ + this->nu()/2.0
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    this->nut_/sigmaEps_ + this->nu()/2.0
                )
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

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> v2k() const
        {
            return v2k_;
        }

        //- Return the damping function
        virtual tmp<volScalarField> a() const
        {
            return a_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ellipticBlendingLagKE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

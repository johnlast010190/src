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
    (c) 2010-2016 Esi Ltd.
    (c) 2015-2016 OpenCFD Ltd.

Class
    Foam::LESModels::SpalartAllmarasDES

Group
    grpDESIurbulence

Description
    SpalartAllmarasDES DES turbulence model for incompressible and
    compressible flows

    Reference:
    \verbatim
        Spalart, P. R., Jou, W. H., Strelets, M., & Allmaras, S. R. (1997).
        Comments on the feasibility of LES for wings, and on a hybrid
        RANS/LES approach.
        Advances in DNS/LES, 1, 4-8.
    \endverbatim

    Including the low Reynolds number correction documented in
    \verbatim
        Spalart, P. R., Deck, S., Shur, M.L., Squires, K.D., Strelets, M.Kh,
        Travin, A. (2006).
        A new version of detached-eddy simulation, resistant to ambiguous grid
        densities.
        Theor. Comput. Fluid Dyn., 20, 181-195.
    \endverbatim

    Extended according to
    \verbatim
        "An Unstructured Grid Generation and Adaptive Solution Technique
        for High Reynolds Number Compressible Flows"
        G.A. Ashford,
        Ph.D. thesis, University of Michigan, 1996.
    \endverbatim
    by using the optional flag \c ashfordCorrection

SourceFiles
    SpalartAllmarasDES.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmarasDES_H
#define SpalartAllmarasDES_H

#include "DES/DESModel/DESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class SpalartAllmarasDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasDES
:
    public DESModel<BasicTurbulenceModel>
{
    // Private Member Functions

        // Check if is compressible
        bool isIsochoric() const;

        // Disallow default bitwise copy construct and assignment
        SpalartAllmarasDES(const SpalartAllmarasDES&);
        void operator=(const SpalartAllmarasDES&);


protected:

    // Protected data

        // Model constants

            dimensionedScalar sigmaNut_;
            dimensionedScalar sigmaNu_;
            dimensionedScalar kappa_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cs_;
            dimensionedScalar CDES_;
            dimensionedScalar ck_;
            Switch ft2Term_;
            Switch gradUstd_;

            // Low Reynolds number correction

                Switch lowReCorrection_;
                dimensionedScalar Ct3_;
                dimensionedScalar Ct4_;
                dimensionedScalar fwStar_;

        // Fields

            volScalarField nuTilda_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            //const volScalarField& y_; // moved to each function call
            // this is done to auto re-create meshObject should it be deleted
            // i.e. during load balancing


    // Protected Member Functions

        // return velocity gradient with wall function correction
        tmp<volTensorField> gradUnwc() const;

        tmp<volScalarField> chi() const;

        virtual tmp<volScalarField> Cb1
        (
            const volScalarField& Omega
        ) const;
        virtual tmp<volScalarField> Cw1
        (
            const volScalarField& Omega
        ) const;

        tmp<volScalarField> fv1(const volScalarField& chi) const;

        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        tmp<volScalarField> ft2(const volScalarField& chi) const;

        tmp<volScalarField> Omega(const volTensorField& gradU) const;

        tmp<volScalarField> Stilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volScalarField& Omega,
            const volScalarField& dTilda
        ) const;

        tmp<volScalarField> r
        (
            const volScalarField& nur,
            const volScalarField& Stilda,
            const volScalarField& dTilda
        ) const;

        tmp<volScalarField> fw
        (
            const volScalarField& Stilda,
            const volScalarField& dTilda
        ) const;

        tmp<volScalarField> psi
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU
        ) const;

        //- Length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU
        ) const;

        void correctNut(const volScalarField& fv1);
        virtual void correctNut();

        geometricZeroField gradNuTildaDotGradRho
        (
            const volVectorField& gradNuTilda,
            const geometricOneField& rho
        ) const;

        tmp<volScalarField> gradNuTildaDotGradRho
        (
            const volVectorField& gradNuTilda,
            const volScalarField& rho
        ) const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmarasDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasDES
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
    virtual ~SpalartAllmarasDES()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const;

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const;

        tmp<volScalarField> nuTilda() const
        {
            return nuTilda_;
        }

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;

        //- Correct nuTilda and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2018-2019 Esi Ltd.
    (c) 2015-2017 OpenFOAM Foundation
    (c) D. Segersson

    Copyright is held by the author.

Class
    Foam::fv::dalpeMassonCanopySource

Group
    grpFvOptionsSources

Original Author
    David Segersson,
    Swedish Meteorological and Hydrological Institute,
    david.segersson@smhi.se

Description

    Momentum source term is expressed as:
        \f[
            S_u = -\rho C_d \alpha u^2
        \f]

    Source term for k is expressed as:

        \f[
            S_k = \rho C_d \alpha \left(\beta_p u^3 - \beta_d k u \right)
        \f]

    Source term for epsilon is expressed as:

        \f[
            S_{\epsilon}=\rho C_d \alpha \frac{\epsilon}{k} \left(C_{\epsilon 4}
        \beta_p u^3 - C_{\epsilon 5} \beta_d k u \right)
        \f]

    whereas for omega

        \f[
            S_{\omega}=S_\frac{{\epsilon}{k \beta_star}
        \f]

    where
    \beta_p (1.0), \betad_d (5.03), C_{\epsilon 4} (0.78), C_{epsilon 5} (0.78)
    and \beta_star (1) are dimensionless coefficients

    \vartable
        \rho           | density
        u              | magnitude of velocity
        k              | turbulent kinetic energy
        \alpha         | leaf area density [m]
        C_d            | tree canopy drag coefficient
    \endvartable

    Use with any variant of the k-\epsilon and k-\omega turbulence model.

    Reference:
       Momentum and turbulence source/sink-terms for tree canopy according to
       B. Dalpé and C. Masson. Numerical simulation of wind flow near a forest
       edge. J. of Wind Eng. and Ind. Aerodyn, 2006


Usage
    Example usage:
    \verbatim
    canopy
    {
        type                       dalpeMassonCanopySource;
        active                     true;
    dalpeMassonCanopySourceCoeffs
    {
            selectionMode   cellZone;
        cellZone        tree;
        fields          ( U k epsilon );
        origin          ( 0 0 1 ); // at the tree base
        direction       ( 0 0 1 ); // along the tree height
        Cd              0.2;   // drag coefficient of trees
        LAI             1; // Leaf Area Index (used if LADmax is not specified)
        LADmax          4; Maximum Leaf Area Density [m ^-1]
        LADProfile ( 0.05 0.1 0.15 0.35 1.1 0.9 0.5 0.2 0.15 0.05 0.01 );
    }
    }

    // Cd drag coefficient
    // LADProfile is the vertical profile of Leaf Area Density
    // first value is closest to ground.
    // Each value represents an equal share of the tree height
    // e.g. for a 4 m tree and 4 values, each value will represent 1 m.
    // Values are scaled so that the highest will correspond to LADMax.

    // The user needs to set:
    // 1) LADProfile and LADmax, or
    // 2) LADProfile and LAI (LADmax calculated accordingly)

    // The origin is a point at the base of the tree, with the
    // tree height measured along the direction

    \endverbatim

SourceFiles
    dalpeMassonCanopySource.C

\*---------------------------------------------------------------------------*/

#ifndef dalpeMassonCanopySource_H
#define dalpeMassonCanopySource_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "dimensionedTypes/dimensionedTypes.H"
#include "sources/derived/vegetationSource/vegetationSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class dalpeMassonCanopySource Declaration
\*---------------------------------------------------------------------------*/

class dalpeMassonCanopySource
:
    public vegetationSource
{
private:

        // Canopy coefficients
        dimensionedScalar betaP_;
        dimensionedScalar betaD_;
        dimensionedScalar C4_;
        dimensionedScalar C5_;

        //- Disallow default bitwise copy construct
        dalpeMassonCanopySource(const dalpeMassonCanopySource&);

        //- Disallow default bitwise assignment
        void operator=(const dalpeMassonCanopySource&);


        //- Source term to momentum equation
        template<class RhoFieldType>
        void addSup
        (
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        template<class RhoFieldType>
        void addSup
        (
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );


        // Private Member Functions


public:


    //- Runtime type information
    TypeName("dalpeMassonCanopySource");


    // Constructors

        //- Construct from components
        dalpeMassonCanopySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~dalpeMassonCanopySource()
    {}


    // Member Functions
        // access

       dimensionedScalar betaP() const
        {
            return betaP_;
        }

       dimensionedScalar betaD() const
        {
            return betaD_;
        }

       dimensionedScalar C4() const
        {
            return C4_;
        }

       dimensionedScalar C5() const
        {
            return C5_;
        }

        // add source terms to momentum equations
        void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );


        //- Source term to scalar equations (k, epsilon, omega, T, w )
        void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Source term to compressible scalar equation (k, epsilon, omega, h, w )
        void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

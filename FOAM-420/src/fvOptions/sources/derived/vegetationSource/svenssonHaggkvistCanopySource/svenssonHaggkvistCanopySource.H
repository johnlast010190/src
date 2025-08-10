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
    (c) 2015-2017 OpenFOAM Foundation D. Segersson is held by the author.

    Copyright is held by the author.

Class
    Foam::fv::svenssonHaggkvistCanopySource

Group
    grpFvOptionsSources

Original Author
    David Segersson,
    Swedish Meteorological and Hydrological Institute,
    david.segersson@smhi.se

Description
    Momentum and turbulence source/sink-terms for tree canopy according
    to Svensson & Häggkvist (1990), with parameters adjusted
    according to Mochida et al. (2008).

    Momentum source term is expressed as:
        \f[
            S_u = -\rho C_d \alpha u^2
        \f]

    Source term for k is expressed as:

        \f[
            S_k = \rho C_d \alpha u^3
        \f]

    Source term for epsilon is expressed as:

        \f[
            S_{\epsilon} = \rho C_d \alpha C_{p \epsilon 1} \frac{\epsilon}{k} u^3
        \f]

    whereas for omega

        \f[
            S_{\omega}=S_\frac{{\epsilon}{k \beta_star}
        \f]

    where

    \vartable
        \rho             | density
        u                | magnitude of velocity
        k                | turbulent kinetic energy
        \alpha           | leaf area density [m]
        C_d              | tree canopy drag coefficient
        C_{p \epsilon 1} | coefficient with value 1.8 (Mochida et al. 2008)
    \endvartable

    Use with any variant of the k-\epsilon and k-\omega turbulence model.


    Reference:
        Mochida A., Tabata Y., Iwata T., Yoshino H., Examining tree canopy models
    for CFD prediction of wind environment at pedestrian level, J. Wind Eng. Ind.
        Aerodyn. 96 (2008) 1667-1677.

        Svensson, U., Häggkvist, K. A two-equation turbulence model for canopy flows.
    J. Wind Eng. Ind. Aerodyn. 35 (1990) 201–211.

Usage
    Example usage:
    \verbatim
    canopy
    {
        type                       svenssonHaggkvistCanopySource;
        active                     true;
    svenssonHaggkvistCanopySourceCoeffs
    {
            selectionMode   cellZone;
        cellZone        tree;
        fields          ( U k epsilon );
        origin          ( 0 0 1 );
        direction       ( 0 0 1 );
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
    svenssonHaggkvistCanopySource.C

\*---------------------------------------------------------------------------*/

#ifndef svenssonHaggkvistCanopySource_H
#define svenssonHaggkvistCanopySource_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "dimensionedTypes/dimensionedTypes.H"
#include "sources/derived/vegetationSource/vegetationSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class svenssonHaggkvistCanopySource Declaration
\*---------------------------------------------------------------------------*/

class svenssonHaggkvistCanopySource
:
    public vegetationSource
{
private:

        // Canopy coefficients
        dimensionedScalar CpEps1_;

        //- Disallow default bitwise copy construct
        svenssonHaggkvistCanopySource(const svenssonHaggkvistCanopySource&);

        //- Disallow default bitwise assignment
        void operator=(const svenssonHaggkvistCanopySource&);


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
    TypeName("svenssonHaggkvistCanopySource");


    // Constructors

        //- Construct from components
        svenssonHaggkvistCanopySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~svenssonHaggkvistCanopySource()
    {}


    // Member Functions
        // access

       dimensionedScalar CpEps1() const
        {
            return CpEps1_;
        }

        //- Source term to momentum equation
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


        //- Source term to compressible momentum equation
        void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
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

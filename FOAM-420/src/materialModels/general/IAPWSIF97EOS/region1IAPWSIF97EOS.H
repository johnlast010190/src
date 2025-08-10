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
    (c) 2022 Esi Ltd.

Class
    Foam::region1IAPWSIF97EOS

SourceFiles
    region1IAPWSIF97EOS.C

\*---------------------------------------------------------------------------*/

#ifndef region1IAPWSIF97EOS_H
#define region1IAPWSIF97EOS_H

#include "memory/autoPtr/autoPtr.H"
#include "db/regIOobject/regIOobject.H"
#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class region1IAPWSIF97EOS Declaration
\*---------------------------------------------------------------------------*/

class region1IAPWSIF97EOS
{
    // Private member data

        //- Ii tabulated data
        const scalarField I_;

        //- Ji tabulated data
        const scalarField J_;

        //- ni tabulated data
        const scalarField n_;

        //- Reference pressure [Pa]
        const scalar pStar_;

        //- Reference temperature [K]
        const scalar TStar_;

        //- Specific gas constant for ordinary water [J*kg^-1*K^-1]
        const scalar R_;


public:

    // Constructors

        //- Construct form the components
        region1IAPWSIF97EOS(const dictionary& dict);


    // Member functions

        //- Pressure ratio
        scalar phiRatio(const scalar p) const;

        //- Temperature ratio
        scalar tauRatio(const scalar T) const;

        //- The dimensionless Gibbs free energy
        scalar gamma(const scalar phi, const scalar tau) const;

        //- Gibbs energy derivitive d(gamma)/d(phi)
        scalar gammaPhi(const scalar phi, const scalar tau) const;

        //- Gibbs energy derivitive d^2(gamma)/d(phi^2)
        scalar gammaPhiPhi(const scalar phi, const scalar tau) const;

        //- Gibbs energy derivitive d(gamma)/d(tau)
        scalar gammaTau(const scalar phi, const scalar tau) const;

        //- Gibbs energy derivitive d^2(gamma)/d(tau^2)
        scalar gammaTauTau(const scalar phi, const scalar tau) const;

        //- Gibbs energy derivitive d^2(gamma)/d(tau)d(phi)
        scalar gammaPhiTau(const scalar phi, const scalar tau) const;

        //- Specific volume [m^3/kg]
        scalar v(const scalar p, const scalar T) const;

        //- Specific internal energy [J/kg]
        scalar u(const scalar p, const scalar T) const;

        //- Specific entropy [J/kgK]
        scalar s(const scalar p, const scalar T) const;

        //- Specific enthalpy [J/kg]
        scalar h(const scalar p, const scalar T) const;

        //- Specific isochoric heat capacity [J/kgK]
        scalar Cv(const scalar p, const scalar T) const;

        //- Specific isobaric heat capacity [J/kg*K]
        scalar Cp(const scalar p, const scalar T) const;

        //- Speed of sound [m/s]
        scalar w(const scalar p, const scalar T) const;

        //- Density [kg/m^3]
        scalar rho(const scalar p, const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

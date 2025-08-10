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
    Foam::region3IAPWSIF97EOS

SourceFiles
    region3IAPWSIF97EOS.C

\*---------------------------------------------------------------------------*/

#ifndef region3IAPWSIF97EOS_H
#define region3IAPWSIF97EOS_H

#include "memory/autoPtr/autoPtr.H"
#include "db/regIOobject/regIOobject.H"
#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class region3IAPWSIF97EOS Declaration
\*---------------------------------------------------------------------------*/

class region3IAPWSIF97EOS
{
    // Private member data

        //- Ii tabulated data
        const scalarField I_;

        //- Ji tabulated data
        const scalarField J_;

        //- ni tabulated data
        const scalarField n_;

        //- Critical temperature [K]
        const scalar Tc_;

        //- Critical density [kg/m^3]
        const scalar rhoc_;

        //- Specific gas constant for ordinary water [J*kg^-1*K^-1]
        const scalar R_;


public:

    // Constructors

        //- Construct form the components
        region3IAPWSIF97EOS(const dictionary& dict);


    // Member functions

        //- Compute density corresponing to pressure/temperature point
        scalar computeRho(const scalar p, const scalar T);

        //- Density ratio
        scalar deltaRatio(const scalar rho) const;

        //- Temperature ratio
        scalar tauRatio(const scalar T) const;

        //- The dimensionless Gibbs free energy
        scalar PHI(const scalar delta, const scalar tau) const;

        //- Gibbs energy derivitive d(PHI)/d(delta)
        scalar PHIDelta(const scalar delta, const scalar tau) const;

        //- Gibbs energy derivitive d^2(PHI)/d(delta^2)
        scalar PHIDeltaDelta(const scalar delta, const scalar tau) const;

        //- Gibbs energy derivitive d(PHI)/d(tau)
        scalar PHITau(const scalar delta, const scalar tau) const;

        //- Gibbs energy derivitive d^2(PHI)/d(tau^2)
        scalar PHITauTau(const scalar delta, const scalar tau) const;

        //- Gibbs energy derivitive d^2(PHI)/d(tau)d(delta)
        scalar PHIDeltaTau(const scalar delta, const scalar tau) const;

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

        //- Joule-Thomson coefficient
        scalar mu(const scalar p, const scalar T) const;

        //- Isothermal throttling coefficient
        scalar deltaT(const scalar p, const scalar T) const;

        //- Isentropic temperature-pressure coefficient
        scalar betaS(const scalar p, const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

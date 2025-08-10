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
    Foam::IAPWSIF97EOS

Description
    This calculation page is based on the Thermodynamic Properties of Water
    and Steam (Reference) [1] (named as IAPWS-IF97) provided by IAPWS.
    Detailed information about used equations, constants, range of validity
    etc is presented in PDF version of IAPWS.

    Release which can be downloaded from IAPWS web site www.iapws.org

Reference
    Revised Release on the IAPWS Industrial Formulation 1997 for
    the Thermodynamic Properties of Water and Steam.
    Lucerne, Switzerland, August 2007, available at http://www.iapws.org

SourceFiles
    IAPWSIF97EOS.C

\*---------------------------------------------------------------------------*/

#ifndef IAPWSIF97EOS_H
#define IAPWSIF97EOS_H

#include "fvMesh/fvMesh.H"
#include "region1IAPWSIF97EOS.H"
#include "region2IAPWSIF97EOS.H"
#include "region3IAPWSIF97EOS.H"
#include "region5IAPWSIF97EOS.H"
#include "global/etcFiles/etcFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IAPWSIF97EOS Declaration
\*---------------------------------------------------------------------------*/

class IAPWSIF97EOS
{
    // Private member data

        //- Data tables
        const dictionary dict_;

        //- Coefficients for saturation curve
        const scalarField n_;

        //- Coefficients for boundary between regions 2 and 3
        const scalarField n23_;

        //- Region1 function set
        region1IAPWSIF97EOS region1_;

        //- Region2 function set
        region2IAPWSIF97EOS region2_;

        //- Region3 function set
        region3IAPWSIF97EOS region3_;

        //- Region5 function set
        region5IAPWSIF97EOS region5_;


        // Viscosity coefficients

            //- Reference to mu dictionary
            const dictionary& muDict_;

            //- Maximum temperature up to 100MPa
            const scalar muTmax1_;

            //- Mximum temperature up to 50MPa
            const scalar muTmax2_;

            //- n0 coefficients
            const scalarField muN0_;

            //- I exponents
            const scalarField muI_;

            //- J exponents
            const scalarField muJ_;

            //- n coefficients
            const scalarField muN_;


        // Conductivity coefficients

            //- Reference to kappa dictionary
            const dictionary& kappaDict_;

            //- n coefficients
            const scalarField lambda0n_;

            //- n coefficients
            const Field<scalarField> lambda1n_;

            //- n coefficients
            const scalarField lambda2n_;

            //- n coefficients
            const scalarField lambda2n1_;

            //- Validity delta1
            const scalar delta1_;

            //- n coefficients
            const scalarField lambda2n2_;

            //- Validity delta2
            const scalar delta2_;

            //- n coefficients
            const scalarField lambda2n3_;

            //- Validity delta3
            const scalar delta3_;

            //- n coefficients
            const scalarField lambda2n4_;

            //- Validity delta4
            const scalar delta4_;

            //- n coefficients
            const scalarField lambda2n5_;

            //- Critical pressure
            const scalar pc_;

            //- Gas constant for kappa model
            const scalar kappaR_;


    // Private member functions

        //- Function loading dictionary from etc directory
        dictionary loadDataDict(const objectRegistry& obr) const;

        //- Load conductivity coefficients
        Field<scalarField> loadLambda1Coeffs(const dictionary& dict) const;


public:

    // Constructors

        //- Construct form the components
        IAPWSIF97EOS(const objectRegistry& obr);

        //- Quick model access
        enum {region1, region2, region3, region4, region5, nRegions_};


    // Member functions

        //- Calculate temperature or pressure on saturation curve
        scalar saturation
        (
            scalar ps,
            scalar Ts = 273.15,
            bool calcTemp = true
        ) const;

        //- Calculate boundary between regions 2 and 3
        scalar boundary23
        (
            scalar p,
            scalar T = 273.15,
            bool calcTemp = true
        ) const;

        //- Return which region we are computing based on pressure/temperature
        //  point
        label whichRegion(const scalar p, const scalar T) const;

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

        //- Isobaric compressibility [Pa^-1]
        scalar KT(const scalar p, const scalar T) const;

        //- Dynamic viscosity [kg/ms]
        scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        scalar kappa(const scalar p, const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

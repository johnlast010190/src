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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

Class
    Foam::functionTIsoSolidTransport

Description
    Non-constant properties (w.r.t. temperature) Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    functionTIsoSolidTransportI.H
    functionTIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionTIsoSolidTransport_H
#define functionTIsoSolidTransport_H

#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class functionTIsoSolidTransport;

template<class Thermo>
inline functionTIsoSolidTransport<Thermo> operator*
(
    const scalar,
    const functionTIsoSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const functionTIsoSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                   Class functionTIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class functionTIsoSolidTransport
:
    public Thermo
{
    // Private data

        //- Isotropic thermal conductivity
        autoPtr<Function1<scalar>> kappa_;


    // Private Member Functions

        //- Construct from components
        inline functionTIsoSolidTransport
        (
            const Thermo& t,
            const Function1<scalar>& kappa
        );


public:

    // Constructors

        //- Construct as named copy
        inline functionTIsoSolidTransport
        (
            const word&,
            const functionTIsoSolidTransport&
        );

        //- Construct from dictionary
        functionTIsoSolidTransport(const objectRegistry& obr, const dictionary& dict);

        // Selector from dictionary
        inline static autoPtr<functionTIsoSolidTransport> New
        (
            const objectRegistry& obr,
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "functionTIso<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

        //- Isotropic thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Un-functionTIsotropic thermal conductivity [W/mK]
        inline vector Kappa(const scalar p, const scalar T) const;

         //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;


        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const functionTIsoSolidTransport&);
        inline void operator+=(const functionTIsoSolidTransport&);


    // Friend operators

        friend functionTIsoSolidTransport operator* <Thermo>
        (
            const scalar,
            const functionTIsoSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const functionTIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transport/functionT/functionTIsoSolidTransportI.H"

#ifdef NoRepository
    #include "transport/functionT/functionTIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

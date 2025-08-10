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
    Foam::functionTTransport

Group
    grpSpecieTransport

Description
    Non-constant properties (w.r.t. temperature) Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    functionTTransportI.H
    functionTTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionTTransport_H
#define functionTTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class functionTTransport;

template<class Thermo>
inline functionTTransport<Thermo> operator+
(
    const functionTTransport<Thermo>&,
    const functionTTransport<Thermo>&
);

template<class Thermo>
inline functionTTransport<Thermo> operator*
(
    const scalar,
    const functionTTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const functionTTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class functionTTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class functionTTransport
:
    public Thermo
{
    // Private data

        //- Dynamic viscosity [Pa.s]
        autoPtr<Function1<scalar>> mu_;

        //- Thermal conductivity
        autoPtr<Function1<scalar>> kappa_;


    // Private Member Functions

        //- Construct from components
        inline functionTTransport
        (
            const Thermo& t,
            const Function1<scalar>& mu,
            const Function1<scalar>& kappa
        );


public:

    // Constructors

        //- Construct as named copy
        inline functionTTransport(const word&, const functionTTransport&);

        //- Construct from dictionary
        functionTTransport(const objectRegistry& mesh, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<functionTTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<functionTTransport> New
        (
            const objectRegistry& mesh,
            const dictionary& dict
        );

        //- Return a pointer to a new phaseModel created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew
            (
                const fvMesh& mesh
            )
            :
                mesh_(mesh)
            {}

            autoPtr<functionTTransport> operator()(const dictionary& dict) const
            {
                return New(mesh_, dict);
            }
        };

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "functionT<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const functionTTransport&);

        inline void operator+=(const functionTTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend functionTTransport operator+ <Thermo>
        (
            const functionTTransport&,
            const functionTTransport&
        );

        friend functionTTransport operator* <Thermo>
        (
            const scalar,
            const functionTTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const functionTTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transport/functionT/functionTTransportI.H"

#ifdef NoRepository
    #include "transport/functionT/functionTTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

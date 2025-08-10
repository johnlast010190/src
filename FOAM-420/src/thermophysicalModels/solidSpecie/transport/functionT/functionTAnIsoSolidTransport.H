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
    Foam::functionTAnIsoSolidTransport

Description
    Non-constant properties (w.r.t. temperature) Transport package.
    Templated into a given Thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    functionTAnIsoSolidTransportI.H
    functionTAnIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionTAnIsoSolidTransport_H
#define functionTAnIsoSolidTransport_H
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class functionTAnIsoSolidTransport;

template<class Thermo>
inline functionTAnIsoSolidTransport<Thermo> operator*
(
    const scalar,
    const functionTAnIsoSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const functionTAnIsoSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                   Class functionTAnIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class functionTAnIsoSolidTransport
:
    public Thermo
{
    // Private data

        //- Constant anisotropic thermal conductivity.
        autoPtr<Function1<vector>> kappa_;


    // Private Member Functions

        //- Construct from components
        inline functionTAnIsoSolidTransport
        (
            const Thermo& t,
            const Function1<vector>& kappa
        );


public:

    // Constructors

        //- Copy constructor
        functionTAnIsoSolidTransport(const functionTAnIsoSolidTransport&) = default;

        //- Construct as named copy
        inline functionTAnIsoSolidTransport
        (
            const word&,
            const functionTAnIsoSolidTransport&
        );

        //- Construct from dictionary
        functionTAnIsoSolidTransport(const objectRegistry& objectRegistry, const dictionary&);

        // Selector from dictionary
        inline static autoPtr<functionTAnIsoSolidTransport> New
        (
            const objectRegistry& objectRegistry,
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

            autoPtr<functionTAnIsoSolidTransport> operator()(const dictionary& dict) const
            {
                return New(mesh_, dict);
            }
        };

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "functionTAnIso<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = false;

        //- Isotropic thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Un-isotropic thermal conductivity [W/mK]
        inline vector Kappa(const scalar p, const scalar T) const;

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline vector alphah(const scalar p, const scalar T) const;


        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const functionTAnIsoSolidTransport&);
        inline void operator+=(const functionTAnIsoSolidTransport&);


    // Friend operators

        friend functionTAnIsoSolidTransport operator* <Thermo>
        (
            const scalar,
            const functionTAnIsoSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const functionTAnIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transport/functionT/functionTAnIsoSolidTransportI.H"

#ifdef NoRepository
    #include "transport/functionT/functionTAnIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Class
    Foam::constIsoSolidTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constIsoSolidTransportI.H
    constIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constIsoSolidTransport_H
#define constIsoSolidTransport_H

#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class constIsoSolidTransport;

template<class Thermo>
inline constIsoSolidTransport<Thermo> operator*
(
    const scalar,
    const constIsoSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const constIsoSolidTransport<Thermo>&
);

class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class constIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constIsoSolidTransport
:
    public Thermo
{
    // Private data

        //- Constant isotropic thermal conductivity
        scalar kappa_;


    // Private Member Functions

        //- Construct from components
        inline constIsoSolidTransport(const Thermo& t, const scalar kappa);


public:

    // Constructors

        //- Copy constructor
        constIsoSolidTransport(const constIsoSolidTransport&) = default;

        //- Construct as named copy
        inline constIsoSolidTransport
        (
            const word&,
            const constIsoSolidTransport&
        );

        //- Construct from dictionary
        constIsoSolidTransport(const objectRegistry& obr, const dictionary& dict);

        // Selector from dictionary
        inline static autoPtr<constIsoSolidTransport> New
        (
            const objectRegistry& mesh,
            const dictionary& dict
        );

        //- Return a pointer to a new phaseModel created on freestore
        //  from Istream
        class iNew
        {
            const objectRegistry& mesh_;

        public:

            iNew
            (
                const objectRegistry& mesh
            )
            :
                mesh_(mesh)
            {}

            autoPtr<constIsoSolidTransport> operator()(const dictionary& dict) const
            {
                return New(mesh_, dict);
            }
        };

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "constIso<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

        //- Isotropic thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Un-isotropic thermal conductivity [W/mK]
        inline vector Kappa(const scalar p, const scalar T) const;

         //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;


        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const constIsoSolidTransport&);
        inline void operator+=(const constIsoSolidTransport&);


    // Friend operators

        friend constIsoSolidTransport operator* <Thermo>
        (
            const scalar,
            const constIsoSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const constIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transport/const/constIsoSolidTransportI.H"

#ifdef NoRepository
    #include "transport/const/constIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

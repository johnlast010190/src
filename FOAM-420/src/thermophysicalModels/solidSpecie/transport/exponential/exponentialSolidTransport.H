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
    Foam::exponentialSolidTransport

Description
    Exponential properties for solid heat transport
    Templated into a given thermodynamics package.

SourceFiles
    exponentialSolidTransportI.H
    exponentialSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialSolidTransport_H
#define exponentialSolidTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class exponentialSolidTransport;

template<class Thermo>
inline exponentialSolidTransport<Thermo> operator*
(
    const scalar,
    const exponentialSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const exponentialSolidTransport<Thermo>&
);

class fvMesh;

/*---------------------------------------------------------------------------*\
                  Class exponentialSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class exponentialSolidTransport
:
    public Thermo
{
    // Private data

        //- Constant thermal coefficient.
        scalar kappa0_;

        //- Exponent coefficient
        scalar n0_;

        //- Reference temperature
        scalar Tref_;


    // Private Member Functions

        //- Construct from components
        inline exponentialSolidTransport
        (
            const Thermo& t,
            const scalar kappa0,
            const scalar n0,
            const scalar Tref
        );


public:

    // Constructors

        //- Copy constructor
        exponentialSolidTransport(const exponentialSolidTransport&) = default;


        //- Construct as named copy
        inline exponentialSolidTransport
        (
            const word&,
            const exponentialSolidTransport&
        );

        //- Construct from dictionary
        exponentialSolidTransport(const objectRegistry& obr, const dictionary&);

        // Selector from dictionary
        inline static autoPtr<exponentialSolidTransport> New
        (
            const objectRegistry& obr,
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

            autoPtr<exponentialSolidTransport> operator()(const dictionary& dict) const
            {
                return New(mesh_, dict);
            }
        };

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "exponential<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline vector Kappa(const scalar p, const scalar T) const;

         //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const exponentialSolidTransport&);
        inline void operator+=(const exponentialSolidTransport&);


    // Friend operators

        friend exponentialSolidTransport operator* <Thermo>
        (
            const scalar,
            const exponentialSolidTransport&
        );

    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const exponentialSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transport/exponential/exponentialSolidTransportI.H"

#ifdef NoRepository
    #include "transport/exponential/exponentialSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

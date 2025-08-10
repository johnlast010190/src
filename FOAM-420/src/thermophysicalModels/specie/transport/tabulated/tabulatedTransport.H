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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::tabulatedTransport

Group
    grpSpecieTransport

Description

    Transport properties using tabulated mu (Pa.s) and kappa (
    defined as:

    transport
    {
        //- required
        muTableCoeffs
        {
            fileName        "$FOAM_CASE/constant/mu.dat";
            outOfBounds     clamp;
        }

        //- required only if kappaLookup true
        kappaTableCoeffs
        {
            fileName        "$FOAM_CASE/constant/kappa.dat";
            outOfBounds     clamp;
        }
    }

    The mu and kappa tables are described in terms of
    the pressure and temperature
    For example, the table:

                           temperature (K)
                       |  300   450   600
                  -----+-----------------
                  1e5  |   A     B     C
    pressure (Pa) 1.1e5|   D     E     F
                  2e5  |   G     H     I


    Is specified by the following:

       (
            (
                1e5
                (
                    (300    A)
                    (450    B)
                    (600    C)
                )
            )
            (
                1.1e5
                (
                    (300    D)
                    (450    E)
                    (600    F)
                )
            )
            (
                2e5
                (
                    (300    G)
                    (450    H)
                    (600    I)
                )
            )
        )

    All pressures and temperatures must be in ascending order
    in each list.

SourceFiles
    tabulatedTransportI.H
    tabulatedTransport.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedTransport_H
#define tabulatedTransport_H

#include "interpolations/interpolation2DTable/interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class tabulatedTransport;

template<class Thermo>
inline tabulatedTransport<Thermo> operator+
(
    const tabulatedTransport<Thermo>&,
    const tabulatedTransport<Thermo>&
);

template<class Thermo>
inline tabulatedTransport<Thermo> operator*
(
    const scalar,
    const tabulatedTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const tabulatedTransport<Thermo>&
);

class fvMesh;

/*---------------------------------------------------------------------------*\
                    Class tabulatedTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulatedTransport
:
    public Thermo
{
    // Private data


    //- dictionary of mu table settings
    dictionary muDict_;

    //- dictionary of kappa table settings
    dictionary kappaDict_;

    //- 2D look up table viscosity = function of pressure and temperature
    autoPtr<interpolation2DTable<scalar>> muTable_;

    //- 2D look up table kappa = function of pressure and temperature
    autoPtr<interpolation2DTable<scalar>> kappaTable_;

    // Private Member Functions

        //- Construct from components
        inline tabulatedTransport
        (
            const Thermo& t,
            const dictionary&  muDict_,
            const dictionary&  kappaDict_,
            const autoPtr<interpolation2DTable<scalar>>& muTable_,
            const autoPtr<interpolation2DTable<scalar>>& kappaTable_
        );


public:

    // Constructors

        //- Copy constructor
        tabulatedTransport(const tabulatedTransport&) = default;

        //- Construct as named copy
        inline tabulatedTransport(const word&, const tabulatedTransport&);

        //- Construct from dictionary
        tabulatedTransport(const objectRegistry& obr, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<tabulatedTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedTransport> New(const objectRegistry& obr, const dictionary& dict);

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

            autoPtr<tabulatedTransport> operator()(const dictionary& dict) const
            {
                return New(mesh_, dict);
            }
        };

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulated<" + Thermo::typeName() + '>';
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

        inline void operator=(const tabulatedTransport&);

        inline void operator+=(const tabulatedTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend tabulatedTransport operator+ <Thermo>
        (
            const tabulatedTransport&,
            const tabulatedTransport&
        );

        friend tabulatedTransport operator* <Thermo>
        (
            const scalar,
            const tabulatedTransport&
        );

    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const tabulatedTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transport/tabulated/tabulatedTransportI.H"

#ifdef NoRepository
    #include "transport/tabulated/tabulatedTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

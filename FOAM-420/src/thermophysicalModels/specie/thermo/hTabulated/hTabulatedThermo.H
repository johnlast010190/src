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
    Foam:: hTabulatedThermo

Group
    grpSpecieThermo

Description
    Tabulated properties thermodynamics package
    templated into the EquationOfState.

    thermodynamics
    {
        //- Heat capacity table as function of P and T (required)
        CpTableCoeffs
        {
            fileName        "$FOAM_CASE/constant/CpTableName.dat";
            outOfBounds     clamp;
        }

        //- Heat of formation (required)
        Hf 1901;
    }

    The Cp tables are described in terms of
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
     hTabulatedThermoI.H
     hTabulatedThermo.C

\*---------------------------------------------------------------------------*/

#ifndef  hTabulatedThermo_H
#define  hTabulatedThermo_H

#include "interpolations/interpolation2DTable/interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class  hTabulatedThermo;

template<class EquationOfState>
inline  hTabulatedThermo<EquationOfState> operator+
(
    const  hTabulatedThermo<EquationOfState>&,
    const  hTabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
inline  hTabulatedThermo<EquationOfState> operator*
(
    const scalar,
    const  hTabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
inline  hTabulatedThermo<EquationOfState> operator==
(
    const  hTabulatedThermo<EquationOfState>&,
    const  hTabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const  hTabulatedThermo<EquationOfState>&
);

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class  hTabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class  hTabulatedThermo
:
    public EquationOfState
{
    // Private data

        //- dictionary of Cp table settings
        dictionary CpDict_;

        //- 2D look up table Cp = function of pressure and temperature
        //- Specific heat at constant pressure [J/(kg.K)]
        autoPtr<interpolation2DTable<scalar>> CpTable_;

        //- Heat of formation
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline hTabulatedThermo
        (
            const EquationOfState& st,
            const dictionary& CpDict,
            const autoPtr<interpolation2DTable<scalar>>& CpTable,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from dictionary
        hTabulatedThermo(const objectRegistry& obr, const dictionary& dict);

        //- Construct as named copy
        inline hTabulatedThermo(const word&, const  hTabulatedThermo&);

        //- Construct and return a clone
        inline autoPtr<hTabulatedThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hTabulatedThermo> New(const objectRegistry& obr, const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hTabulated<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const hTabulatedThermo&);


    // Friend operators

        friend hTabulatedThermo operator+ <EquationOfState>
        (
            const hTabulatedThermo&,
            const hTabulatedThermo&
        );

        friend hTabulatedThermo operator* <EquationOfState>
        (
            const scalar,
            const hTabulatedThermo&
        );

        friend hTabulatedThermo operator== <EquationOfState>
        (
            const hTabulatedThermo&,
            const hTabulatedThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hTabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermo/hTabulated/hTabulatedThermoI.H"

#ifdef NoRepository
    #include "thermo/hTabulated/hTabulatedThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

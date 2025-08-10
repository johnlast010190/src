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
    Foam::hFunctionTThermo

Group
    grpSpecieThermo

Description
    Thermodynamics package templated on the equation of state, using Function1
    functions for \c cp, \c h and \c s.

    Functions for \c h and \c s derived from \c cp.

Usage

    \table
        Property     | Description
        Hf           | Heat of formation
        Sf           | Standard entropy
        Cp           | Specific heat at constant pressure function
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        Hf              0;
        Sf              0;
        Cp              polynomial ( ( 1000 0 ) (-0.05 1) (0.003 2) );
    }
    \endverbatim

    The polynomial expression is evaluated as so:

        \f[
            Cp = 1000 - 0.05 T + 0.003 T^2
        \f]

Note
    - Heat of formation is inputted in [J/kg], but internally uses [J/kmol]
    - Standard entropy is inputted in [J/kg/K], but internally uses [J/kmol/K]
    - Specific heat at constant pressure polynomial coefficients evaluate to an
      expression in [J/(kg.K)].

SourceFiles
    hFunctionTThermoI.H
    hFunctionTThermo.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef hFunctionTThermo_H
#define hFunctionTThermo_H

#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState>
class hFunctionTThermo;

template<class EquationOfState>
inline hFunctionTThermo<EquationOfState> operator+
(
    const hFunctionTThermo<EquationOfState>&,
    const hFunctionTThermo<EquationOfState>&
);

template<class EquationOfState>
inline hFunctionTThermo<EquationOfState> operator*
(
    const scalar,
    const hFunctionTThermo<EquationOfState>&
);

template<class EquationOfState>
inline hFunctionTThermo<EquationOfState> operator==
(
    const hFunctionTThermo<EquationOfState>&,
    const hFunctionTThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hFunctionTThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                      Class hFunctionTThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hFunctionTThermo
:
    public EquationOfState
{
    // Private data

        //- Heat of formation
        scalar Hf_;

        //- Standard entropy
        scalar Sf_;

        //- Specific heat at constant pressure
        autoPtr<Function1<scalar>> Cp_;

    // Private Member Functions

        //- Construct from components
        inline hFunctionTThermo
        (
            const EquationOfState& pt,
            const scalar Hf,
            const scalar Sf,
            const Function1<scalar>& Cp
        );


public:

    // Constructors

        //- Construct from dictionary
        hFunctionTThermo(const objectRegistry& obr, const dictionary& dict);

        //- Construct as a named copy
        inline hFunctionTThermo(const word&, const hFunctionTThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hFunctionT<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


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

        inline void operator=(const hFunctionTThermo&);
        inline void operator+=(const hFunctionTThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend hFunctionTThermo operator+ <EquationOfState>
        (
            const hFunctionTThermo&,
            const hFunctionTThermo&
        );

        friend hFunctionTThermo operator* <EquationOfState>
        (
            const scalar,
            const hFunctionTThermo&
        );

        friend hFunctionTThermo operator== <EquationOfState>
        (
            const hFunctionTThermo&,
            const hFunctionTThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hFunctionTThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermo/hFunctionT/hFunctionTThermoI.H"

#ifdef NoRepository
    #include "thermo/hFunctionT/hFunctionTThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::constAnIsoSolidTransport

Description
    Constant properties Transport package.
    Templated into a given Thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constAnIsoSolidTransportI.H
    constAnIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constAnIsoSolidTransport_H
#define constAnIsoSolidTransport_H
#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class constAnIsoSolidTransport;

template<class Thermo>
inline constAnIsoSolidTransport<Thermo> operator*
(
    const scalar,
    const constAnIsoSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const constAnIsoSolidTransport<Thermo>&
);

class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class constAnIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constAnIsoSolidTransport
:
    public Thermo
{
    // Private data

        //- Constant anisotropic thermal conductivity.
        vector kappa_;


    // Private Member Functions

        //- Construct from components
        inline constAnIsoSolidTransport(const Thermo& t, const vector kappa);


public:

    // Constructors

        //- Copy constructor
        constAnIsoSolidTransport(const constAnIsoSolidTransport&) = default;

        //- Construct as named copy
        inline constAnIsoSolidTransport
        (
            const word&,
            const constAnIsoSolidTransport&
        );

        //- Construct from dictionary
        constAnIsoSolidTransport(const objectRegistry& obr, const dictionary&);

        // Selector from dictionary
        inline static autoPtr<constAnIsoSolidTransport> New
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

            autoPtr<constAnIsoSolidTransport> operator()(const dictionary& dict) const
            {
                return New(mesh_, dict);
            }
        };

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "constAnIso<" + Thermo::typeName() + '>';
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

        inline void operator=(const constAnIsoSolidTransport&);
        inline void operator+=(const constAnIsoSolidTransport&);


    // Friend operators

        friend constAnIsoSolidTransport operator* <Thermo>
        (
            const scalar,
            const constAnIsoSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const constAnIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transport/const/constAnIsoSolidTransportI.H"

#ifdef NoRepository
    #include "transport/const/constAnIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

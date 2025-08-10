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
    Foam::dummyThermo

\*---------------------------------------------------------------------------*/

#ifndef dummyThermo_H
#define dummyThermo_H

#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class dummyThermo;

dummyThermo operator+
(
    const dummyThermo&,
    const dummyThermo&
);

dummyThermo operator==
(
    const dummyThermo&,
    const dummyThermo&
);

dummyThermo operator*
(
    const scalar,
    const dummyThermo&
);

/*---------------------------------------------------------------------------*\
                      Class dummyThermo Declaration
\*---------------------------------------------------------------------------*/

class dummyThermo
{

protected:

    // Protected member data

        //- Dictionary
        dictionary dict_;

        //- Object registry
        const objectRegistry& obr_;

        //- Phase name
        word phaseName_;

        //- Specie name
        word specieName_;

        //- Specie mol weight
        scalar W_;


public:

    // Constructors

        //- Copy constructor
        dummyThermo(const dummyThermo&) = default;


        //- Construct from object registry, dictionary, phase name
        // specie name and object name
        dummyThermo
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        )
        :
            dict_(dict),
            obr_(obr),
            phaseName_(phaseName),
            specieName_(specieName),
            W_(dict_.optionalSubDict("specieCoeffs").lookup<scalar>("molWeight"))
        {}

        dummyThermo
        (
            const objectRegistry& obr,
            const dictionary& dict
        )
        :
            dict_(dict),
            obr_(obr),
            phaseName_(word::null),
            specieName_(dict.dictName()),
            W_(dict_.optionalSubDict("specieCoeffs").lookup<scalar>("molWeight"))
        {}

        dummyThermo
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName
        )
        :
            dict_(dict),
            obr_(obr),
            phaseName_(phaseName),
            specieName_(dict.dictName()),
            W_(dict_.optionalSubDict("specieCoeffs").lookup<scalar>("molWeight"))
        {}

        //- Construct as named copy
        dummyThermo(const word&, const dummyThermo& s1)
        :
            dict_(s1.dict_),
            obr_(s1.obr_),
            phaseName_(s1.phaseName_),
            specieName_(s1.specieName_),
            W_(dict_.optionalSubDict("specieCoeffs").lookup<scalar>("molWeight"))
        {}

        const fvMesh& mesh() const
        {
            if (isA<fvSolutionRegistry>(obr_))
            {
                return dynamic_cast<const fvSolutionRegistry&>(obr_).mesh();
            }
            return dynamic_cast<const fvMesh&>(obr_);
        }

        //- Construct and return a clone
        autoPtr<dummyThermo> clone() const
        {
            return autoPtr<dummyThermo>
            (
                new dummyThermo(*this)
            );
        }

        //- Selector from dictionary
        static autoPtr<dummyThermo> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        )
        {
            return autoPtr<dummyThermo>
            (
                new dummyThermo(obr, dict, phaseName, specieName, name)
            );
        }

    virtual ~dummyThermo(){};


    // Member Functions

        //- Is this correct?
        typedef dummyThermo thermoType;

        //- Return a pointer to a new phaseModel created on freestore
        //  from Istream
        class iNew
        {
            const objectRegistry& mesh_;

        public:

            iNew(const objectRegistry& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<dummyThermo> operator()(const dictionary& dict) const
            {
                return New(mesh_, dict, word::null, dict.dictName(), word::null);
            }
        };

        //- Return the instantiated type name
        // static word typeName;
        static word typeName()
        {
            return "dummyThermo";
        }

        word phaseName()
        {
            return phaseName_;
        }

        //- Molecular weight [kg/kmol]
        scalar W() const
        {
            return W_;
        }

        //- Molecular weight [kg/kmol]
        scalar Y() const
        {
            NotImplemented;
        }


    // Dummy member functions

        //- Temperature from absolute enthalpy
        //  given an initial temperature T0 [K]
        scalar THa(const scalar ha, const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Chemical enthalpy [J/kg]
        scalar Hc() const
        {
            NotImplemented;
        }

        //- Absolute Enthalpy [J/kg]
        scalar Ha(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Absolute Enthalpy [J/kmol]
        scalar ha(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Heat capacity at constant pressure [J/(kmol K)]
        scalar cp(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Heat capacity at constant pressure [J/(kg K)]
        scalar Cp(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Heat capacity at constant volume [J/(kg K)]
        scalar Cv(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Sensible enthalpy [J/kg]
        scalar Hs(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Entropy [J/(kg K)]
        scalar S(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Sensible internal energy [J/kg]
        scalar Es(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Gibbs free energy [J/kg]
        scalar G(const scalar p, const scalar T) const
        {
            NotImplemented;
        }
        scalar K(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Helmholtz free energy [J/kg]
        scalar A(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Enthalpy/Internal energy [J/kg]
        scalar HE(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Dynamic viscosity [kg/ms]
        scalar mu(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Thermal conductivity [W/mK]
        scalar kappa(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Thermal diffusivity of enthalpy [kg/ms]
        scalar alphah(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Return density [kg/m^3]
        scalar rho(const scalar p, const scalar T) const
        {
            NotImplemented;
        }

        //- Equilibrium constant i.t.o. molar concentration
        //  = PIi(ci/cstd)^nui
        //  For low pressures (where the gas mixture is near perfect)
        //  Kc = Kp(pstd/(RR*T))^nu
        scalar Kc(const scalar p, const scalar T) const
        {
            NotImplemented;
        }


    // Member operators

        //- Assignemt operator
        void operator=(const dummyThermo& c)
        {
            NotImplemented;
        }

        void operator+=(const dummyThermo&)
        {
            NotImplemented;
        }


    // Friend operators

        friend dummyThermo operator==
        (
            const dummyThermo& c1,
            const dummyThermo& c2
        )
        {
            NotImplemented;
        }

        friend dummyThermo operator+
        (
            const dummyThermo& c1,
            const dummyThermo& c2
        )
        {
            NotImplemented;
        }

        friend dummyThermo operator*
        (
            const scalar,
            const dummyThermo& c1
        )
        {
            NotImplemented;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2024 Esi Ltd.

Class
    Foam::basicThermo

Description
    Abstract base-class for fluid and solid thermodynamic properties

SourceFiles
    basicThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicThermo_H
#define basicThermo_H

#include "fields/volFields/volFields.H"
#include "db/typeInfo/typeInfo.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "primitives/strings/word/wordIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class materialTables;

/*---------------------------------------------------------------------------*\
                        Class basicThermo Declaration
\*---------------------------------------------------------------------------*/
class basicThermo
:
    public IOdictionary
{

protected:

    // Protected data

        //- Phase-name
        const word& phaseName_;

        //- Material type
        const word materialType_;


        // Fields

            //- Gauge pressure [Pa]
            volScalarField& p_;

            //- Reference pressure for gauge pressure
            scalar pRef_;

            //- Temperature [K]
            volScalarField& T_;

            //- Laminar thermal diffusuvity [kg/m/s]
            volScalarField alpha_;

        //- Should the dpdt term be included in the enthalpy equation
        Switch dpdt_;


    // Protected Member Functions

        //- Construct as copy (not implemented)
        basicThermo(const basicThermo&);

        //- Look up or construct field
        volScalarField& lookupOrConstruct
        (
            const objectRegistry& mesh,
            const word& name,
            const dimensionSet& dims = dimless,
            const scalar& defaultValue = 0
        ) const;

        volScalarField& lookupOrConstruct
        (
            const objectRegistry& mesh,
            const char* name,
            const dimensionSet& dims = dimless,
            const scalar& defaultValue = 0
        ) const;

        //- Look up or construct field; if phaseName is valid, try to construct
        //  from the field without phase name if phasic one is not found
        volScalarField& lookupOrConstructPhasic
        (
            const objectRegistry& mesh,
            const word& name,
            const dimensionSet& dims = dimless,
            const scalar& defaultValue = 0
        ) const;

        //- Return list of enthalpy/internal energy field boundary types
        //  by interrogating the temperature field boundary types
        wordList heBoundaryTypes();

        //- Return list of enthalpy/internal energy field boundary base types
        //  by interrogating the temperature field boundary types
        wordList heBoundaryBaseTypes();

        //- Lookup and check out field
        void lookupAndCheckout(const char* name) const;

        //- Load pRef from dict
        scalar readPref() const;

        const fvMesh& mesh() const;
public:

    // Static member functions

        //- Return the enthalpy/internal energy field boundary type
        //  by interrogating the temperature field boundary type
        static word heBoundaryType(const fvPatchScalarField& tpf);

        //- Return the enthalpy/internal energy field boundary base types
        //  by interrogating the temperature field boundary type
        static word heBoundaryBaseType(const fvPatchScalarField& tpf);

public:

    //- Runtime type information
    TypeName("basicThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicThermo,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );


    // Constructors

        //- Construct from obr and phase name
        basicThermo
        (
            const objectRegistry&,
            const word& phaseName
        );

        //- Construct from obr, dictionary and phase name
        basicThermo
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName
        );


    // Selectors

        //- Generic lookup for thermodynamics package thermoTypeName
        template<class Thermo, class Table>
        static typename Table::iterator lookupThermo
        (
            const dictionary& thermoTypeDict,
            Table& tablePtr,
            const int nCmpt,
            const char* cmptNames[],
            const word& thermoTypeName,
            const word& phaseName = word::null
        );

        //- Generic lookup for each of the related thermodynamics packages
        template<class Thermo, class Table>
        static typename Table::iterator lookupThermo
        (
            const dictionary& thermoDict,
            Table& tablePtr,
            const word& phaseName = word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const objectRegistry&,
            const word& phaseName=word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName=word::null
        );

        //- Specialisation of the Generic New for basicThermo
        static autoPtr<basicThermo> New
        (
            const objectRegistry&,
            const word& phaseName=word::null
        );

        //- Generic lookup for any thermo type
        static basicThermo* lookupPtr
        (
            const objectRegistry&, const word& phaseName=word::null
        );

        //- Generic lookup or create for any thermo type
        static basicThermo& lookupOrCreate
        (
            const objectRegistry&, const word& phaseName=word::null
        );

    //- Destructor
    ~basicThermo() override;


    // Member functions

        //- Dictionary name thermophysicalProperties/materialProperties
        static word dictName;

        //- Materials dict name
        static const word matDictName;

        //- Loads either thermophysical or material properties dict
        static IOobject loadDictionary
        (
            const objectRegistry& obr,
            const word& phaseName,
            const bool isRegistred
        );

        static word phasePropertyName
        (
            const word& name,
            const word& phaseName
        )
        {
            // Apparently there exhist Euler models that can contain different
            // pressures for each phase. By including all variables in
            // phasePropertyName function we can handle naming globally.
            if (name == "p" || name == matDictName)
            {
                return name;
            }
            return IOobject::groupName(name, phaseName);
        }

        word phasePropertyName(const word& name) const
        {
            return basicThermo::phasePropertyName(name, phaseName_);
        }

        static const basicThermo& lookupThermo(const fvPatchScalarField& pf);

        //- Return dictionary for phase
        const dictionary& phaseDict() const
        {
            if (phaseName_ == word::null)
            {
                return *this;
            }
            else
            {
                return this->subDict(phaseName_);
            }
        }

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&,
            const word&,
            const word&
        ) const;

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const word& thermoName,
            const int nCmpt
        );

        //- Override to perform additional init after heThermo's init is run
        virtual void init() {};

        //- Update properties
        virtual void correct() = 0;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const = 0;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const = 0;

        //- Return true if the Cp or Cv value (depending on selected energy
        //  type) is constant
        virtual bool isCpvConst() const
        {
            FatalErrorInFunction
                << "This solver requires the new material properties library."
                << nl
                << "Ensure 'system/materialProperties' is present."
                << endl << exit(FatalError);
            return false;
        }

        //- It shouldn't be used from here. It is overriden in
        // new material library
        // (but we don't want the function to be forced in old thermo derived classes)
        virtual bool buoyant() const
        {
            FatalErrorInFunction
                << "This solver requires the new material properties library."
                << nl
                << "Ensure 'system/materialProperties' is present."
                << endl << exit(FatalError);
            return false;
        }

        //- It the buoyant rho same as main rho?
        virtual bool distinctBuoyancy() const
        {
            FatalErrorInFunction
                << "This solver requires the new material properties library."
                << "Ensure 'system/materialProperties' is present."
                << endl << exit(FatalError);
            return false;
        }

        //- Should the dpdt term be included in the enthalpy equation
        Switch dpdt() const
        {
            return dpdt_;
        }

        //- Do we solve for T from he or vice-versa?
        virtual bool solvesTFromhe() const
        {
            return true;
        }

        //- Set to solve for T from he
        virtual void setSolveTFromhe()
        {}

        //- Set to solve for he from T
        virtual void setSolveheFromT()
        {
            FatalErrorInFunction
                << "This solver requires the new material properties library."
                << "Ensure 'system/materialProperties' is present."
                << endl << exit(FatalError);
        }


        // Access to thermodynamic state variables

            //- Gauge pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Return reference pressure for gauge pressure
            virtual dimensionedScalar pRef() const
            {
                return dimensionedScalar("pRef", p().dimensions(), pRef_);
            }

            virtual inline scalar pRefValue() const
            {
                return pRef_;
            }

            //- Pressure [Pa]
            virtual const volScalarField& p() const;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density [kg/m^3] (density for use in buoyancy term)
            virtual tmp<volScalarField> buoyantRho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he() = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const = 0;

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const = 0;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const = 0;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const = 0;


        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            virtual const volScalarField& T() const;

            //- Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& T();

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const = 0;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const = 0;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const = 0;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity ratio []
            virtual tmp<volScalarField> CpByCpv() const = 0;

            //- Heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;


        // Access to transport state variables

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual const volScalarField& alpha() const;

            //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
            virtual const scalarField& alpha
            (
                const label patchi
            ) const;


        // Fields derived from transport state variables

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappa() const = 0;

            //- Thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappa
            (
                const label patchi
            ) const = 0;

            //- Effective thermal diffusivity for temperature
            //  of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const = 0;

            //- Effective thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;

            //- Effective thermal diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Effective thermal diffusivity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;


        //- Read thermophysical properties dictionary
        virtual bool read() override;


    // Member operators

        //- Return material tables
        virtual materialTables& materials() const
        {
            FatalErrorInFunction
                << "This solver requires the new material properties library."
                << nl
                << "Ensure 'system/materialProperties' is present."
                << endl << exit(FatalError);
            NotImplemented;
        }

        //- Return calculated model
        tmp<volScalarField> operator[]
        (
            const word& modelName
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicThermo/basicThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

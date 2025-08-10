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
    (c) 2019 Esi Ltd.

Class
    Foam::multiComponentMixtureBase

Group
    grpReactionThermophysicalMixtures

Description
    Foam::multiComponentMixtureBase

SourceFiles
    multiComponentMixtureBase.C

\*---------------------------------------------------------------------------*/

#ifndef multiComponentMixtureBase_H
#define multiComponentMixtureBase_H

#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "speciesTable/speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiComponentMixtureBase Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multiComponentMixtureBase
{
protected:

    // Protected data

        //- Table of specie names
        const speciesTable& species_;

        //- Species mass fractions
        PtrList<volScalarField>& Y_;

        //- Species data
        PtrList<ThermoType> speciesData_;

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable ThermoType mixtureVol_;


    // Private Member Functions

        //- Construct as copy (not implemented)
        //multiComponentMixtureBase(const multiComponentMixtureBase<ThermoType>&);

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        virtual const ThermoType& constructSpeciesData(const objectRegistry& obr, const dictionary& thermoDict);

        //- Correct the mass fractions to sum to 1
        virtual void correctMassFractions();

public:

    //- Runtime type information
    TypeName("multiComponentMixtureBase");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            multiComponentMixtureBase,
            dictionary,
            (
                const dictionary& dict,
                const objectRegistry& obr,
                const speciesTable& species,
                PtrList<volScalarField>& Y,
                const word& phaseName
            ),
            (dict, obr, species, Y, phaseName)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            multiComponentMixtureBase,
            objectRegistry,
            (
                const dictionary& dict,
                const wordList& specieNames,
                const HashPtrTable<ThermoType>& thermoData,
                const objectRegistry& obr,
                const speciesTable& species,
                PtrList<volScalarField>& Y,
                const word& phaseName
            ),
            (dict, specieNames, thermoData, obr, species, Y, phaseName)
        );


    // Selectors

        //- Standard selection based onobjectRegistry
        static autoPtr<multiComponentMixtureBase> New
        (
            const dictionary& dict,
            const objectRegistry& obr,
            const speciesTable& species,
            PtrList<volScalarField>& Y,
            const word& phaseName=word::null
        );

        //- Selection based on components
        static autoPtr<multiComponentMixtureBase> New
        (
            const dictionary& dict,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const objectRegistry& obr,
            const speciesTable& species,
            PtrList<volScalarField>& Y,
            const word& phaseName=word::null
        );


    // Constructors

        //- Construct from dictionary, specie names, thermo database,
        //  obr and phase name
        multiComponentMixtureBase
        (
            const dictionary&,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const objectRegistry&,
            const speciesTable& species,
            PtrList<volScalarField>& Y,
            const word&
        );

        //- Construct from dictionary, obr and phase name
        multiComponentMixtureBase
        (
            const dictionary&,
            const objectRegistry&,
            const speciesTable& species,
            PtrList<volScalarField>& Y,
            const word&
        );


    //- Destructor
    virtual ~multiComponentMixtureBase()
    {}


    // Member functions

        virtual const ThermoType& cellMixture(const label celli) const = 0;

        virtual const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const = 0;

        virtual const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const = 0;

        virtual const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const = 0;

        //- Return the raw specie thermodynamic data
        const PtrList<ThermoType>& speciesData() const
        {
            return speciesData_;
        }

        //- Read dictionary
        virtual void read(const objectRegistry&, const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& getLocalThermo(const label speciei) const
        {
            return speciesData_[speciei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixtures/multiComponentMixtureBase/multiComponentMixtureBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeMultiComponentMixture(Thermo)                                               \
    defineTemplateTypeNameAndDebug(multiComponentMixtureBase<Thermo>, 0);               \
    defineTemplateRunTimeSelectionTable(multiComponentMixtureBase<Thermo>, dictionary)  \
    defineTemplateRunTimeSelectionTable(multiComponentMixtureBase<Thermo>, objectRegistry)

#define makeMultiComponentMixtureType(Mixture, Thermo)                                  \
    defineNamedTemplateTypeNameAndDebug(Mixture<Thermo>, 0);                            \
    multiComponentMixtureBase<Thermo>::adddictionaryConstructorToTable<Mixture<Thermo>> \
        add##Mixture##Thermo##ConstructorToTable_;                                      \
    multiComponentMixtureBase<Thermo>::addobjectRegistryConstructorToTable<Mixture<Thermo>>    \
        add##Mixture##Thermo##objectRegistryConstructorToTable_;

// for templated multi component mixture
#define addTemplateMultiComponentMixtureType(Mixture, Thermo)                           \
    defineNamedTemplateTypeNameAndDebug(Mixture, 0);                                    \
    multiComponentMixtureBase<Thermo>::adddictionaryConstructorToTable<Mixture>         \
        add##Mixture##Thermo##ConstructorToTable_;                                      \
    multiComponentMixtureBase<Thermo>::addobjectRegistryConstructorToTable<Mixture<Thermo>>    \
        add##Mixture##Thermo##objectRegistryConstructorToTable_;


#endif

// ************************************************************************* //

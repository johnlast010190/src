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
    (c) 2021-2023 Esi Ltd.

Class
    Foam::materialTables

SourceFiles
    materialTables.C
    materialTablesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef materialTables_H
#define materialTables_H

#include "materialModels/baseModels/materialModels.H"
#include "materialModels/materialTypedefs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class materialTables Declaration
\*---------------------------------------------------------------------------*/

class materialModel;

class materialTables
:
    public regIOobject
{
// Private data memmbers

    //- Helper object for collecting information
    //  for linking models
    struct linkingInfo
    {
        word funcName;
        word phaseName;
        word specieName;
        word linkToModel;
        word modelCastType;
    };

    //- List of model references
    HashTable<linkingInfo> linkingModelsInfo_;

    //- Links to models allowed?
    Switch enableLinkingModels_;

    //- Struct for the solver flags
    struct solverFlags
    {
        word funcName;
        word phaseName;
        word specieName;
        bool incompressible;
        bool isochoric;
        bool constCpv;
        bool isotropic;
        bool mixture;
    };

    //- Information about rho and kappa model flags
    //  (can't be accessed through model pointers at the moment)
    HashTable<solverFlags> solverFlagsInfo_;



// Private member functions

    //- Update table pointers
    template<class Type>
    void updateTable(const HashTable<Type>& table);

    //- Collect information about solver flags
    void collectSolverFlagsInfo
    (
        const word& phaseName,
        const word& specieName,
        const word& funcName,
        const materialModel* modelPtr
    );

    //- Set solver flags
    void setSolverFlags();

    //- Models for mixture
    wordList listModels
    (
        const word& phaseName,
        const wordList& groupModelNames = wordList::null()
    );

    //- Group error message
    void groupNotFoundError
    (
        const word& modelName,
        const word& noDictDepGroup
    );


protected:

    //- Scalar material models
    HashTable<matScalarTable> scalarModels_;

    //- Vector material models
    HashTable<matVectorTable> vectorModels_;

    //- Tensor material models
    HashTable<matTensorTable> tensorModels_;

    //- Object registry
    const objectRegistry& obr_;

    //- Materials object registry
    const objectRegistry& matObr_;

    //- Dictionary
    const dictionary& dict_;

    //- Is the equation of state is incompressible i.e. rho != f(p)
    HashTable<bool> incompressible_;

    //- Is the equation of state is isochoric i.e. rho = const
    HashTable<bool> isochoric_;

    //- Return true if Cp or Cv is uniform, depending on enthalpy/energy type
    HashTable<bool> constCpv_;

    //- Return true if thermal conductivity is isotropic
    HashTable<bool> isotropic_;

    //- Are material tables multi-phase
    bool isMultiphase_;

    //- Are material tables multi-species
    HashTable<bool> isMultiSpecies_;

    //- Type of the model
    enum typeOfTheModel
    {
        sType, // scalar model type
        vType, // vector model type
        tType  // tensor model type
    };


public:

    //- Runtime type information
    TypeName("materialTables");

    // Constructors

        //- Construct from components
        materialTables
        (
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~materialTables()
    {};


    // Public member functions

        //- Species key word
        static const word speciesListName;

        //- Phases key word
        static const word phasesListName;

        //- Return sub-table name
        word tableName
        (
            const word& phaseName,
            const word& specieName
        ) const;

        //- Check if model name found in the table
        template<class Type>
        bool foundModel
        (
            const word& tName,
            const word& modelName
        );

        //- Check if model is found in any of the tables
        bool foundModel
        (
            const word& tName,
            const word& modelName
        )
        {
            return
                foundModel<scalar>(tName, modelName)
             || foundModel<vector>(tName, modelName)
             || foundModel<tensor>(tName, modelName);
        }

        //- Table name of first species/phase in the mixture
        word firstSpeciesPhase(const word& phaseName);

        //- Extract names from the dict
        const dictionary& processDict
        (
            const word& phaseName,
            const word& specieName,
            const word& typeName,
            const word& funcName,
            word& modelType,
            word& obrModelName,
            const word& noDictModelType = word::null
        );

        //- Return reference to top level dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Are links between models enabled?
        bool linksEnabled() const
        {
            return enableLinkingModels_;
        }

        //- Return const reference to scalar table
        const HashTable<matScalarTable>& sTable() const;

        //- Return vector table
        const HashTable<matVectorTable>& vTable() const;

        //- Return tensor table
        const HashTable<matTensorTable>& tTable() const;

        //- Return reference to scalar table
        HashTable<matScalarTable>& sTable();

        //- Return reference to vector table
        HashTable<matVectorTable>& vTable();

        //- Return reference to tensor table
        HashTable<matTensorTable>& tTable();

        //- Return reference to scalar table based on phase/specie
        matScalarTable& sTable(const word& phaseName, const word& specieName);

        //- Return reference to vector table based on phase/specie
        matVectorTable& vTable(const word& phaseName, const word& specieName);

        //- Return reference to tensor table based on phase/specie
        matTensorTable& tTable(const word& phaseName, const word& specieName);

        //- Helper function to add scalar models into the hash table
        bool addScalarModel
        (
            const word& funcName,
            const word& phaseName,
            const word& specieName,
            materialModel* modelPtr,
            const word& castType
        );

        //- Helper function to add vector models into the hash table
        bool addVectorModel
        (
            const word& funcName,
            const word& phaseName,
            const word& specieName,
            materialModel* modelPtr,
            const word& castType
        );

        //- Helper function to add tensor models into the hash table
        bool addTensorModel
        (
            const word& funcName,
            const word& phaseName,
            const word& specieName,
            materialModel* modelPtr,
            const word& castType
        );

        //- Adding model
        void add
        (
            const word& typeName,
            const word& funcName,
            const word& phaseName = word::null,
            const word& specieName = word::null,
            const word& noDictModelType = word::null,
            const word& modelCastType = word::null
        );

        //- Group of models
        void addGroup
        (
            const wordList& modelNames,
            const wordList& funcNames,
            const word& phaseName = word::null,
            const word& specieName = word::null,
            const word& noDictModelType = word::null,
            bool isMixture = false,
            const word& mixModelName = word::null,
            const wordList& modelCastTypes = wordList::null()
        );

        //- Add models on species level
        bool addSpecies
        (
            const word& modelName,
            const wordList& funcNames,
            const word& phaseName,
            const word& noDictModelType,
            const word& noDictDepGroup,
            const wordList& modelCastTypes = wordList::null()
        );

        //- Add models on phase level
        bool addPhases
        (
            const word& modelName,
            const wordList& funcNames,
            const word& noDictModelType,
            const word& noDictDepGroup,
            const wordList& modelCastTypes = wordList::null()
        );

        //- Adding compositions of species
        bool addSpeciesMixtures
        (
            const word& modelName,
            const wordList& funcNames,
            const word& phaseName,
            const word& noDictModelType,
            const word& noDictDepGroup,
            const wordList& modelCastTypes = wordList::null(),
            bool groupExhist = false
        );

        //- Adding compositions of phases
        bool addPhaseMixtures
        (
            const word& modelName,
            const wordList& funcNames,
            const word& noDictModelType,
            const word& noDictDepGroup,
            const wordList& modelCastTypes = wordList::null(),
            bool groupExhist = false
        );

        //- Adding scalar models
        // Note:
        // Function is by default looking up entries from the dictionary.
        // If dictionary entry isn't specified and noDictModelType
        // it is picked up and added on the level of noDictDepGroup.
        template<class objectType>
        void addSpeciesAndSpeciesMixtures
        (
            const word& phaseName,
            const word& noDictModelType = word::null,
            const word& noDictDepGroup = word::null
        );

        //- Test for new adding structure
        void addSpeciesAndSpeciesMixtures
        (
            const word& modelName,
            const wordList& funcNames,
            const word& phaseName,
            const word& noDictModelType = word::null,
            const word& noDictDepGroup = word::null,
            const wordList& modelCastTypes = wordList::null()
        );

        //- Default mixture definition
        word standardMixture
        (
            const word& funcName,
            const word& phaseName,
            const word& defaultType = word::null
        );

        //- Check if model exists in the phase dictionary
        bool isModelInDict
        (
            const word& modelName, const word& phaseName=word::null
        ) const;

        //- Check model dependencies
        virtual void checkDepenencies();

        //- Add references to the models
        //  and add mixtures if they are not present through linking
        void addModelReferences();

        //- Link models and update tables all the dependent pointers
        void linkModelsAndUpdateTables();

        //- Update field
        void updateScalarField
        (
            const word& fieldName,
            const volScalarField& field
        );

        void updateVectorField
        (
            const word& fieldName,
            const volVectorField& field
        );

        void updateTensorField
        (
            const word& fieldName,
            const volTensorField& field
        );

        void updateSubDictsPtrs();

        //- Is the phase incompressible?
        bool incompressible(const word& phaseName)
        {
            if (incompressible_.found(phaseName))
            {
                return incompressible_[phaseName];
            }

            return false;
        }

        //- Is the mixture incompressible?
        bool incompressible()
        {
            bool incompressible = true;
            for
            (
                HashTable<bool>::iterator iter = incompressible_.begin();
                iter != incompressible_.end();
                ++iter
            )
            {
                incompressible &= iter();
            }
            return incompressible;
        }

        //- Is the phase isochoric?
        bool isochoric(const word& phaseName)
        {
            if (isochoric_.found(phaseName))
            {
                return isochoric_[phaseName];
            }

            return false;
        }

        //- Is the mixture isochoric?
        bool isochoric()
        {
            // If more than one material, we can assume their rho's are not
            // equal and therefore not constant
            return !(isochoric_.size() > 1);
        }


        //- Is the phase Cpv constant?
        bool isCpvConst(const word& phaseName)
        {
            if (constCpv_.found(phaseName))
            {
                return constCpv_[phaseName];
            }

            return false;
        }

        //- Is does the phase mixture have const Cpv?
        bool isCpvConst()
        {
            // If more than one material, we can assume their Cpv is not equal
            // and therefore not constant
            return !(constCpv_.size() > 1);
        }


        //- Is the phase isotropic?
        bool isotropic(const word& phaseName)
        {
            if (isotropic_.found(phaseName))
            {
                return isotropic_[phaseName];
            }

            return false;
        }

        //- Is the mixture isotropic?
        bool isotropic()
        {
            bool isotropic = true;
            for
            (
                HashTable<bool>::iterator iter = isotropic_.begin();
                iter != isotropic_.end();
                ++iter
            )
            {
                isotropic &= iter();
            }
            return isotropic;
        }


        //- Checks if the name is in the mixture models list
        bool isMixture(const word& funcName)
        {
            return
                (
                   baseModels<scalar>::isMixture(funcName)
                || baseModels<vector>::isMixture(funcName)
                || baseModels<tensor>::isMixture(funcName)
                );
        }


    // IO

        //- Write for data
        virtual bool writeData(Ostream& os) const override
        {
            NotImplemented;
        }

        //- Write
        virtual bool write(const bool valid) const override
        {
            NotImplemented;
        }

        //- Write
        virtual bool read() const;

        //- Print models information
        virtual void reportModelsInfo() const;

        //- Print model information for each type
        template<class Type>
        void reportModelsInfo(const HashTable<Type>& table) const;


    // Operators available for scalar quantities
    // Note: template call operators would look too complicated
    // const baseModels<scalar>& model =
    //     materials.operator()<scalar>(modelName);
    // Use vTable, tTable instead

        //- Access operator
        const baseModels<scalar>& operator()
        (
            const word& modelName,
            const word& phaseName,
            const word& specieName
        ) const;

        //- Specie mixtures acces operator
        const baseModels<scalar>& operator()
        (
            const word& modelName,
            const word& phaseName
        ) const;

        //- Phase mixtures acces operator
        const baseModels<scalar>& operator()
        (
            const word& modelName
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "materialTablesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

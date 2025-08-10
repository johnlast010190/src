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
    (c) 2016-2023 Esi Ltd.

Class
    Foam::stateFunction

Description
    Abstract base class for run-time selectable state specific functions
    Also contains functional code for state specialisation and assembly

SourceFiles
    stateFunction.C

\*---------------------------------------------------------------------------*/

#ifndef stateFunction_H
#define stateFunction_H

#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "primitives/enums/NamedEnum.H"
#include "modificationSwitches/modificationSwitches.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class stateIndex;

/*---------------------------------------------------------------------------*\
                          Class stateFunction Declaration
\*---------------------------------------------------------------------------*/

class stateFunction
{
public:

    class interfacePair
    :
        public Pair<word>
    {
    public:

        class hash
        :
            public Hash<interfacePair>
        {
        public:

            hash()
            = default;

            label operator()(const interfacePair& key) const
            {
                return word::hash()(key.first()) + word::hash()(key.second());
            }
        };


        // Constructors

            interfacePair()
            = default;

            interfacePair(const word& alpha1Name, const word& alpha2Name)
            :
                Pair<word>(alpha1Name, alpha2Name)
            {}

        // Friend Operators

            friend bool operator==
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return
                (
                    ((a.first() == b.first()) && (a.second() == b.second()))
                 || ((a.first() == b.second()) && (a.second() == b.first()))
                );
            }

            friend bool operator!=
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return (!(a == b));
            }
    };


    enum regionType
    {
        rtFluid, rtSolid, rtSurface, rtNone, rtMultipointAdjoint
    };
    enum timeType
    {
        ttSteady, ttTrans, ttNone
    };
    enum turbulenceType
    {
        tuRAS, tuLES, tuLam, tuNone
    };
    enum compressibilityType
    {
        ctIncomp, ctComp, ctNone
    };
    enum customSolver
    {
        csOff, csMixing, csTemp, csNone
    };
    enum matrixType
    {
        maSegregated, maCoupled
    };
    enum meshType
    {
        meVolume, mePoint, meSurface
    };

    enum fieldType
    {
        ftScalar, ftVector, ftTensor, ftSymmTensor, ftSphericalTensor
    };

    static const NamedEnum<regionType, 5> regionTypeNames_;
    static const NamedEnum<timeType, 3> timeTypeNames_;
    static const NamedEnum<turbulenceType, 4> turbulenceTypeNames_;
    static const NamedEnum<compressibilityType, 3> compressibilityTypeNames_;
    static const NamedEnum<customSolver, 4> customSolverNames_;
    static const NamedEnum<matrixType, 2> matrixTypeNames_;
    static const NamedEnum<meshType, 3> meshTypeNames_;
    static const NamedEnum<fieldType, 5> fieldTypeNames_;



    // Static member functions

        //- directory file path
        static fileName fileDir
        (
            const Time&,
            const fileName&,
            const word& regions = word::null
        );

        //- write directories to file
        static void writeDirDicts
        (
            const objectRegistry& obr,
            const fileName& region,
            const fileName& local,
            const dictionary& dicts,
            const wordList& excludedDicts = wordList::null(),
            const wordList& includedDicts = wordList::null()
        );

        //- find element in list
        static bool wordListFind(const wordList& wl, const word& w);

        //- check if two word lists contain the same entries
        static bool wordListMatch(const wordList& a, const wordList& b);

        //- read file from defaults and return dictionary
        static Xfer<dictionary> etcDictionary(const fileName&, bool);

        //- remove module names from state list and return module list
        static  wordList extractModules
        (
            const dictionary& modules,
            wordList& stateList
        );

        //- return state identifier
        static word stateIdentifier
        (
            const dictionary& input,
            const dictionary& defaults
        );

        //- find state identifier string from register of all possible states
        static word matchState
        (
            const PtrList<entry>& stateRegister,
            const wordList& stateList
        );


private:

    // Private Data

        //- region name
        word regionName_;

        //- mesh name
        word meshName_;

        //- reference to input dictionary
        dictionary input_;

        //- reference to defaults dictionary
        dictionary& defaults_;

        //- state identifier
        // can be different from type to provide
        // reuse potential via stateType entry in state defaults
        word stateId_;

        //- Modification switches
        modificationSwitches modSw_;

        // list of field names (store for reuse)
        autoPtr<wordList> fieldNamesPtr_;

protected:

    // Protected Data

        //- assembly dictionary
        dictionary stateDict_;

        //- dictionary IO options
        const bool distributed_;
        const bool collated_;

        const stateIndex& index_;

private:

    // Private Functions

        //- inject module states into stateDict
        void insertModules
        (
            const dictionary& defaults,
            const wordList& moduleIds,
            word stateId,
            dictionary& targetDict
        );

        void insertProgramaticSettings
        (
            const dictionary& defaults,
            const dictionary& input
        );

        Xfer<dictionary> readDictionaries
        (
            const dictionary&,
            const stateFunction&
        ) const;


        //- Disallow default bitwise copy construct
        stateFunction(const stateFunction&);

        //- Disallow default bitwise assignment
        void operator=(const stateFunction&);


protected:

    // Protected member functions


public:

    //- Runtime type information
        TypeName("stateFunction");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            stateFunction,
            dictionary,
            (
                word region,
                const dictionary& input,
                const dictionary& defaults,
                const stateFunction& master,
                const stateIndex& index,
                word meshName
            ),
            (region, input, defaults, master, index, meshName)
        );

    // Constructors

        //- Construct abstract master region from components
        stateFunction
        (
            const dictionary& input,
            const dictionary& defaults,
            const bool distributed,
            const bool collated,
            const stateIndex& index
        );

        //- Construct meshed region from components
        stateFunction
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

    // Selectors

        //- Return an autoPtr to the selected stateFunction
        static autoPtr<stateFunction> New
        (
            word region,
            const dictionary& input,
            dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

    //- Destructor
    virtual ~stateFunction()
    {}


    // Member Functions

        // Access

            // state property indicators

                //- region (fluid/solid/surface/none)
                virtual regionType region() const = 0;

                //- time (transient/steady/none)
                virtual timeType time() const = 0;

                //- turbulence (RANS/LES/none)
                virtual turbulenceType turbulence() const = 0;

                //- compressibility (compressible/incompressible/none)
                virtual compressibilityType compressibility() const = 0;

                //- coupling type (segregated/coupled)
                virtual matrixType matrix() const
                {
                    return maSegregated;
                };

                //- USF vs legacy solver state
                virtual bool usf() const
                {
                    return false;
                }


            // privat member access

                //- Const access to region name
                inline word regionName() const
                {
                    return regionName_;
                }

                inline word meshName() const
                {
                    return meshName_;
                }

                //- Const access to input dictionary
                inline const dictionary& input() const
                {
                    return input_;
                }

                //- Non-const access to input dictionary
                inline dictionary& input()
                {
                    return input_;
                }

                //- Const access to input dictionary
                inline const dictionary& defaults() const
                {
                    return defaults_;
                }

                //- Const access to input dictionary
                inline const word stateType() const
                {
                    return stateId_;
                }

                //- Access to modification switches
                inline const modificationSwitches& switches() const
                {
                    return modSw_;
                }



            //- const access to constant dictionaries
            const dictionary& constant() const;

            //- non-const access to constant dictionaries
            dictionary& constant();


            //- const access to system dictionaries
            const dictionary& system() const;

            //- non-const access to system dictionaries
            dictionary& system();

            //- const access to system dictionaries
            const dictionary& uniform() const;

            //- non-const access to system dictionaries
            dictionary& uniform();

            //- const access to function object list
            const dictionary& functions() const;

            //- const access fieldDefinitions
            const dictionary& fieldDefinitions() const;

            //- return list of fields that will be initialised
            // (requires fieldDefinitions to be fully populated)
            const wordList& fieldNames();


        // Check

        // Edit

            // Dictionary constructions

                //- assemble field definitions
                virtual void initialise();

                //- correct stateDict
                virtual void correct();

                //- inject const/system user input
                virtual void finalise();

                //- merge dictionary into stateDict
                void mergeToState(const dictionary&);

                //- edit controlDict and other global dictionaries
                void updateMaster(stateFunction&);

            // Update Switches

                //- reset boundary switch when only initialise fields
                void resetBoundarySwitch(const bool& reset);

            // Mesh and field construction

                //- generate meshes and store old fields
                virtual void createObjects(const Time& time) = 0;

                //- modify mesh
                virtual void modifyMesh(scalar writePause) = 0;

                //- create fields
                virtual void createFields() = 0;

        // Write

            //- write system and constant output dictionaries
            void writeAllDictionaries
            (
                const Time& runTime,
                const wordList& uniqueSystemDicts = wordList::null(),
                bool excludeUniqueDicts = false
            ) const;

            //- write system output dictionaries
            void writeSystemDictionaries
            (
                const Time& runTime,
                const wordList& uniqueDicts = wordList::null(),
                bool excludeUniqueDicts = false
            ) const;

            //- write constant output dictionaries
            void writeConstantDictionaries(const Time& runTime) const;

            //- write constant output dictionaries
            void writeUniformDictionaries(const Time& runTime) const;
    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

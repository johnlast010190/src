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

Class
    Foam::chemistryReader

Group
    grpReactionThermophysicalChemistryReaders

Description
    Abstract class for reading chemistry

SourceFiles
    chemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryReader_H
#define chemistryReader_H

#include "db/typeInfo/typeInfo.H"
#include "specieElement/specieElement.H"
#include "reaction/Reactions/Reaction/Reaction.H"
#include "reaction/Reactions/ReactionList/ReactionList.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "db/objectRegistry/objectRegistry.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef HashTable<List<specieElement>> speciesCompositionTable;

class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class chemistryReader Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistryReader
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        chemistryReader(const chemistryReader&);

        //- Disallow default bitwise assignment
        void operator=(const chemistryReader&);
protected:

        const fvMesh& getMesh(const objectRegistry& obr) const;
public:

    //- Runtime type information
    TypeName("chemistryReader");

    //- The type of thermo package the reader was instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct null
        chemistryReader()
        {}


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            chemistryReader,
            dictionary,
            (
                const dictionary& thermoDict,
                speciesTable& species,
                const objectRegistry& obr
            ),
            (thermoDict, species, obr)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<chemistryReader> New
        (
            const dictionary& thermoDict,
            speciesTable& species,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~chemistryReader()
    {}


    // Member Functions

        //- Return access to the list of species
        virtual const speciesTable& species() const = 0;

        //- Table of species composition
        virtual const speciesCompositionTable& specieComposition() const = 0;

        //- Return access to the thermo packages
        virtual const HashPtrTable<ThermoType>& speciesThermo() const = 0;

        //- Return access to the list of reactions
        virtual const ReactionList<ThermoType>& reactions() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryReaders/chemistryReader/chemistryReader.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#define makeChemistryReader(Thermo)                                            \
    defineTemplateTypeNameAndDebug(chemistryReader<Thermo>, 0);                \
    defineTemplateRunTimeSelectionTable(chemistryReader<Thermo>, dictionary)


#define makeChemistryReaderType(Reader, Thermo)                                \
    defineNamedTemplateTypeNameAndDebug(Reader<Thermo>, 0);                    \
    chemistryReader<Thermo>::adddictionaryConstructorToTable<Reader<Thermo>> \
        add##Reader##Thermo##ConstructorToTable_


// for non-templated chemistry readers
#define addChemistryReaderType(Reader, Thermo)                                 \
    defineTypeNameAndDebug(Reader, 0);                                         \
    chemistryReader<Thermo>::adddictionaryConstructorToTable<Reader>           \
        add##Reader##Thermo##ConstructorToTable_


// for templated chemistry readers
#define addTemplateChemistryReaderType(Reader, Thermo)                         \
    defineNamedTemplateTypeNameAndDebug(Reader, 0);                            \
    chemistryReader<Thermo>::adddictionaryConstructorToTable<Reader>           \
        add##Reader##Thermo##ConstructorToTable_


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

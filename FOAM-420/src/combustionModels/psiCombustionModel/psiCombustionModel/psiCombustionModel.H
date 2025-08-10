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
    Foam::psiCombustionModel

Description
    Combustion models for compressibility-based thermodynamics

SourceFiles
    psiCombustionModelI.H
    psiCombustionModel.C
    psiCombustionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef psiCombustionModel_H
#define psiCombustionModel_H

#include "combustionModel/combustionModel.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "psiReactionThermo/psiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                     class psiCombustionModel Declaration
\*---------------------------------------------------------------------------*/

class psiCombustionModel
:
    public combustionModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        psiCombustionModel(const psiCombustionModel&);

        //- Disallow default bitwise assignment
        void operator=(const psiCombustionModel&);


public:

    typedef psiReactionThermo ReactionThermo;


    //- Runtime type information
    TypeName("psiCombustionModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiCombustionModel,
        dictionary,
        (
            const word& modelType,
            const objectRegistry& obr,
            const word& combustionProperties,
            const word& phaseName
        ),
        (modelType, obr, combustionProperties, phaseName)
    );


    // Constructors

        //- Construct from components
        psiCombustionModel
        (
            const word& modelType,
            const objectRegistry& obr,
            const word& combustionProperties,
            const word& phaseName
        );


        //- Selector
        static autoPtr<psiCombustionModel> New
        (
            const objectRegistry& obr,
            const word& combustionProperties=combustionPropertiesName,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~psiCombustionModel();


    // Member Functions

        //- Return access to the thermo package
        virtual psiReactionThermo& thermo() = 0;

        //- Return const access to the thermo package
        virtual const psiReactionThermo& thermo() const = 0;


    // IO

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

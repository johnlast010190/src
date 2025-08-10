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
    Foam::rhoCombustionModel

Description
    Combustion models for rho-based thermodynamics

SourceFiles
    rhoCombustionModelI.H
    rhoCombustionModel.C
    rhoCombustionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef rhoCombustionModel_H
#define rhoCombustionModel_H

#include "combustionModel/combustionModel.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "rhoReactionThermo/rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                     class rhoCombustionModel Declaration
\*---------------------------------------------------------------------------*/

class rhoCombustionModel
:
    public combustionModel
{

    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoCombustionModel(const rhoCombustionModel&);

        //- Disallow default bitwise assignment
        void operator=(const rhoCombustionModel&);


public:

    typedef rhoReactionThermo ReactionThermo;


    //- Runtime type information
    TypeName("rhoCombustionModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoCombustionModel,
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
        rhoCombustionModel
        (
            const word& modelType,
            const objectRegistry& obr,
            const word& combustionProperties,
            const word& phaseName
        );



    //- Selector
    static autoPtr<rhoCombustionModel> New
    (
        const objectRegistry& obr,
        const word& combustionProperties=combustionPropertiesName,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~rhoCombustionModel();


    // Member Functions

        // Access functions

            //- Access combustion dict
             inline const dictionary& coeff() const;


        //- Return access to the thermo package
        virtual rhoReactionThermo& thermo() = 0;

        //- Return const access to the thermo package
        virtual const rhoReactionThermo& thermo() const = 0;


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

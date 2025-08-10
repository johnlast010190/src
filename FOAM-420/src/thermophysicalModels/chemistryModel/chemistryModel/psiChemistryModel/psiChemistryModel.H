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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::psiChemistryModel

Description
    Chemistry model for compressibility-based thermodynamics

SourceFiles
    psiChemistryModelI.H
    psiChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef psiChemistryModel_H
#define psiChemistryModel_H

#include "chemistryModel/basicChemistryModel/basicChemistryModel.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "psiReactionThermo/psiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class psiChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class psiChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        psiChemistryModel(const psiChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const psiChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<psiReactionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("psi");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiChemistryModel,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        psiChemistryModel(const objectRegistry& obr, const word& phaseName);


    //- Selector
    static autoPtr<psiChemistryModel> New
    (
        const objectRegistry& obr,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~psiChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline psiReactionThermo& thermo();

        //- Return const access to the thermo package
        inline const psiReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryModel/psiChemistryModel/psiChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

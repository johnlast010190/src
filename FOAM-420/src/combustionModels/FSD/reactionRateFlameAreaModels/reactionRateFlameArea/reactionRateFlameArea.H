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
    Foam::reactionRateFlameArea

Description
    Abstract class for reaction rate per flame area unit

SourceFiles
    reactionRateFlameArea.C
    reactionRateFlameAreaNew.C

\*---------------------------------------------------------------------------*/

#ifndef reactionRateFlameArea_H
#define reactionRateFlameArea_H

#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "db/dictionary/dictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "fields/volFields/volFields.H"
#include "combustionModel/combustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class reactionRateFlameArea Declaration
\*---------------------------------------------------------------------------*/

class reactionRateFlameArea
{

protected:

    // Protected data

        //- Dictionary
        dictionary coeffDict_;

        //- Solution Registry Reference
        const objectRegistry& obr_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Combstion model owner
        const combustionModel& combModel_;

        //- Fuel name
        word fuel_;

        //- Fuel consumption rate per unit of flame area
        volScalarField omega_;

        const fvMesh& getMesh(const objectRegistry& obr) const;

private:

    // Private member functions

         //- Disallow copy construct
        reactionRateFlameArea(const reactionRateFlameArea&);

         //- Disallow default bitwise assignment
        void operator=(const reactionRateFlameArea&);


public:

    //- Runtime type information
    TypeName("reactionRateFlameArea");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            reactionRateFlameArea,
            dictionary,
            (
                const word modelType,
                const dictionary& dict,
                const objectRegistry& obr,
                const combustionModel& combModel
            ),
            (modelType, dict, obr, combModel)
        );


    // Constructors

        //- Construct from components
        reactionRateFlameArea
        (
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr,
            const combustionModel& combModel
        );


    // Selector

        static autoPtr<reactionRateFlameArea> New
        (
            const dictionary& dict,
            const objectRegistry& obr,
            const combustionModel& combModel
        );


    // Destructor

        virtual ~reactionRateFlameArea();


    // Member functions

        //- Access functions

            //- Return omega
            const volScalarField& omega() const
            {
                return omega_;
            }


        //- Correct omega
        virtual void correct(const volScalarField& sigma) = 0;

        //- Update from dictionary
        virtual bool read(const dictionary& dictProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

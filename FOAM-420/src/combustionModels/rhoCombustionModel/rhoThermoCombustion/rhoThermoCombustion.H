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
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::rhoThermoCombustion

Description
    Density-based thermo model wrapper for combustion models

SourceFiles
    rhoThermoCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef rhoThermoCombustion_H
#define rhoThermoCombustion_H

#include "memory/autoPtr/autoPtr.H"
#include "rhoCombustionModel/rhoCombustionModel/rhoCombustionModel.H"
#include "chemistryModel/rhoChemistryModel/rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    class rhoThermoCombustion Declaration
\*---------------------------------------------------------------------------*/

class rhoThermoCombustion
:
    public rhoCombustionModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoThermoCombustion(const rhoThermoCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const rhoThermoCombustion&);


protected:

    // Protected data

        //- Pointer to thermo model
        rhoReactionThermo *thermoPtr_;


public:

    // Constructors

        //- Construct from components
        rhoThermoCombustion
        (
            const word& modelType,
            const objectRegistry& obr,
            const word& phaseName
        );


    //- Destructor
    virtual ~rhoThermoCombustion();


    // Member Functions

        //- Return access to the thermo package
        virtual rhoReactionThermo& thermo();

        //- Return const access to the thermo package
        virtual const rhoReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

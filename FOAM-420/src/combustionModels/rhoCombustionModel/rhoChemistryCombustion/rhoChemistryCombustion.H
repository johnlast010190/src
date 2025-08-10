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
    Foam::rhoChemistryCombustion

Description
    Density-based chemistry model wrapper for combustion models

SourceFiles
    rhoChemistryCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef rhoChemistryCombustion_H
#define rhoChemistryCombustion_H

#include "memory/autoPtr/autoPtr.H"
#include "rhoCombustionModel/rhoCombustionModel/rhoCombustionModel.H"
#include "chemistryModel/rhoChemistryModel/rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    class rhoChemistryCombustion Declaration
\*---------------------------------------------------------------------------*/

class rhoChemistryCombustion
:
    public rhoCombustionModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoChemistryCombustion(const rhoChemistryCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const rhoChemistryCombustion&);


protected:

    // Protected data

        //- Pointer to chemistry model
        autoPtr<rhoChemistryModel> chemistryPtr_;


public:

    // Constructors

        //- Construct from components and thermo
        rhoChemistryCombustion
        (
            const word& modelType,
            const objectRegistry& obr,
            const word& combustionProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~rhoChemistryCombustion();


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

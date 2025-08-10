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
    (c) 2013-2013 Esi Ltd.

Class
    Foam::convergenceTermination

Description
    Exit calculation when convergence criteria satsified.

SourceFiles
    convergenceTermination.C

\*---------------------------------------------------------------------------*/

#ifndef convergenceTermination_H
#define convergenceTermination_H

#include "db/objectRegistry/objectRegistry.H"
#include "db/dictionary/dictionary.H"
#include "containers/LinkedLists/user/FIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class convergenceTermination Declaration
\*---------------------------------------------------------------------------*/

class convergenceTermination
{
public:

    // Public data

        //- Enumeration defining the type of action
        enum actionType
        {
            noWriteNow,    /*!< stop immediately without writing data */
            writeNow,      /*!< write data and stop immediately */
            nextWrite      /*!< stop the next time data are written */
        };

private:

    // Private data

        const objectRegistry& obr_;

        //- Action type names
        static const NamedEnum<actionType, 3> actionTypeNames_;

        //- Name of field
        const word fieldName_;

        //- The type of action
        actionType action_;

        const scalar stdDevTol_;

        const label sampleSize_;

        FIFOStack<scalar> storedData_;

public:

    // Constructors

        //- Construct from description and refinemenet surfaces
        convergenceTermination
        (
            const objectRegistry& obr,
            const dictionary& dict
        );

      // Member Functions

        //- Read the dictionary settings
        void read(const dictionary&);

        // Access

            const word fieldName() const
            {
                return fieldName_;
            };

       // Query

            bool calculate(const scalar value);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::mixture

Description

SourceFiles
    mixtureI.H
    mixture.C
    mixtureIO.C

\*---------------------------------------------------------------------------*/

#ifndef mixture_H
#define mixture_H

#include "db/error/error.H"

#include "containers/Lists/List/List.H"
#include "substance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixture Declaration
\*---------------------------------------------------------------------------*/

class mixture
:
    public List<substance>
{
    // Private data

        word name_;


public:

    // Constructors

        //- Construct from Istream
        mixture(Istream& is)
        :
            List<substance>(is),
            name_(is)
        {
            scalar volTot = 0;

            for (label i = 0; i < size(); i++)
            {
                volTot += operator[](i).volFrac();
            }

            if (volTot > 1.001 || volTot < 0.999)
            {
                FatalErrorInFunction
                    << "Sum of volume fractions for Mixture " << name_
                    << " = " << volTot << endl
                    << "Should equal one."
                    << abort(FatalError);
            }
        }


    // Member Functions

        // Access

            const word& name() const
            {
                return name_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

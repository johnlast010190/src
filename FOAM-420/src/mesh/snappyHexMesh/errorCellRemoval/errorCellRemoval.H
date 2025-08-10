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
    (c) 2018-2020 Esi Ltd.

Class
    Foam::errorCellRemoval

Description
    Removal of cells with mesh quality errors

SourceFiles
    errorCellRemoval.C

\*---------------------------------------------------------------------------*/

#ifndef errorCellRemoval_H
#define errorCellRemoval_H

#include "meshRefinement/meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class errorCellRemoval Declaration
\*---------------------------------------------------------------------------*/

class errorCellRemoval
{
private:

    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- removal dictionary
        const dictionary& removalDict_;

public:

    enum errorType
    {
        NOERROR  = 0,  //!< No errors
        ORTHO = 1,    //!< Orthogonality error
        VOL = 2,     //!< Volume error
        PYR = 4,      //!< Pyramid error
        WEIGHTS = 8,      //!< weights error
        AR = 16,      //!< aspect ratio error
        SKEW = 32      //!< skewness error
    };

    //- Runtime type information
    TypeName("errorCellRemoval");

    //- Perform cell removal
    void remove(const bool updateIntersections = true);

    // Constructors

        //- Construct from components
        errorCellRemoval
        (
            meshRefinement& meshRefiner,
            const dictionary& dict
        );


    // Destructor

    virtual ~errorCellRemoval(){};

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

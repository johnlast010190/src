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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::radiation::cloudScatter

Description
    Cloud radiation scatter model

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef radiationCloudScatter_H
#define radiationCloudScatter_H

#include "submodels/scatterModel/scatterModel/scatterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                         Class cloudScatter Declaration
\*---------------------------------------------------------------------------*/

class cloudScatter
:
    public scatterModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Cloud names
        wordList cloudNames_;


public:

    //- Runtime type information
    TypeName("cloudScatter");


    // Constructors

        //- Construct from components
        cloudScatter
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~cloudScatter();


    // Member Operators

        // Access

            //- Return scatter coefficient
            tmp<volScalarField> sigmaEff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2015 OpenCFD Ltd.

Class
    Foam::radiation::constantTransmissivity

Group
    grpRadiationTransmissivitySubModels

Description
    Constant radiation transmissivity coefficient

SourceFiles
    constantTransmissivity.C

\*---------------------------------------------------------------------------*/

#ifndef radiationConstantTransmissivity_H
#define radiationConstantTransmissivity_H

#include "submodels/transmissivityModel/transmissivityModel/transmissivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                   Class constantTransmissivity Declaration
\*---------------------------------------------------------------------------*/

class constantTransmissivity
:
    public transmissivityModel
{

    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Transmissivity coefficient / []
        scalar tau_;


public:

    //- Runtime type information
    TypeName("constantTransmissivity");


    // Constructors

        //- Construct from components
        constantTransmissivity(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~constantTransmissivity();


    // Member Functions

        //- Return scatter coefficient
        tmp<volScalarField> tauEff(const label bandI = 0) const;


        //- Is Grey
        inline bool isGrey() const
        {
            return true;
        }


        //- Number of bands
        inline label nBands() const
        {
            return 1;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

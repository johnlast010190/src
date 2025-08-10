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
    (c) 2016 OpenFOAM Foundation

Class
    Foam::wallBoilingModels::linear

Description
    Linear wall heat flux partitioning model.

    Proposed threshold liquid fractions:
      - alphaLiquid1 0.1
      - alphaLiquid0 0.05

    \verbatim
        Ioilev, A., Samigulin, M., Ustinenko (2007).
        Advances in the modeling of cladding heat transfer
        and critical heat flux in boiling water reactor fuel assemblies.
        In Proc. 12th International Topical Meeting on
        Nuclear Reactor Thermal Hydraulics (NURETH-12),
        Pittsburgh, Pennsylvania, USA.
    \endverbatim

SourceFiles
    linear.C

\*---------------------------------------------------------------------------*/

#ifndef linear_H
#define linear_H

#include "derivedFvPatchFields/wallBoilingSubModels/partitioningModels/partitioningModel/partitioningModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{

/*---------------------------------------------------------------------------*\
                         Class linear Declaration
\*---------------------------------------------------------------------------*/

class linear
:
    public partitioningModel
{
    // Private data

        //- Model parameters, threshold liquid phase fractions
        scalar alphaLiquid1_;
        scalar alphaLiquid0_;


public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from a dictionary
        linear(const dictionary& dict);


    //- Destructor
    virtual ~linear();


    // Member Functions

        //- Calculate and return the wall heat-flux partioning
        virtual tmp<scalarField> fLiquid(const scalarField& alphaLiquid) const;

        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace partitioningModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

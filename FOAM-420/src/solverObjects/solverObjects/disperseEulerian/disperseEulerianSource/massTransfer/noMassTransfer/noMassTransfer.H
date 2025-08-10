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
    Foam::massTransferModels::noMassTransfer

Description
    Ranz-Marshall correlation for turbulent heat transfer from the surface of a
    sphere to the surrounding fluid.

SourceFiles
    noMassTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef noMassTransfer_H
#define noMassTransfer_H

#include "disperseEulerian/disperseEulerianSource/massTransfer/massTransferModel/massTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decoupledEulerian
{

class phase;

namespace massTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class noMassTransfer Declaration
\*---------------------------------------------------------------------------*/

class noMassTransfer
:
    public massTransferModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noMassTransfer
        (
            const dictionary& dict,
            const phase& phase
        );


    //- Destructor
    virtual ~noMassTransfer();


    // Member Functions

        //- Compute the mass transfer rate
        virtual tmp<volScalarField> mDot() const;

        //- Compute the Sherwood number
        virtual tmp<volScalarField> Sh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massTransferModels
} // End namespace decoupledEulerian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

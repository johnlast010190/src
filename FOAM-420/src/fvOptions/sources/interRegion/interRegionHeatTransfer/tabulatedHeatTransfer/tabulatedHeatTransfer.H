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
    Foam::fv::tabulatedHeatTransfer

Group
    grpFvOptionsSources

Description
    Tabulated heat transfer model. The heat exchange area per unit volume
    must be provided.  The 2D table returns the heat transfer coefficient
    by querying the local and neighbour region velocities

\*---------------------------------------------------------------------------*/

#ifndef tabulatedHeatTransfer_H
#define tabulatedHeatTransfer_H

#include "sources/interRegion/interRegionHeatTransfer/interRegionHeatTransferModel/interRegionHeatTransferModel.H"
#include "memory/autoPtr/autoPtr.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class tabulatedHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class tabulatedHeatTransfer
:
    public interRegionHeatTransferModel
{
private:

    // Private data

        //- Name of velocity field; default = U
        word UName_;

        //- Name of neighbour velocity field; default = U
        word UNbrName_;

        //- 2D look up table
        autoPtr<interpolation2DTable<scalar>> hTable_;

        //- Area per unit volume of heat exchanger
        autoPtr<volScalarField> AoV_;

        //- Heat transfer coefficient table
        const interpolation2DTable<scalar>& hTable();

        //- Field of area divided by volume
        const volScalarField& AoV();

        //- Start time name
        const word startTimeName_;


public:

    //- Runtime type information
        TypeName("tabulatedHeatTransfer");


    // Constructors

        //- Construct from dictionary
        tabulatedHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~tabulatedHeatTransfer();


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc(scalar deltaT);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

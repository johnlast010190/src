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
    (c) 2018 Esi Ltd.
    (c) 2011 OpenFOAM Foundation

Class
    Foam::extrudeModels::linearRadialCylindrical

Description

\*---------------------------------------------------------------------------*/

#ifndef linearRadialCylindrical_H
#define linearRadialCylindrical_H

#include "extrudeModel/extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                           Class linearRadialCylindrical Declaration
\*---------------------------------------------------------------------------*/

class linearRadialCylindrical
:
    public extrudeModel
{
    // Private data

        scalar R_;
        scalar Rsurface_;

        //- optional origin
        vector origin_;

        //- cylinder axis
        vector axis_;

public:

    //- Runtime type information
    TypeName("linearRadialCylindrical");

    // Constructors

        //- Construct from dictionary
        linearRadialCylindrical(const dictionary& dict);


    //- Destructor
    virtual ~linearRadialCylindrical();


    // Member Operators

        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

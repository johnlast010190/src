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
    Foam::SpecularReflection

Description
    Reversing the wall-normal component of the particle velocity

\*---------------------------------------------------------------------------*/

#ifndef SpecularReflection_H
#define SpecularReflection_H

#include "submodels/WallInteractionModel/WallInteractionModel/WallInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class SpecularReflection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SpecularReflection
:
    public WallInteractionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("SpecularReflection");


    // Constructors

        //- Construct from dictionary
        SpecularReflection
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~SpecularReflection();


    // Member Functions

        //- Apply wall correction
        virtual void correct
        (
            typename CloudType::parcelType& p,
            const wallPolyPatch& wpp
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/WallInteractionModel/SpecularReflection/SpecularReflection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

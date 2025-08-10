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
    Foam::Rebound

Group
    grpLagrangianIntermediatePatchInteractionSubModels

Description
    Simple rebound patch interaction model

\*---------------------------------------------------------------------------*/

#ifndef Rebound_H
#define Rebound_H

#include "submodels/Kinematic/PatchInteractionModel/PatchInteractionModel/PatchInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class Rebound Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Rebound
:
    public PatchInteractionModel<CloudType>
{
    // Private data

       //- Factor applied to velocity on rebound
       //  Normal rebound = 1
       scalar UFactor_;


public:

    //- Runtime type information
    TypeName("rebound");


    // Constructors

        //- Construct from dictionary
        Rebound(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        Rebound(const Rebound<CloudType>& pim);

        //- Construct and return a clone
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new Rebound<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Rebound();


    // Member Functions
        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle,
            const scalar trackFraction,
            const tetIndices& tetIs
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/Kinematic/PatchInteractionModel/Rebound/Rebound.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

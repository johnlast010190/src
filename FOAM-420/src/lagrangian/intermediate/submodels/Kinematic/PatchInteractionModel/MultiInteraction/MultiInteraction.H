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
    Foam::MultiInteraction

Group
    grpLagrangianIntermediatePatchInteractionSubModels

Description
    Runs multiple patch interaction models in turn. Takes dictionary
    where all the subdictionaries are the interaction models.

        // Exit upon first successful interaction or continue doing other
        // models. Returned nteraction status will be true if there has been any
        // interaction (so logical or)
        oneInteractionOnly true;

        model1
        {
            patchInteractionModel coincidentBaffleInteraction;
            coincidentBaffleInteractionCoeffs
            {
                coincidentPatches
                (
                    (pipetteWall_A pipetteCyclic_half0)
                    (pipetteWall_B pipetteCyclic_half1)
                );
            }
        }
        model2
        {
            patchInteractionModel localInteraction;
            localInteractionCoeffs
            {
                patches
                (
                    cWall
                    {
                        type rebound;
                    }
                    pipetteWall_A
                    {
                        type rebound;
                    }
                    pipetteWall_B
                    {
                        type rebound;
                    }
                );
            }
        }


\*---------------------------------------------------------------------------*/

#ifndef MultiInteraction_H
#define MultiInteraction_H

#include "submodels/Kinematic/PatchInteractionModel/PatchInteractionModel/PatchInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class MultiInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MultiInteraction
:
    public PatchInteractionModel<CloudType>
{
    // Private data

        Switch oneInteractionOnly_;

        //- Submodels
        PtrList<PatchInteractionModel<CloudType>> models_;


    // Private Member Functions

        //- Read settings
        bool read(const dictionary&);

public:

    //- Runtime type information
    TypeName("multiInteraction");


    // Constructors

        //- Construct from dictionary
        MultiInteraction(const dictionary& dict, CloudType& cloud);

        //- Construct copy from owner cloud and patch interaction model
        MultiInteraction(const MultiInteraction<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new MultiInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MultiInteraction();


    // Member Functions

        //- Flag to indicate whether model activates patch interaction model
        virtual bool active() const;

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
    #include "submodels/Kinematic/PatchInteractionModel/MultiInteraction/MultiInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

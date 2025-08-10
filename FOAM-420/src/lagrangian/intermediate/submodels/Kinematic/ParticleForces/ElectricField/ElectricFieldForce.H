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
    (c) 2010-2012 Esi Ltd.
    (c) 2011 OpenFOAM Foundation

Class
    Foam::ElectricFieldForce

Description
    Calculates particle Electric field force (qE), for particles with
    identical charges q in an electric field E.  Does not include
    particle-particle charge interaction.

Author:
    Dan Combest

SourceFiles
    ElectricFieldForceI.H
    ElectricFieldForce.C

\*---------------------------------------------------------------------------*/

#ifndef ElectricFieldForce_H
#define ElectricFieldForce_H

#include "submodels/Kinematic/ParticleForces/ParticleForce/ParticleForce.H"
#include "interpolation/interpolation/interpolation/interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class ElectricFieldForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ElectricFieldForce
:
    public ParticleForce<CloudType>
{
    // Private data

        //- Name of ElectricField field
        const word EFieldName_;

        //- EField interpolator
        const interpolation<vector>* EFieldInterpPtr_;

        //- Particle Charge
        const scalar particleCharge_;

public:

    //- Runtime type information
    TypeName("ElectricField");


    // Constructors

        //- Construct from mesh
        ElectricFieldForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        ElectricFieldForce(const ElectricFieldForce& gf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new ElectricFieldForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ElectricFieldForce();


    // Member Functions

        // Access

            //- Return the name of ElectricField field
            const word& EFieldName() const;

            //- Return the particle charge
            scalar particleCharge() const;

        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Kinematic/ParticleForces/ElectricField/ElectricFieldForceI.H"

#ifdef NoRepository
    #include "submodels/Kinematic/ParticleForces/ElectricField/ElectricFieldForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

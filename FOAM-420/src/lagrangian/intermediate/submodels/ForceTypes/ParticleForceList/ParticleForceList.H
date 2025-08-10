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
    Foam::ParticleForceList

Description
    List of particle forces

SourceFiles
    ParticleForceListI.H
    ParticleForceList.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleForceList_H
#define ParticleForceList_H

#include "submodels/Kinematic/ParticleForces/ParticleForce/ParticleForce.H"
#include "submodels/Kinematic/ParticleForces/forceSuSp/forceSuSp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ParticleForceList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleForceList
:
    public PtrList<ParticleForce<CloudType>>
{
    // Private data

        //- Reference to the owner cloud
        CloudType& owner_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Forces dictionary
        const dictionary dict_;

        //- Calculate coupled forces flag
        bool calcCoupled_;

        //- Calculate non-coupled forces flag
        bool calcNonCoupled_;


public:

    // Constructors

        //- Null constructor
        ParticleForceList(CloudType& owner, const fvMesh& mesh);

        //- Construct from mesh
        ParticleForceList
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool readFields
        );

        //- Construct copy
        ParticleForceList(const ParticleForceList& pfl);


    //- Destructor
    virtual ~ParticleForceList();


    // Member Functions

        // Access

            //- Return const access to the cloud owner
            inline const CloudType& owner() const;

            //- Return refernce to the cloud owner
            inline CloudType& owner();

            //- Return the mesh database
            inline const fvMesh& mesh() const;

            //- Return the forces dictionary
            inline const dictionary& dict() const;

            //- Set the calcCoupled flag
            inline void setCalcCoupled(bool flag);

            //- Set the calcNonCoupled flag
            inline void setCalcNonCoupled(bool flag);


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            //- Return the effective mass
            virtual scalar massEff
            (
                const typename CloudType::parcelType& p,
                const scalar mass
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/ForceTypes/ParticleForceList/ParticleForceListI.H"

#ifdef NoRepository
    #include "submodels/ForceTypes/ParticleForceList/ParticleForceList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::CloudFunctionObjectList

Description
    List of cloud function objects

SourceFiles
    CloudFunctionObjectListI.H
    CloudFunctionObjectList.C

\*---------------------------------------------------------------------------*/

#ifndef CloudFunctionObjectList_H
#define CloudFunctionObjectList_H

#include "containers/Lists/PtrList/PtrList.H"
#include "submodels/CloudFunctionObjects/CloudFunctionObject/CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class CloudFunctionObjectList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudFunctionObjectList
:
    public PtrList<CloudFunctionObject<CloudType>>
{
protected:

    // Protected Data

        //- Reference to the owner cloud
        const CloudType& owner_;

        //- Dictionary
        const dictionary dict_;


public:

    // Constructors

        //- Null constructor
        CloudFunctionObjectList(CloudType& owner);

        //- Construct from mesh
        CloudFunctionObjectList
        (
            CloudType& owner,
            const dictionary& dict,
            const bool readFields
        );

        //- Construct copy
        CloudFunctionObjectList(const CloudFunctionObjectList& ppml);


    //- Destructor
    virtual ~CloudFunctionObjectList();


    // Member Functions

        // Access

            //- Return const access to the cloud owner
            inline const CloudType& owner() const;

            //- Return refernce to the cloud owner
            inline CloudType& owner();

            //- Return the forces dictionary
            inline const dictionary& dict() const;


        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-evolve hook
            virtual void postEvolve();

            //- Post-move hook
            virtual void postMove
            (
                typename CloudType::parcelType& p,
                const label celli,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );

            //- Post-patch hook
            virtual void postPatch
            (
                const typename CloudType::parcelType& p,
                const polyPatch& pp,
                const scalar trackFraction,
                const tetIndices& tetIs,
                bool& keepParticle
            );

            //- Post-face hook
            virtual void postFace
            (
                const typename CloudType::parcelType& p,
                const label facei,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/CloudFunctionObjects/CloudFunctionObjectList/CloudFunctionObjectListI.H"

#ifdef NoRepository
    #include "submodels/CloudFunctionObjects/CloudFunctionObjectList/CloudFunctionObjectList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

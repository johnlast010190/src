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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::areaAveragedForce

Description
    Records particle face quantities on used-specified face zone

    Currently supports:
        accummulated mass
        average mass flux

SourceFiles
    areaAveragedForce.C

\*---------------------------------------------------------------------------*/

#ifndef areaAveragedForce_H
#define areaAveragedForce_H

#include "submodels/CloudFunctionObjects/CloudFunctionObject/CloudFunctionObject.H"
#include "meshes/polyMesh/zones/faceZone/faceZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class areaAveragedForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class areaAveragedForce
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Face zone IDs
        labelList faceZoneIDs_;

        //- Surface output format
        const word surfaceFormat_;

        //- tangential force magnitude
        List<vectorField> Fw_;

        //- field of total tangential force magnitude in zone
        List<vector> Fwt_;

        //- field of total tangential pressure magnitude in zone
        List<scalar> Pwt_;

        //- normal force magnitude
        List<vectorField> Fn_;

        //- field of total normal force magnitude in zone
        List<vector> Fnt_;

        //- field of total normal pressure magnitude in zone
        List<scalar> Pnt_;

        //- total area
        List<scalar> totalArea_;

        //- bulk direction used to determine tangential direction
        vector bulkDir_;

        //- Flag to indicate whether data should be written to file
        Switch log_;

        //- Output file pointer per zone
        PtrList<OFstream> outputFilePtr_;

        //- Last calculation time
        scalar timeOld_;


    // Private Member Functions

        //- Helper function to create log files
        void makeLogFile
        (
            const word& zoneName,
            const label zoneI,
            const label nFaces,
            const scalar totArea
        );


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("areaAveragedForce");


    // Constructors

        //- Construct from dictionary
        areaAveragedForce
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        areaAveragedForce(const areaAveragedForce<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new areaAveragedForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~areaAveragedForce();


    // Member Functions

        // Evaluation

            //- Post-face hook
            virtual void postFace
            (
                const parcelType& p,
                const label faceI,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "submodels/CloudFunctionObjects/areaAveragedForce/areaAveragedForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

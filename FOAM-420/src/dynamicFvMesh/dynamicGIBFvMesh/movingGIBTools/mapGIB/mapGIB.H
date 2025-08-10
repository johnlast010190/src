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
    (c) 2016 Esi Ltd.

Class
    Foam::deformingBodyMotionFunctions::mapGIB

Description
    SolidBodyMotionFvMesh 6DoF motion function. Constant velocity displacement.

SourceFiles
    mapGIB.C

\*---------------------------------------------------------------------------*/

#ifndef mapGIB_H
#define mapGIB_H

#include "db/typeInfo/typeInfo.H"
#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"
#include "interpolations/patchToPatchInterpolation/patchToPatchInterpolation.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "interpolations/patchToPatchInterpolation/PatchToPatchInterpolation.H"
#include "fields/Fields/scalarField/scalarIOField.H"
#include "AMIInterpolation/AMIInterpolation/AMIPatchToPatchInterpolation.H"
#include "primitives/Tuple2/Tuple2.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dynamicGIBFvMesh;

/*---------------------------------------------------------------------------*\
                          Class mapGIB Declaration
\*---------------------------------------------------------------------------*/

class mapGIB
{
    // Private data
    const dynamicGIBFvMesh& mesh_;

    bool mapSlave_;

    // AMI matching settings
    const scalar AMICosMatchAngle_;
    const scalar maxAMIWeightScale_;

    List<Tuple2< word, Field<scalar >>> gmfvScalars_;
    List<Tuple2< word, Field<vector>>> gmfvVectors_;
    List<Tuple2< word, Field<sphericalTensor>>> gmfvSpTensors_;
    List<Tuple2< word, Field<symmTensor>>> gmfvSymmTensors_;
    List<Tuple2< word, Field<tensor>>> gmfvTensors_;

    List<Tuple2< word, Field<scalar>>> gmfvsScalars_;
    List<Tuple2< word, Field<vector>>> gmfvsVectors_;
    List<Tuple2< word, Field<sphericalTensor>>> gmfvsSpTensors_;
    List<Tuple2< word, Field<symmTensor>>> gmfvsSymmTensors_;
    List<Tuple2< word, Field<tensor>>> gmfvsTensors_;

    List<Tuple2< word, Field<scalar>>> gsfvScalars_;
    List<Tuple2< word, Field<vector>>> gsfvVectors_;
    List<Tuple2< word, Field<sphericalTensor>>> gsfvSpTensors_;
    List<Tuple2< word, Field<symmTensor>>> gsfvSymmTensors_;
    List<Tuple2< word, Field<tensor>>> gsfvTensors_;

    List<Tuple2< word, Field<scalar>>> gsfvsScalars_;
    List<Tuple2< word, Field<vector>>> gsfvsVectors_;
    List<Tuple2< word, Field<sphericalTensor>>> gsfvsSpTensors_;
    List<Tuple2< word, Field<symmTensor>>> gsfvsSymmTensors_;
    List<Tuple2< word, Field<tensor>>> gsfvsTensors_;

    mutable AMIPatchToPatchInterpolation* mAMIInterPtr_;
    mutable AMIPatchToPatchInterpolation* sAMIInterPtr_;

    //-
    const primitivePatch& sourcePatch_;


    //- local patch data
    mutable primitivePatch* sourcePatchm_;
    mutable primitivePatch* sourcePatchs_;

    Field<vector> sourcePointsm_;
    Field<vector> sourcePointss_;
    faceList sourceFacesm_;
    faceList sourceFacess_;

    mutable primitivePatch* targetPatchm_;
    mutable primitivePatch* targetPatchs_;
    Field<vector> targetPointsm_;
    Field<vector> targetPointss_;
    faceList targetFacesm_;
    faceList targetFacess_;


    //- global patch data
    Field<vector> gppPoints_;
    faceList gppFaces_;

    // Private Member Functions

        void storeFields();
        void mapSerial();
        void combinePolyPatch();
        void mapMaster();
        void mapSlave();

        void mapNonOverlappingTargetFaces();

        labelList calcNearestPatchFaceMapping
        (
            const labelList&,
            const primitivePatch&,
            const primitivePatch&,
            const labelList&
        ) const;

        void makeMasterAMIPatchToPatchInterpolation() const;
        void makeSlaveAMIPatchToPatchInterpolation() const;

        //- Clear all demand-driven data
        virtual void clearOutGIBData();

        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh
        >
        void StoreOldFieldsToPatch
        (
            List< Tuple2<word, Field<Type>>>&,
            const label&
        );


        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh
        >
        void MapGIBField
        (
            List<Tuple2<word, Field<Type>>>&,
            const label&,
            const AMIPatchToPatchInterpolation&
        );

        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh
        >
        void Map1to1GIBField
        (
            const List<Tuple2<word, Field<Type>>>&,
            const List<Tuple2<word, Field<Type>>>&
        );


        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh
        >
        void MapNonOverlapFaces
        (
            const List<Tuple2<word, Field<Type>>>&,
            const label,
            const labelList&,
            const labelList&
        );


public:

    //- Runtime type information
    TypeName("mapGIB");


    // Constructors

        //- Construct from components
        mapGIB
        (
            const dynamicGIBFvMesh&,
            const primitivePatch&,
            bool mapSlave = true,
            // Default to the most permissive possible AMI settings in the GIB
            const scalar& AMICosMatchAngle = -1,
            const scalar& maxAMIWeightScale = GREAT
        );


    //- Destructor
    virtual ~mapGIB()
    {
        clearOutGIBData();
    };


    // Member Functions


        //- update

        virtual void mapBcs();

        virtual void map1To1();

        virtual void mapNonOverlapFaces
        (
            const labelList&,
            const labelList&,
            const label
        );

        virtual void mapInterfaceVelocity(scalarField*);

        virtual triSurface triS();

        const AMIPatchToPatchInterpolation& masterAMIInter() const
        {
            if (!mAMIInterPtr_)
            {
                makeMasterAMIPatchToPatchInterpolation();
            }
            return *mAMIInterPtr_;
        }

        const AMIPatchToPatchInterpolation& slaveAMIInter() const
        {
            if (!sAMIInterPtr_)
            {
                makeSlaveAMIPatchToPatchInterpolation();
            }
            return *sAMIInterPtr_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "dynamicGIBFvMesh/movingGIBTools/mapGIB/mapGIBTemplates.C"
#endif

#endif

// ************************************************************************* //

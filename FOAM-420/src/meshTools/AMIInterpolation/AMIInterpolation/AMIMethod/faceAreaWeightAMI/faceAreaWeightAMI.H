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
    (c) 2016,2020 OpenCFD Ltd.
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

Class
    Foam::faceAreaWeightAMI

Description
    Face area weighted Arbitrary Mesh Interface (AMI) method

SourceFiles
    faceAreaWeightAMI.C

\*---------------------------------------------------------------------------*/

#ifndef faceAreaWeightAMI_H
#define faceAreaWeightAMI_H

#include "AMIInterpolation/AMIInterpolation/AMIMethod/AMIMethod/AMIMethod.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class faceAreaWeightAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class faceAreaWeightAMI
:
    public AMIMethod<SourcePatch, TargetPatch>
{
private:

    // Private data

        //- Flag to restart uncovered source faces
        bool restartUncoveredSourceFace_;

        //- If cos(angle). If inner product between two faces that overlap is
        //  larger than cosAngle_  --> add it to the normal overlap list
        //  otherwise add it to the symmetry overlap
        const scalar cosAngle_;


protected:

    // Protected Member Functions

        //- Disallow default bitwise copy construct
        faceAreaWeightAMI(const faceAreaWeightAMI&);

        //- Disallow default bitwise assignment
        void operator=(const faceAreaWeightAMI&);

        // Marching front

            //- Calculate addressing and weights using temporary storage
            virtual void calcAddressing
            (
                List<DynamicList<label>>& srcAddress,
                List<DynamicList<scalar>>& srcWeights,
                List<scalar>& srcSymWeights,
                List<DynamicList<point>>& srcCentroids,
                List<DynamicList<label>>& tgtAddress,
                List<DynamicList<scalar>>& tgtWeights,
                List<scalar>& tgtSymWeights,
                label srcFacei,
                label tgtFacei
            );

            //- Determine overlap contributions for source face srcFacei
            virtual bool processSourceFace
            (
                const label srcFacei,
                const label tgtStartFacei,
                DynamicList<label>& nbrFaces,
                DynamicList<label>& visitedFaces,
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<scalar>& srcSymWght,
                List<DynamicList<point>>& srcCtr,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght,
                List<scalar>& tgtSymWght
            );

            //- Attempt to re-evaluate source faces that have not been included
            virtual void restartUncoveredSourceFace
            (
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<scalar>& srcSymWght,
                List<DynamicList<point>>& srcCtr,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght,
                List<scalar>& tgtSymWght
            );

            //- Set the source and target seed faces
            virtual void setNextFaces
            (
                label& startSeedI,
                label& srcFacei,
                label& tgtFacei,
                const boolList& mapFlag,
                labelList& seedFaces,
                const DynamicList<label>& visitedFaces,
                bool errorOnNotFound = true
            ) const;


        // Evaluation

            //- Area of intersection between source and target faces
            virtual Tuple2<scalar, scalar> interArea
            (
                const label srcFacei,
                const label tgtFacei,
                point&
            ) const;


public:

    //- Runtime type information
    TypeName("faceAreaWeightAMI");


    // Constructors

        //- Construct from components
        faceAreaWeightAMI
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const scalarField& srcMagSf,
            const scalarField& tgtMagSf,
            const faceAreaIntersect::triangulationMode& triMode,
            const scalar& cosMatchAngle,
            const bool reverseTarget = false,
            const bool requireMatch = true,
            const bool restartUncoveredSourceFace = true
        );


    //- Destructor
    virtual ~faceAreaWeightAMI();


    // Member Functions

            //- Flag to indicate that interpolation patches are conformal
            virtual bool conformal() const;

            virtual bool restartUncoveredFaces() const
            {
                return restartUncoveredSourceFace_;
            }

        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                scalarList& symSrcWeights,
                pointListList& srcCentroids,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                scalarList& symTgtWeights,
                label srcFacei = -1,
                label tgtFacei = -1
            );

            //- Normalise the weight. Can optionally subset addressing
            //  (e.g. for mapNearest)
            virtual void normaliseWeights
            (
                const bool verbose,
                AMIInterpolation<SourcePatch, TargetPatch>& inter
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AMIInterpolation/AMIInterpolation/AMIMethod/faceAreaWeightAMI/faceAreaWeightAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

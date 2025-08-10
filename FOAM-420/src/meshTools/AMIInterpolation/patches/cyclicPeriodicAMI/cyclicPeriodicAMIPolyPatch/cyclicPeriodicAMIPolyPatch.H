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
    (c) 2015 OpenCFD Ltd.
    (c) 2018-2020 OpenFOAM Foundation
    (c) 2021-2022 Esi Ltd.

Class
    Foam::cyclicPeriodicAMIPolyPatch

Description
    Cyclic patch for periodic Arbitrary Mesh Interface (AMI)

SourceFiles
    cyclicPeriodicAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPeriodicAMIPolyPatch_H
#define cyclicPeriodicAMIPolyPatch_H

#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class Type>
class keepIfTrueOp
{
    public:

        Tuple2<bool, Type> operator()
        (
            const Tuple2<bool, Type>& x,
            const Tuple2<bool, Type>& y
        ) const
        {
            if (x.first())
            {
                return x;
            }
            else if (y.first())
            {
                return y;
            }
            else
            {
                return Tuple2<bool, Type>(false, Type());
            }
        }
};


class OBJstream;

/*---------------------------------------------------------------------------*\
                    Class cyclicPeriodicAMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicPeriodicAMIPolyPatch
:
    public cyclicAMIPolyPatch
{
private:

    // Private data

        //- Periodic patch name
        mutable word periodicPatchName_;

        //- Periodic patch ID
        mutable label periodicPatchID_;

        //- Current number of transformations (+ve forward, -ve backward)
        mutable label nTransforms_;

        //- Number of sectors in a rotationally periodic geometry (optional)
        mutable label nSectors_;

        //- Maximum number of attempts to match the AMI geometry
        const label maxIter_;

        mutable List<transformer> srcAMITransforms_;

        mutable List<transformer> tgtAMITransforms_;


    // Private Member Functions

        void computeSector() const;

        //- Synchronise the periodic transformations
        void syncTransforms() const;

        transformer getTransform(const coupledPolyPatch&) const;

        //- Debug: write obj files of patch (collected on master)
        void writeOBJ(const primitivePatch& p, OBJstream& str) const;

        //- Reset the AMI interpolator
        virtual void resetAMI
        (
            const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod =
                AMIPatchToPatchInterpolation::imFaceAreaWeight
        ) const;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("cyclicPeriodicAMI");


    // Constructors

        //- Construct from (base coupled patch) components
        cyclicPeriodicAMIPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        cyclicPeriodicAMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicPeriodicAMIPolyPatch
        (
            const cyclicPeriodicAMIPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicPeriodicAMIPolyPatch
        (
            const cyclicPeriodicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct given the original patch and a map
        cyclicPeriodicAMIPolyPatch
        (
            const cyclicPeriodicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicPeriodicAMIPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<directPolyPatch>
            (
                new cyclicPeriodicAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName_
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<directPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<directPolyPatch>
            (
                new cyclicPeriodicAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~cyclicPeriodicAMIPolyPatch();


    // Member Functions

        virtual bool implicitNonCovered() const
        {
            return false;
        }

        //- Number of sectrors
        virtual label nSectors() const;

        //- Periodic patch ID
        virtual label periodicPatchID() const;

        //- Neighbour periodic patch ID
        virtual label neiPeriodicPatchID() const;

        //- patch owner or neighbour that has more sectors (less angle)
        virtual label smallerSectorPatchID() const;

        //- true on | false nb fewer blades (more sectors)
        virtual bool onNbSmallerSector() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

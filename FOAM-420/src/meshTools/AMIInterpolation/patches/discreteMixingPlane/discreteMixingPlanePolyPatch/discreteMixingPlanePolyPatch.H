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
    Foam::discreteMixingPlanePolyPatch

Description
    This is a generalised frozen rotor approach.
    The complete 360 from both sides and all the possible sector-based
    combinations are constructed.
    The "classic" mixingPlane finds the average on different circumference
    levels/bands. In this approach the sectors are rotated in a discrete way
    depending on the sectors. If you have for example a 10-6 rotor stator
    sectors, the rotor will be rotated 10 times and the stator 6 times.
    All the possible combinations will be 60.
    Thus, it does not depend on the sector initial position as the the frozen
    rotor  approach (frozen rotor in literature).
    It has a small dependency.
    It is called “discrete” because discrete type of averaging is applied
    and not the band-based averaging as the classic mixingPlane.

SourceFiles
    discreteMixingPlanePolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef discreteMixingPlanePolyPatch_H
#define discreteMixingPlanePolyPatch_H

#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OBJstream;

/*---------------------------------------------------------------------------*\
                    Class discreteMixingPlanePolyPatch Declaration
\*---------------------------------------------------------------------------*/

class discreteMixingPlanePolyPatch
:
    public cyclicAMIPolyPatch
{
private:

    // Private data

        //- Periodic patch name
        mutable word periodicPatchName_;

        //- Periodic patch ID
        mutable label periodicPatchID_;

        //- Number of sectors in a rotationally periodic geometry (optional)
        mutable label nSectors_;

        //- Transformations for each rotation
        mutable List<transformer> srcAMITransforms_;

        //- Transformations for each rotation
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
    TypeName("discreteMixingPlane");


    // Constructors

        //- Construct from (base coupled patch) components
        discreteMixingPlanePolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        discreteMixingPlanePolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        discreteMixingPlanePolyPatch
        (
            const discreteMixingPlanePolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        discreteMixingPlanePolyPatch
        (
            const discreteMixingPlanePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct given the original patch and a map
        discreteMixingPlanePolyPatch
        (
            const discreteMixingPlanePolyPatch& pp,
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
                new discreteMixingPlanePolyPatch(*this, bm)
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
                new discreteMixingPlanePolyPatch
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
                new discreteMixingPlanePolyPatch
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
    virtual ~discreteMixingPlanePolyPatch();


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

        virtual label onSectors() const;

        virtual label nbSectors() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2019 Esi Ltd.

Class
    Foam::mappedFvPatch

Description
    A region-coupled patch, which acts as an lduInterface for
    implicit coupling if used with a monolithicMatrix

SourceFiles
    mappedFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFvPatch_H
#define mappedFvPatch_H

#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyPatches/polyPatch/polyPatch.H"
#include "mappedPatches/mappedPolyPatch/mappedPolyPatch.H"
#include "mappedPatches/mappedLduInterface/mappedLduInterface.H"
#include "regionCoupled/patches/regionCoupledPolyPatch/regionCoupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mappedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedFvPatch
:
    virtual public fvPatch,
    public mappedLduInterface,
    public lduInterface
{
    // Private data

        //- Reference to mappedPolyPatch
        const mappedPatchBase& mappedPolyPatch_;


protected:

    // Protected members

        //- Returns fvMesh
        const fvMesh& nbrFvMesh() const
        {
            return refCast<const fvMesh>(nbrMesh());
        }

        virtual void initMovePoints()
        {
            // If this mesh moves but the other doesn't, the other
            // polyPatch does not clear its data as it needs to, e.g. to
            // recompute the AMI if it owns the AMI. Hence, force clear.
            const_cast<mappedPatchBase&>
            (
                dynamicCast<const mappedPatchBase>
                (
                    mappedPolyPatch_.samplePolyPatch()
                )
            ).clearOut();

            fvPatch::initMovePoints();
        }

public:

    //- Runtime type information
    TypeName(mappedPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        mappedFvPatch
        (
            const polyPatch& pp,
            const fvBoundaryMesh& bm
        )
        :
            fvPatch(pp, bm),
            mappedLduInterface(),
            mappedPolyPatch_
            (
                refCast<const mappedPatchBase>(pp)
            )
        {}


    //- Destructor
    virtual ~mappedFvPatch()
    {}


    // Member Functions

        // Access

            //- Whether it is able to act as an implicitly coupled interface
            virtual bool isInterface() const
            {
                return
                (
                    (
                        this->mappedPolyPatch_.mode()
                        == mappedPolyPatch::NEARESIPATCHFACE
                     || this->mappedPolyPatch_.mode()
                        == mappedPolyPatch::NEARESIPATCHFACEAMI
                    )
                 && nbrMeshFound()
                 && isA<mappedFvPatch>(nbrFvMesh().boundary()[nbrPatchID()])
                );
            }

            //- Return neighbour
            virtual label nbrPatchID() const
            {
                return mappedPolyPatch_.samplePolyPatch().index();
            }

            //- Return mappedFvPatch neighbour patch
            virtual const mappedFvPatch& nbrPatch() const
            {
                return refCast<const mappedFvPatch>
                (
                    nbrFvMesh().boundary()[nbrPatchID()]
                );
            }

            //- Does neighbour mesh exist?
            virtual bool nbrMeshFound() const
            {
                return mappedPolyPatch_.sampleMeshFound();
            }

            //- Returns neighbour polyMesh
            virtual const polyMesh& nbrMesh() const
            {
                return mappedPolyPatch_.sampleMesh();
            }

            //- Return fvPatch
            const fvPatch& patch() const
            {
                return *this;
            }

            //- Is it the owner of the AMI?
            virtual bool owner() const
            {
                if (mappedPolyPatch_.sameRegion())
                {
                    return this->index() < nbrPatchID();
                }
                return
                (
                    this->boundaryMesh().mesh().name()
                  < this->nbrFvMesh().name()
                );
            }

            //- Return a reference to the AMI interpolator
            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                return mappedPolyPatch_.AMI();
            }

            //- Return reference to direct map
            virtual const mapDistribute& map() const
            {
                return mappedPolyPatch_.map();
            }

            //- Whether we are using AMI or direct map
            virtual bool usingAMI() const
            {
                return
                (
                    mappedPolyPatch_.mode()
                 == mappedPolyPatch::NEARESIPATCHFACEAMI
                );
            }


            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return fvPatch::faceCells();
            }

            //- Returns false - no interpolation across patch
            virtual bool coupled() const
            {
                return false;
            }

            virtual const word& patchName() const
            {
                return this->name();
            }

            virtual const word& regionName() const
            {
                return this->boundaryMesh().mesh().name();
            }

            virtual const word& nbrPatchName() const
            {
                return nbrPatch().name();
            }

            virtual const word& nbrRegionName() const
            {
                return nbrPatch().boundaryMesh().mesh().name();
            }

            //- Return the interface type
            const word& regionCoupleType() const
            {
                return patch().type();
            }

            //- Name of the sample/neighbour region
            const word& sampleRegion() const
            {
                return mappedPolyPatch_.sampleRegion();
            }

            //- Map or use the AMI to perform interpolation from neighbour
            // defaultValues are used for AMI low-weight correction
            template <class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues
            ) const;

            //- Map or use the AMI to perform interpolation from neighbour
            // defaultValues are used for AMI low-weight correction
            template <class Type>
            tmp<Field<Type>> interpolate
            (
                tmp<Field<Type>> tfld,
                const UList<Type>& defaultValues
            ) const;

    // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                labelUList& iF
            ) const
            {}

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh/fvPatches/derived/mapped/mappedFvPatchTemplates.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

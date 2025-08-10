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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::nonConformalErrorPolyPatch

Description
    Non-conformal error poly patch. As nonConformalPolyPatch. This patch is a
    non-coupled non-conformal patch which is used to manage the errors created
    during the construction of a non-conformal coupled interface. Every patch
    used as the original patch of the owner side of a non-conformal coupled
    interface must also have an associated error patch.

See also
    Foam::nonConformalPolyPatch

SourceFiles
    nonConformalErrorPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalErrorPolyPatch_H
#define nonConformalErrorPolyPatch_H

#include "meshes/polyMesh/polyPatches/directPolyPatch/directPolyPatch.H"
#include "nonConformal/polyPatches/nonConformal/nonConformalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalCyclicPolyPatch;

/*---------------------------------------------------------------------------*\
                  Class nonConformalErrorPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalErrorPolyPatch
:
    public directPolyPatch,
    public nonConformalPolyPatch
{
protected:

    // Protected Member Functions

        //- Reset the patch name
        virtual void rename(const wordList& newNames) override;

        //- Reset the patch index
        virtual void reorder(const labelUList& oldToNewIndex) override;


public:

    //- Runtime type information
    TypeName("nonConformalError");


    // Constructors

        //- Construct from components
        nonConformalErrorPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from components
        nonConformalErrorPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& origPatchName
        );

        //- Construct from dictionary
        nonConformalErrorPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        nonConformalErrorPolyPatch
        (
            const nonConformalErrorPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        nonConformalErrorPolyPatch
        (
            const nonConformalErrorPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& origPatchName
        );

        //- Construct given the original patch and a map
        nonConformalErrorPolyPatch
        (
            const nonConformalErrorPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm
        ) const override
        {
            return autoPtr<polyPatch>
            (
                new nonConformalErrorPolyPatch(*this, bm)
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
        ) const override
        {
            return autoPtr<directPolyPatch>
            (
                new nonConformalErrorPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    origPatchName()
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
        ) const override
        {
            return autoPtr<directPolyPatch>
            (
                new nonConformalErrorPolyPatch
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
    virtual ~nonConformalErrorPolyPatch();


    // Member Functions

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

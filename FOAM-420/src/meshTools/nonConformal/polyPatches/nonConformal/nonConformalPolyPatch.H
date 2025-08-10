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
    Foam::nonConformalPolyPatch

Description
    Non-conformal poly patch. This patch is a placeholder and must have no
    faces. This patch is linked to an "original" patch. Other parts of the code
    (i.e., finite volume) will use this patch to create non-empty patches that
    do not conform to the poly mesh.

SourceFiles
    nonConformalPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalPolyPatch_H
#define nonConformalPolyPatch_H

#include "meshes/polyMesh/polyPatches/constraint/cyclic/cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class nonConformalPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalPolyPatch
{
    // Private data

        //- Reference to the polyPatch
        const polyPatch& patch_;

        //- Name of the patch from which this patch originates
        mutable word origPatchName_;

        //- Index of the patch from which this patch originates
        mutable label origPatchID_;


    // Private Member Functions

        //- Validate the size of the patch
        void validateSize() const;


protected:

    // Protected Member Functions

        //- Reset the patch name
        virtual void rename(const wordList& newNames);

        //- Reset the patch index
        virtual void reorder(const labelUList& oldToNewIndex);


public:

    //- Runtime type information
    TypeName("nonConformalBase");


    // Constructors

        //- Construct from a patch
        nonConformalPolyPatch(const polyPatch& patch);

        //- Construct from a patch and an original patch name
        nonConformalPolyPatch
        (
            const polyPatch& patch,
            const word& origPatchName
        );

        //- Construct from a patch and a dictionary
        nonConformalPolyPatch(const polyPatch& patch, const dictionary& dict);

        //- Construct from a non-conformal patch
        nonConformalPolyPatch(const nonConformalPolyPatch& patch);


    //- Destructor
    virtual ~nonConformalPolyPatch();


    // Member Functions

        //- Original patch name
        const word& origPatchName() const
        {
            return origPatchName_;
        }

        //- Original patch
        const polyPatch& origPatch() const
        {
            return patch_.boundaryMesh()[origPatchID()];
        }

        //- Original patchID
        label origPatchID() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

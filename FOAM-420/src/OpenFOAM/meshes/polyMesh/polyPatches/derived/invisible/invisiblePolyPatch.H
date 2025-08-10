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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

Class
    Foam::invisiblePolyPatch

Description
    Foam::invisiblePolyPatch

SourceFiles
    invisiblePolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef invisiblePolyPatch_H
#define invisiblePolyPatch_H

#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class invisiblePolyPatch Declaration
\*---------------------------------------------------------------------------*/

class invisiblePolyPatch
:
    public indirectPolyPatch
{
public:

    //- Runtime type information
    TypeName("invisible");


    // Constructors

        //- Construct from components
        invisiblePolyPatch
        (
            const word&,
            const label,
            const label,
            const word,
            const label,
            const polyBoundaryMesh&,
            const word&
        );

        //- Construct from dictionary
        invisiblePolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        invisiblePolyPatch(const invisiblePolyPatch&, const polyBoundaryMesh&);

        //- Construct as copy, resetting the boundary mesh
        invisiblePolyPatch
        (
            const invisiblePolyPatch&,
            const polyBoundaryMesh&,
            const label&
        );

        //- Construct as copy
        invisiblePolyPatch(const invisiblePolyPatch&);

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new invisiblePolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting thefaceZone
        //  and boundary mesh
        virtual autoPtr<indirectPolyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label& fzI
        ) const
        {
            return autoPtr<indirectPolyPatch>
            (
                new invisiblePolyPatch(*this, bm, fzI)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

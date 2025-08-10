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
    (c) 2021 ESI Ltd.

Class
    Foam::parhipDecomp

\*---------------------------------------------------------------------------*/

#ifndef parhipDecomp_H
#define parhipDecomp_H

#include "decompositionMethod/decompositionMethod.H"

#include "parhip_interface.h"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class parhipDecomp Declaration
\*---------------------------------------------------------------------------*/

class parhipDecomp
:
    public decompositionMethod
{
    typedef unsigned long long idxtype;

    // Private Member Functions

        //- Dump graph for debug
        void writeGraphSimple
        (
            const fileName& meshPath,
            idxtype* vtxdist,
            idxtype* xadj,
            idxtype* adjncy,
            MPI_Comm* comm
        ) const;

        //- Decompose. Handles size 0 arrays
        label decompose
        (
            const fileName& meshPath,
            const List<label>& adjncy,
            const List<label>& xadj,
            const scalarField& cWeights,
            std::vector<idxtype>& finalDecomp
        ) const;

        //- Low level decompose
        label decompose
        (
            const fileName& meshPath,
            const label adjncySize,
            const label adjncy[],
            const label xadjSize,
            const label xadj[],
            const scalarField& cWeights,
            std::vector<idxtype>& finalDecomp
        ) const;

        //- Disallow default bitwise copy construct and assignment
        void operator=(const parhipDecomp&);
        parhipDecomp(const parhipDecomp&);


public:

    //- The predefined KaHIP configuration types
    enum configs
    {
        FAST, //!< default
        ECO,
        STRONG,
        FASTSOCIAL,
        ECOSOCIAL,
        STRONGSOCIAL
    };


    //- The selection names for predefined KaHIP configurations
    static const NamedEnum<configs, 6> configNames;


    //- Runtime type information
    TypeName("parhip");


    // Constructors

        //- Construct given the decomposition dictionary and mesh
        parhipDecomp(const dictionary& decompositionDict);


    //- Destructor
    virtual ~parhipDecomp()
    {}


    // Member Functions

        virtual bool parallelAware() const
        {
            // parhip does not know about proc boundaries
            return true;
        }

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed). See note on weights in scotchDecomp.H
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Gets
        //  passed agglomeration map (from fine to coarse cells) and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively. See note on weights in scotchDecomp.H
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const labelList& agglom,
            const pointField& regionPoints,
            const scalarField& regionWeights
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided mesh connectivity.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        //  See note on weights in scotchDecomp.H
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

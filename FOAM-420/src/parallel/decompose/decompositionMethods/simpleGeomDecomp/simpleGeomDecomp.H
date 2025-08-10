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

Class
    Foam::simpleGeomDecomp

Description

SourceFiles
    simpleGeomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef simpleGeomDecomp_H
#define simpleGeomDecomp_H

#include "geomDecomp/geomDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleGeomDecomp Declaration
\*---------------------------------------------------------------------------*/

class simpleGeomDecomp
:
    public geomDecomp
{
    // Private Member Functions

        void assignToProcessorGroup(labelList&, const label) const;

        void assignToProcessorGroup
        (
            labelList& processorGroup,
            const label nProcGroup,
            const labelList& indices,
            const scalarField& weights,
            const scalar summedWeights
        ) const;

        labelList decomposeOneProc(const pointField& points) const;

        labelList decomposeOneProc
        (
            const pointField& points,
            const scalarField& weights
        ) const;

        //- Disallow default bitwise copy construct and assignment
        void operator=(const simpleGeomDecomp&);
        simpleGeomDecomp(const simpleGeomDecomp&);


public:

    //- Runtime type information
    TypeName("simple");


    // Constructors

        //- Construct given the decomposition dictionary
        simpleGeomDecomp(const dictionary& decompositionDict);


    //- Destructor
    virtual ~simpleGeomDecomp()
    {}


    // Member Functions

        virtual bool parallelAware() const
        {
            // simpleDecomp sends all points to the master which does
            // the decomposition.
            return true;
        }

        virtual labelList decompose(const pointField&);

        virtual labelList decompose(const pointField&, const scalarField&);

        virtual labelList decompose(const polyMesh&, const pointField& points)
        {
            return decompose(points);
        }

        virtual labelList decompose
        (
            const polyMesh&,
            const pointField& points,
            const scalarField& pointWeights
        )
        {
            return decompose(points, pointWeights);
        }

        //- Explicitly provided connectivity
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        )
        {
            return decompose(cc, cWeights);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2014-2016 OpenFOAM Foundation
    (c) 2019-2021 OpenCFD Ltd

Class
    Foam::blockEdges::BSplineEdge

Description
    A blockEdge interface for B-splines.

SourceFiles
    BSplineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_BSplineEdge_H
#define blockEdges_BSplineEdge_H

#include "blockEdges/blockEdge/blockEdge.H"
#include "blockEdges/BSplineEdge/BSpline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                      Class BSplineEdge Declaration
\*---------------------------------------------------------------------------*/

class BSplineEdge
:
    public blockEdge,
    public BSpline
{
    // Private Member Functions

        //- No copy construct
        BSplineEdge(const BSplineEdge&) = delete;

        //- No copy assignment
        void operator=(const BSplineEdge&) = delete;


public:

    //- Runtime type information
    TypeName("BSpline");


    // Constructors

        //- Construct from components
        BSplineEdge
        (
            const pointField& points,   //!< Referenced point field
            const edge& fromTo,         //!< Start/end in point field
            const pointField& internalPoints
        );

        //- Construct from components
        BSplineEdge
        (
            const pointField& points,   //!< Referenced point field
            const label from,           //!< Start point in point field
            const label to,             //!< End point in point field
            const pointField& internalPoints
        );

        //- Construct from Istream and point field.
        BSplineEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces&  /*unused*/,
            const pointField& points,  //!< Referenced point field
            Istream& is
        );


    //- Destructor
    virtual ~BSplineEdge() = default;


    // Member Functions

        //- The point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar) const;

        //- The length of the spline curve
        //  \note NotImplemented
        virtual scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

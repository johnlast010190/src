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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    lineRefinement

Description
    Refine objects intersected by a line given by two points

SourceFiles
    lineRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef lineRefinement_H
#define lineRefinement_H

#include "utilities/octrees/meshOctree/refinementControls/objectRefinement/objectRefinement.H"
#include "meshes/primitiveShapes/point/point.H"
#include "db/typeInfo/typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class lineRefinement Declaration
\*---------------------------------------------------------------------------*/

class lineRefinement
:
    public objectRefinement
{
    // Private data
        //- start point of the line
        point p0_;

        //- end point of the line
        point p1_;

public:

    //- Runtime type information
    TypeName("line");


    // Constructors

        //- Null construct
        lineRefinement();

        //- Construct from name, cell size, start point and end point
        lineRefinement
        (
            const word& name,
            const scalar cellSize,
            const direction additionalRefLevels,
            const point& p0,
            const point& p1
        );

        //- Construct from dictionary
        lineRefinement(const word& name, const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<objectRefinement> clone
        (
            const lineRefinement& sr
        ) const
        {
            return autoPtr<objectRefinement>
            (
                new lineRefinement
                (
                    sr.name(),
                    sr.cellSize(),
                    sr.additionalRefinementLevels(),
                    sr.p0_,
                    sr.p1_
                )
            );
        }

    // Member Functions

        //- check if the line intersects the box
        bool intersectsObject(const boundBox&) const;

    //- Return as dictionary of entries
        dictionary dict(bool ignoreType = false) const;

      // Write

        //- Write
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&, bool subDict = true) const;

    // Member Operators

        //- assign from dictionary
        void operator=(const dictionary&);

    // IOstream Operators

        Ostream& operator<<(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

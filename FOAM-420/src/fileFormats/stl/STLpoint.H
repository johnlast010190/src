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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016-2017 OpenCFD Ltd.

Class
    Foam::STLpoint

Description
    A vertex point or facet normal representation for STL files.

\*---------------------------------------------------------------------------*/

#ifndef STLpoint_H
#define STLpoint_H

#include "meshes/primitiveShapes/point/point.H"
#include "primitives/Vector/floatVector/floatVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class STLpoint Declaration
\*---------------------------------------------------------------------------*/

class STLpoint
:
    public floatVector
{

public:

    // Constructors

        //- Construct null
        inline STLpoint()
        {}

        //- Construct from single-precision point
        inline STLpoint(const Vector<float>& p)
        :
            floatVector(p)
        {}

        //- Construct from double-precision point
        inline STLpoint(const Vector<double>& p)
        :
            floatVector(float(p.x()), float(p.y()), float(p.z()))
        {}

        //- Construct from single-precision components
        inline STLpoint(float x, float y, float z)
        :
            floatVector(x, y, z)
        {}

        //- Construct from double-precision components
        inline STLpoint(double x, double y, double z)
        :
            floatVector(float(x), float(y), float(z))
        {}


    // Member Operators

        #ifdef FOAM_DP
        //- Conversion to double-precision point
        inline operator point() const
        {
            return point(x(), y(), z());
        }
        #endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


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
    (c) 2013-2022 Esi Ltd.

Class
    Foam::fieldInstance

Description
    Exit calculation when convergence criteria satsified.

SourceFiles
    fieldInstance.C

\*---------------------------------------------------------------------------*/

#ifndef fieldInstance_H
#define fieldInstance_H

#pragma once

#include "meshes/primitiveShapes/point/point.H"
#include "primitives/Tensor/tensor/tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class fieldInstance;
    Istream& operator>> (Istream& is, fieldInstance&);
    Ostream& operator<< (Ostream& os, const fieldInstance&);

/*---------------------------------------------------------------------------*\
                           Class fieldInstance Declaration
\*---------------------------------------------------------------------------*/

class fieldInstance
{
    // Private data

        scalar val;
        point pos;

public:

        // Construct from components
        fieldInstance(scalar v=0, point p=vector::zero);

        fieldInstance(const fieldInstance&) = default;

        // Access functions
        scalar value() const {return val;};
        point position() const {return pos;};

        // Manipulation functions
        void setValue(scalar v){val=v;}
        void setPosition(point p){pos=p;}

         // Member Operators

         //- Copy assignment
         fieldInstance& operator=(const fieldInstance&) = default;

         bool operator==(const fieldInstance& other) const {
             return val == other.val && pos == other.pos;
         }
         bool operator!=(const fieldInstance& other) const {
             return !(*this == other);
         }
};

template<>
constexpr bool contiguous<fieldInstance>() {return true;}

class maxFIOp
{
public:
    fieldInstance operator()
    (
        const fieldInstance& fI1,
        const fieldInstance& fI2
    ) const
    {
        fieldInstance fI(fI2);
        if (fI1.value() > fI2.value())
        {
            fI = fI1;
        }
        return fI;
    }
};

class minFIOp
{
public:
    fieldInstance operator()
    (
        const fieldInstance& fI1,
        const fieldInstance& fI2
    ) const
    {
        fieldInstance fI(fI2);
        if (fI1.value() < fI2.value())
        {
            fI = fI1;
        }
        return fI;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

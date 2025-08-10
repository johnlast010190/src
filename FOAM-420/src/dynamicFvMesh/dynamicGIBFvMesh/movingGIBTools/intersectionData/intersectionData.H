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
    (c) 2017 Esi Ltd.

Class

Description
    Holds information for points that are exchanged between processors.
    It is used for points leaving the domain in parallel.

SourceFiles
    intersectionData.C

\*---------------------------------------------------------------------------*/

#ifndef intersectionData_H
#define intersectionData_H

#include "db/typeInfo/typeInfo.H"
#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class intersectionData;
Istream& operator>>(Istream&, intersectionData&);
Ostream& operator<<(Ostream&, const intersectionData&);

/*---------------------------------------------------------------------------*\
                          Class intersectionData Declaration
\*---------------------------------------------------------------------------*/

class intersectionData final
{
    // Private data

    // patch point
    label pp_;

    // origin patch label
    label pI_;

    // local patch addressing (global-p.start())
    label pAdd_;

    // from processor
    label fromProc_;

    // to processor
    label toProc_;

    // label of the patch in the neighbour processor
    label neiProcPatch_;

    // starting point
    point pS_;

    // ending point
    point pE_;

    // wall hit
    bool hitWall_;

    // processor hit
    bool hitProc_;

    // processor hit
    bool hitNone_;

    // hitPoint
    point hitP_;


    // Private Member Functions


public:


    // Constructors

        //- Construct from components
        intersectionData
        (
            const label&,
            const label&,
            const label&,
            const label&,
            const label&,
            const label&,
            const point&,
            const point&,
            const bool&,
            const bool&,
            const bool&,
            const point&
        );
        intersectionData();

    // Member Functions


        //- update

        void checkParallelBoundaryIntersections();

        //- Access

        inline label pp() const
        {
            return pp_;
        }

        inline label pI() const
        {
            return pI_;
        }

        inline label pAdd() const
        {
            return pAdd_;
        }

        inline label fromProc() const
        {
            return fromProc_;
        }

        inline label toProc() const
        {
            return toProc_;
        }

        inline label neiProcPatch() const
        {
            return neiProcPatch_;
        }

        inline point pS() const
        {
            return pS_;
        }

        inline point pE() const
        {
            return pE_;
        }

        inline bool hitWall() const
        {
            return hitWall_;
        }

        inline bool hitProc() const
        {
            return hitProc_;
        }

        inline bool hitNone() const
        {
            return hitNone_;
        }

        inline point hitP() const
        {
            return hitP_;
        }

        inline label& pp()
        {
            return pp_;
        }

        inline label& pI()
        {
            return pI_;
        }

        inline label& pAdd()
        {
            return pAdd_;
        }

        inline label& fromProc()
        {
            return fromProc_;
        }

        inline label& toProc()
        {
            return toProc_;
        }

        inline label& neiProcPatch()
        {
            return neiProcPatch_;
        }

        inline point& pS()
        {
            return pS_;
        }

        inline point& pE()
        {
            return pE_;
        }

        inline bool& hitWall()
        {
            return hitWall_;
        }

        inline bool& hitProc()
        {
            return hitProc_;
        }

        inline bool& hitNone()
        {
            return hitNone_;
        }

        inline point& hitP()
        {
            return hitP_;
        }


    // Member Operators

         bool operator==(const intersectionData&) const;
         bool operator!=(const intersectionData&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, intersectionData&);
        friend Ostream& operator<<(Ostream&, const intersectionData&);
};

template<>
constexpr bool contiguous<intersectionData>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //

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
    (c) 2022 Esi Ltd.

Class
    Foam::ignitionSite

Description
    Foam::ignitionSite

SourceFiles
    ignitionSiteI.H
    ignitionSite.C
    ignitionSiteIO.C

\*---------------------------------------------------------------------------*/

#ifndef ignitionSite_H
#define ignitionSite_H

#include "primitives/Vector/vector/vector.H"
#include "primitives/ints/lists/labelList.H"
#include "primitives/Scalar/lists/scalarList.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;
class engineTime;
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class ignitionSite Declaration
\*---------------------------------------------------------------------------*/

class ignitionSite
{
    // Private data

        const Time& db_;
        const fvMesh& mesh_;

        dictionary ignitionSiteDict_;

        vector location_;
        scalar diameter_;
        scalar time_;
        scalar duration_;
        scalar strength_;

        labelList cells_;
        scalarList cellVolumes_;

        //- Current time index.
        //  Used during the update for moving meshes
        mutable label timeIndex_;


    // Private Member Functions

        void findIgnitionCells(const fvMesh&);


public:

    // Public classes

        //- Class used for the read-construction of
        //  PtrLists of ignitionSite
        class iNew
        {
            const Time& db_;
            const fvMesh& mesh_;

        public:

            iNew(const Time& db, const fvMesh& mesh)
            :
                db_(db),
                mesh_(mesh)
            {}

            autoPtr<ignitionSite> operator()(Istream& is) const
            {
                return autoPtr<ignitionSite>(new ignitionSite(is, db_, mesh_));
            }
        };


    // Constructors

        //- Copy construct
        ignitionSite(const ignitionSite&) = default;

        //- Construct from Istream and database
        ignitionSite(Istream&, const Time&, const fvMesh&);

        //- Construct from Istream and engineTime
        ignitionSite(Istream&, const engineTime&, const fvMesh&);

        //- Clone
        autoPtr<ignitionSite> clone() const
        {
            return autoPtr<ignitionSite>(new ignitionSite(*this));
        }


    // Member Functions

        // Access

            const vector& location() const
            {
                return location_;
            }

            scalar diameter() const
            {
                return diameter_;
            }

            scalar time() const
            {
                return time_;
            }

            scalar duration() const
            {
                return duration_;
            }

            scalar strength() const
            {
                return strength_;
            }

            //- Return the ignition cells updated if the mesh moved
            const labelList& cells() const;

            const scalarList& cellVolumes() const
            {
                return cellVolumes_;
            }


        // Check

            bool igniting() const;

            bool ignited() const;


    // Member Operators

        void operator=(const ignitionSite&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2018 Esi Ltd.

Class
    Foam::fieldProfile

Description
    Field profile container

SourceFiles
    fieldProfile.C

\*---------------------------------------------------------------------------*/

#ifndef fieldProfile_H
#define fieldProfile_H

#include "db/typeInfo/typeInfo.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "cfdTools/general/patchDistanceFunction/patchDistanceFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class fieldProfile Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fieldProfile
:
    public tmp<fieldProfile<Type>>::refCount
{
    // Private Data

        //- Patch
        const fvPatch& patch_;

        //- Dictionary
        dictionary dict_;

        //- The distribution being used, including the bounding treatment
        autoPtr<Function1<Type>> distribution_;

        //- Distance field
        autoPtr<patchDistanceFunction> distance_;

        //- Scaling value for input X data
        autoPtr<Function1<scalar>> xscale_;

        //- Scaling value for table Y data
        autoPtr<Function1<scalar>> yscale_;

        //- Offset value for input X data
        scalar xoffset_;

        //- Offset value for table Y data
        autoPtr<Function1<Type>> yoffset_;


    // Private functions

        //- Disallow default bitwise assignment
        void operator=(const fieldProfile&);


public:

    // Constructors

        //- Null constructor
        fieldProfile(const fvPatch&);

        //- Construct from components
        fieldProfile(const fvPatch&, const dictionary&);

        //- copy constructor
        fieldProfile(const fieldProfile&);

        //- Construct and return a clone
        tmp<fieldProfile> clone() const
        {
            return tmp<fieldProfile>
            (
                new fieldProfile(*this)
            );
        }


    //- Destructor
    virtual ~fieldProfile()
    {}


    // Member Functions

        //- update distance class
        bool update();

        //- Distance field
        inline const scalarField& Z() const
        {
            return distance_->value();
        }

        //- Axis
        inline const vector& axis() const
        {
            return distance_->axis();
        }

        //- Origin
        inline const vector& origin() const
        {
            return distance_->origin();
        }

        //- Field profile
        tmp<Field<Type>> value() const;

        //- Field gradient
        tmp<Field<Type>> gradient() const;

        //- Write
        virtual void write(Ostream& os, bool writeFrame = false) const;
};




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "cfdTools/general/fieldProfile/fieldProfile.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

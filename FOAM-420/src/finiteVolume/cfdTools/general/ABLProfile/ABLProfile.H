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
    Foam::ABLProfile

Description
    Abstract base class ABL profiles

SourceFiles
    ABLProfile.C

\*---------------------------------------------------------------------------*/

#ifndef ABLProfile_H
#define ABLProfile_H

#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "cfdTools/general/patchDistanceFunction/patchDistanceFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ABLProfile Declaration
\*---------------------------------------------------------------------------*/

class ABLProfile
{

protected:

        //- Patch
        const fvPatch& patch_;

        //- Distance field
        autoPtr<patchDistanceFunction> Z_;


private:

    // Private functions

        //- Disallow default bitwise assignment
        void operator=(const ABLProfile&);


public:

    //- Runtime type information
        TypeName("ABLProfile");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ABLProfile,
            dictionary,
            (
                const fvPatch& patch,
                const dictionary& ABLProfileDict
            ),
            (patch, ABLProfileDict)
        );

    // Constructors

        //- Null constructor
        ABLProfile(const fvPatch&);

        //- Construct from components
        ABLProfile(const fvPatch&, const dictionary&);

        //- copy constructor
        ABLProfile(const ABLProfile&);

        //- Construct and return a clone
        virtual ABLProfile* clone() const = 0;


    // Selectors

        //- Return an autoPtr to the selected ABLProfile
        static autoPtr<ABLProfile> New
        (
            const fvPatch& patch,
            const dictionary& ABLProfileDict
        );


    //- Destructor
    virtual ~ABLProfile()
    {}


    // Member Functions


        //- return distance field
        inline const Field<scalar> Z() const
        {
            return Z_->value();
        }


        //- Set boundary fields
        virtual void setProfiles
        (
            vectorField& U,
            scalarField& L,
            symmTensorField& R
        ) = 0;

            //- Set boundary fields
        virtual void setProfiles
        (
            vectorField& U,
            vectorField& L,
            symmTensorField& R
        ) = 0;

        //- Write
        virtual void write(Ostream&) const;

    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

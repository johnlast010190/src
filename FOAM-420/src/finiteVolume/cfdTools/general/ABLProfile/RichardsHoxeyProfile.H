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
    Foam::RichardsHoxeyProfile

Description
    RichardsHoxeyProfile for ABL inlet specification

SourceFiles
    RichardsHoxeyProfile.C

\*---------------------------------------------------------------------------*/

#ifndef RichardsHoxeyProfile_H
#define RichardsHoxeyProfile_H

#include "cfdTools/general/ABLProfile/ABLProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ABLProfiles
{

/*---------------------------------------------------------------------------*\
                          Class RichardsHoxeyProfile Declaration
\*---------------------------------------------------------------------------*/

class RichardsHoxeyProfile
:
    public ABLProfile
{
private:

    // Private Data

        //- Cmu constant
        const scalar Cmu_;

        //- Kappa constant
        const scalar Kappa_;

        //- Reference velocity
        const vector Uref_;

        //- Reference height
        const scalar Href_;

        //- Surface roughness height
        const scalar z0_;

        //- scaling for L to compute Lxx, Lyy, Lzz from
        //  L = k*sqrt(k)/epsilon * pow(Cmu_,0.75)
        const vector fLxyz_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const RichardsHoxeyProfile&);


public:

    //- Runtime type information
        TypeName("RichardsHoxey");

    // Constructors

        //- Null constructor
        RichardsHoxeyProfile(const fvPatch&);

        //- Construct from components
        RichardsHoxeyProfile(const fvPatch&, const dictionary&);

        //- copy constructor
        RichardsHoxeyProfile(const RichardsHoxeyProfile&);

        //- return a clone
        virtual ABLProfile* clone() const
        {
            return
            (
                new RichardsHoxeyProfile(*this)
            );
        }

    //- Destructor
    virtual ~RichardsHoxeyProfile(){}


    // Member Functions

        //- Set boundary fields
        virtual void setProfiles
        (
            vectorField& U,
            scalarField& L,
            symmTensorField& R
        );

        //- Set boundary fields
        virtual void setProfiles
        (
            vectorField& U,
            vectorField& L,
            symmTensorField& R
        );

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ABLProfiles
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::tabulatedProfile

Description
    tabulatedProfile

SourceFiles
    tabulatedProfile.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedProfile_H
#define tabulatedProfile_H

#include "cfdTools/general/ABLProfile/ABLProfile.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ABLProfiles
{

/*---------------------------------------------------------------------------*\
                          Class tabulatedProfile Declaration
\*---------------------------------------------------------------------------*/

class tabulatedProfile
:
    public ABLProfile
{
private:

    // Private Data

        //- Cmu constant
        const scalar Cmu_;

        //- Mean velocity profile
        autoPtr<Function1<vector>> Uf1_;

        //- Turbulent kinetic energy profile
        autoPtr<Function1<scalar>> kf1_;

        //- Turbulent intensity profile
        autoPtr<Function1<scalar>> if1_;

        //- Turbulent dissipation profile
        autoPtr<Function1<scalar>> epsilonf1_;

        //- Turbulent length scale profile
        autoPtr<Function1<scalar>> Lf1_;

        //- Turbulent length scale anisotropic profile
        // Lxx, Lyy, Lzz
        autoPtr<Function1<vector>> Lxyzf1_;

        //- Turbulent Reynolds stress
        autoPtr<Function1<symmTensor>> Rf1_;

        //- scaling for L to compute Lxx, Lyy, Lzz
        const vector fLxyz_;

    // Private Member Functions

        //- simple loop applicator
        template<class Type>
        void applyFunction1ToField
        (
            const Function1<Type>& data,
            const scalarField& Z,
            Field<Type>& field
        );

        //- calculate velocity gradient from velocity table
        tmp<tensorField> tableGradient
        (
            const Function1<vector>& U,
            const patchDistanceFunction& Z
        ) const;

        //- Disallow default bitwise assignment
        void operator=(const tabulatedProfile&);


public:

    //- Runtime type information
        TypeName("tabulated");

    // Constructors

        //- Null constructor
        tabulatedProfile(const fvPatch&);

        //- Construct from components
        tabulatedProfile(const fvPatch&, const dictionary&);

        //- copy constructor
        tabulatedProfile(const tabulatedProfile&);

        //- return a clone
        virtual ABLProfile* clone() const
        {
            return
            (
                new tabulatedProfile(*this)
            );
        }

    //- Destructor
    virtual ~tabulatedProfile(){}


    // Member Functions

        //- Set boundary fields
        virtual void setProfiles
        (
            vectorField& U,
            scalarField& L,
            symmTensorField& R
        );

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

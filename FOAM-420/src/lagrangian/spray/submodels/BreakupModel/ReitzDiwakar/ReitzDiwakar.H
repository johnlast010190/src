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

Class
    Foam::ReitzDiwakar

Group
    grpLagrangianIntermediateBreakupSubModels

Description
    secondary breakup model

    @verbatim
    Reitz, R.D.
    "Modelling atomization processes in highpressure vaporizing sprays"
    Atomization and Spray Technology 3 (1987), 309-337
    @endverbatim

    @verbatim
    Reitz, R.D. and Diwakar, R.
    "Effect of drop breakup on fuel sprays"
    SAE Tech. paper series, 860469 (1986)
    @endverbatim

    @verbatim
    Reitz, R.D. and Diwakar, R.
    "Structure of high-pressure fuel sprays"
    SAE Tech. paper series, 870598 (1987)
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef ReitzDiwakar_H
#define ReitzDiwakar_H

#include "submodels/BreakupModel/BreakupModel/BreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class ReitzDiwakar Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReitzDiwakar
:
    public BreakupModel<CloudType>
{
private:

    // Private data

        scalar Cbag_;
        scalar Cb_;
        scalar Cstrip_;
        scalar Cs_;


public:

    //- Runtime type information
    TypeName("ReitzDiwakar");


    // Constructors

        //- Construct from dictionary
        ReitzDiwakar(const dictionary&, CloudType&);

        //- Construct copy
        ReitzDiwakar(const ReitzDiwakar<CloudType>& bum);

        //- Construct and return a clone
        virtual autoPtr<BreakupModel<CloudType>> clone() const
        {
            return autoPtr<BreakupModel<CloudType>>
            (
                new ReitzDiwakar<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ReitzDiwakar();


    // Member Functions

        //- Update the parcel properties
        virtual bool update
        (
            const scalar dt,
            const vector& g,
            scalar& d,
            scalar& tc,
            scalar& ms,
            scalar& nParticle,
            scalar& KHindex,
            scalar& y,
            scalar& yDot,
            const scalar d0,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const vector& U,
            const scalar rhoc,
            const scalar muc,
            const vector& Urel,
            const scalar Urmag,
            const scalar tMom,
            scalar& dChild,
            scalar& massChild
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/BreakupModel/ReitzDiwakar/ReitzDiwakar.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

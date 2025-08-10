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
    Foam::NoAtomization

Group
    grpLagrangianIntermediateAtomizationSubModels

Description
    Dummy phase change model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoAtomization_H
#define NoAtomization_H

#include "submodels/AtomizationModel/AtomizationModel/AtomizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class NoAtomization Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoAtomization
:
    public AtomizationModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoAtomization(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        NoAtomization(const NoAtomization<CloudType>& am);

        //- Construct and return a clone
        virtual autoPtr<AtomizationModel<CloudType>> clone() const
        {
            return autoPtr<AtomizationModel<CloudType>>
            (
                new NoAtomization<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoAtomization();


    // Member Functions

        //- Flag to indicate whether model activates atomization model
        virtual bool active() const;

        //- Initial value of liquidCore
        virtual scalar initLiquidCore() const;

        //- Flag to indicate if chi needs to be calculated
        virtual bool calcChi() const;

        virtual void update
        (
            const scalar dt,
            scalar& d,
            scalar& liquidCore,
            scalar& tc,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const scalar volFlowRate,
            const scalar rhoAv,
            const scalar Urel,
            const vector& pos,
            const vector& injectionPos,
            const scalar pAmbient,
            const scalar chi,
            Random& rndGen
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/AtomizationModel/NoAtomization/NoAtomization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::PackingModel::NoPacking

Group
    grpLagrangianIntermediateMPPICPackingSubModels

Description
    Place holder for the 'none' option

SourceFiles
    NoPacking.C

\*---------------------------------------------------------------------------*/

#ifndef NoPacking_H
#define NoPacking_H

#include "submodels/MPPIC/PackingModels/PackingModel/PackingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PackingModels
{

/*---------------------------------------------------------------------------*\
                         Class NoPacking Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoPacking
:
    public PackingModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        NoPacking(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoPacking(const NoPacking<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<PackingModel<CloudType>> clone() const
        {
            return autoPtr<PackingModel<CloudType>>
            (
                new NoPacking<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoPacking();


    // Member Functions

        //- Calculate the velocity correction
        virtual vector velocityCorrection
        (
            typename CloudType::parcelType& p,
            const scalar deltaT
        ) const;

        //- Return the model 'active' status
        virtual bool active() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PackingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/MPPIC/PackingModels/NoPacking/NoPacking.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2021-2022 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::materialsSpecie

Description
    Perfect gas equation of state.

SourceFiles
    materialsSpecieI.H
    materialsSpecie.C

\*---------------------------------------------------------------------------*/

#ifndef materialsSpecie_H
#define materialsSpecie_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class materialsSpecie Declaration
\*---------------------------------------------------------------------------*/

class materialsSpecie
:
    public materialModel,
    public WModel,
    public YModel,
    public RModel
{
    // Private member data

        //- Name of materialsSpecieModel
        word name_;

        //- Mass fraction []
        scalar Y_;

        //- Molecular weight of materialsSpecieModel [kg/kmol]
        scalar molWeight_;

        scalar R_;

public:

    // Constructors

        //- Construct from dictionary and object registry
        materialsSpecie
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& materialsSpecieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<materialsSpecie> clone() const;

        //- Selector from dictionary
        inline static autoPtr<materialsSpecie> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& materialsSpecieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("specie");

        //- Molecular weight [kg/kmol]
        defineScalarFieldFunctionsZero(W);

        //- Mass fraction []
        defineScalarFieldFunctionsZero(Y);

        //- Gas constant [J/(kg K)]
        defineScalarFieldFunctionsZero(R);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update dependent model table
        virtual void updateTable(const word&) override
        {}

        //- Models casting
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        //- Re-read data
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "materialsSpecieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

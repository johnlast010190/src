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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2021-2023 Esi Ltd.

Class
    Foam::matSensibleEnthalpy

SourceFiles
    matSensibleEnthalpy.C

\*---------------------------------------------------------------------------*/

#ifndef matSensibleEnthalpy_H
#define matSensibleEnthalpy_H

#include "materialModels/materialModel/materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class matSensibleEnthalpy Declaration
\*---------------------------------------------------------------------------*/

class matSensibleEnthalpy
:
    public materialModel,
    public CpvModel,
    public CpByCpvModel,
    public HEModel
{
    // Private data

        //- Quick model access
        enum {Cp, Hs, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary
        matSensibleEnthalpy
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<matSensibleEnthalpy> clone() const;


    //- Destructor
    virtual ~matSensibleEnthalpy()
    {};


    // Member functions

        //- Return the instantiated type name
        TypeName("sensibleEnthalpy");

        //- Sensible internal energy type name
        static word name()
        {
            return "h";
        }

        //- Return whether Cp is uniform
        virtual bool isCpvConst() const override
        {
            return isModelConst<CpModel>();
        }

        virtual bool isHEConst() const override
        {
            return isModelConst<HsModel>();
        }

        //- Heat capacity at constant volume [J/(kg K)]
        declareScalarFieldFunctions(Cpv);

        //- Cp/Cp []
        defineScalarFieldFunctionsZero(CpByCpv);

        //- Sensible internal energy [J/kg]
        declareScalarFieldFunctions(HE);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update table
        virtual void updateTable(const word& modelName) override;

        //- Models casting
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        //- Re-read data
        virtual bool read() override
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "matSensibleEnthalpyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

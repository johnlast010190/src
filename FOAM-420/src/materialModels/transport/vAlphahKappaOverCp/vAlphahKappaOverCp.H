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
    (c) 2022-2023 Esi Ltd.

Class
    Foam::vAlphahKappaOverCp

Description


SourceFiles
    vAlphahKappaOverCp.C

\*---------------------------------------------------------------------------*/

#ifndef vAlphahKappaOverCp_H
#define vAlphahKappaOverCp_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class vAlphahKappaOverCp Declaration
\*---------------------------------------------------------------------------*/

class vAlphahKappaOverCp
:
    public materialModel,
    public vAlphahModel
{
    // Private data

        //- Quick model access
        enum {Cp, modelsEnumSize_};

        enum {vKappaInd, modelsEnumSize2_};


public:

    // Constructors

        //- Construct from dictionary
        vAlphahKappaOverCp
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<vAlphahKappaOverCp> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("vAlphahKappaOverCp");

        bool isvAlphahConst() const override
        {
            return isModelConst<kappaModel>() && isModelConst<CpModel>();
        }

        //- Thermal diffusivity of enthalpy [kg/ms]
        declareVectorFieldFunctions(vAlphah);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update table
        virtual void updateTable(const word& modelName) override;

        //- Models casting
        virtual baseModels<vector>* castVectorModel
        (
            const word& modelName
        ) override;

        //- Re-read data
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

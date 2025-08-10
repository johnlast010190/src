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
    (c) 2017-2023 Esi Ltd.

Class
    Foam::tabulated

Description

    Calculating any rabulated properties:

    //- required
    <model>TableCoeffs
    {
        file        "$FOAM_CASE/constant/mu.dat";
        outOfBounds     clamp;
    }

    The tables are described in terms of
    the field1 and field2

    For example, the table:

                           field2
                       |  300   450   600
                  -----+-----------------
                  1e5  |   A     B     C
        field1  1.1e5  |   D     E     F
                  2e5  |   G     H     I


    Is specified by the following:

       (
            (
                1e5
                (
                    (300    A)
                    (450    B)
                    (600    C)
                )
            )
            (
                1.1e5
                (
                    (300    D)
                    (450    E)
                    (600    F)
                )
            )
            (
                2e5
                (
                    (300    G)
                    (450    H)
                    (600    I)
                )
            )
        )

    All field1 and field2 must be in ascending order
    in each list.

SourceFiles
    tabulatedI.H
    tabulated.C

\*---------------------------------------------------------------------------*/

#ifndef tabulated_H
#define tabulated_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/generalModel.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class tabulated Declaration
\*---------------------------------------------------------------------------*/

class tabulated
:
    public materialModel,
    public generalModel
{
    // Private data

        //- 2D look up table = function of two fields (p, T)
        autoPtr<interpolation2DTable<scalar>> table_;

        //- Name of the first field
        word fieldName1_;

        //- Name of the scond field
        word fieldName2_;

        //- Pointer to the temperature offset by reference
        const refScalarField* field1_;

        //- Pointer to the pressure offset by the reference
        const refScalarField* field2_;

        //- Is rho dependent on f(p)
        bool incompressible_;

        //- Is the equation of state is isochoric i.e. rho = const
        bool isochoric_;

        //- Is the model isotropic?
        bool isotropic_;


public:

    // Constructors

        //- Construct from dictionary
        tabulated
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<tabulated> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("tabulated");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override
        {
            return incompressible_;
        }

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override
        {
            return isochoric_;
        }

        //- Is it isotropic model?
        virtual bool isotropic() const override
        {
            return isotropic_;
        }

        //- Whether the function is uniform
        virtual bool isConst() const override
        {
            return false;
        }

        //- Definition of the function for single value
        virtual inline scalar general
        (
            const scalar val1,
            const scalar val2
        ) const;

        //- Cell calculations
        virtual inline scalar generalCell
        (
            const label celli
        ) const override;

        //- Internal field calculation
        virtual tmp<scalarField> generalInternal() const override;

        //- Calculate values for single patch
        virtual tmp<scalarField> generalPatch
        (
            const label patchi
        ) const override;

        //- Calculation of whole volumetric field
        virtual tmp<volScalarField> generalGeometric() const override;

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update table
        virtual void updateTable(const word& modelName) override
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

#include "tabulatedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
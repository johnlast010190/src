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
    (c) 2022-2023 Esi Ltd.

Class
    Foam::fv::fixedTemperatureConstraint

Group
    grpFvOptionsConstraints

Description
    Fixed temperature equation constraint

Usage
    \verbatim

    // Temperature that can vary in time
    fixedTemperature
    {
        type            fixedTemperatureConstraint;
        active          yes;
        mode            uniform;
        temperature     500; // fixed temperature with time [K]
    }

    // Temperature fixed by lookup from current temperature
    fixedTemperature
    {
        type            fixedTemperatureConstraint;
        active          yes;
        mode            lookup;
        // T            <Tname>;      // optional temperature field name
    }

    // Example spatial equation:
    // T = 300 + 0.5x + 3y + 0.3xy^0.5
    fixedTemperature
    {
        type            fixedTemperatureConstraint;
        active          yes;
        mode            spatial;
        selectionMode   cellZone;
        cellZone        c0;

        temperature xyzPolynomial
        (
            (300 (0 0 0))   // 300
            (0.5 (1 0 0))   // 0.5*x
            (3 (0 1 0))     // 3*y
            (0.3 (1 0.5 0)) // 0.3*x*y^0.5
        );

        // Optional coordinate system
        referenceFrame frameName;
    }
    \endverbatim

SourceFiles
    fvOption.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperatureConstraint_H
#define fixedTemperatureConstraint_H

#include "cellSetOption/cellSetOption.H"
#include "primitives/enums/NamedEnum.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "basicThermo/basicThermo.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class fixedTemperatureConstraint Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperatureConstraint
:
    public cellSetOption
{

public:

    //- Temperature mode
    enum temperatureMode
    {
        tmUniform,
        tmLookup,
        tmSpatial
    };


    //- String representation of temperatureMode enums
    static const NamedEnum<temperatureMode, 3> temperatureModeNames_;

    void setMatrixFieldName
    (
        const word& eqnName,
        label fieldi = 0
    );


protected:

    // Protected data

        //- Operation mode
        temperatureMode mode_;

        //- Uniform temperature [K]
        autoPtr<Function1<scalar>> Tfun1_;

        //- Temperature field name
        word TName_;

        //- Thermo pointer
        const basicThermo* thermoPtr_;

        //- Reference frame
        coordinateFrame* coorFramePtr_;

        //- Fraction of the constraint to apply. Facilitates ramping, or
        //  pulsing, or deactivation after a time. Should take a value between
        //  0 and 1. Defaults to 1 (i.e., apply integral constraint).
        autoPtr<Function1<scalar>> fraction_;


private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const fixedTemperatureConstraint&);


public:

    //- Runtime type information
    TypeName("fixedTemperatureConstraint");


    // Constructors

        //- Construct from components
        fixedTemperatureConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~fixedTemperatureConstraint()
    {}


    // Member Functions

        //- Constrain energy equation to fix the temperature
        virtual void constrain
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) override;

        //- Read dictionary
        virtual bool read(const dictionary& dict) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

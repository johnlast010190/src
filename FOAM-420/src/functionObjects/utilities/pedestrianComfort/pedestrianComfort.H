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
    (c) 2020 Esi Ltd.

Class
    Foam::functionObjects::pedestrianComfort

Group
    grpUtilitiesFunctionObjects

Description
Pedestrian comfort criteria implememnted based on
"Wind Micro-climate guidelines for developments in the City of London"
August 2019 - Authors: City of London & RWDI

Supported only by single mesh multi-region solver

Usage

    Example of function object specification:
    \verbatim
    Lawson
    {
        type pedestrianComfort;
        regions (dir1);
        functionObjectLibs ("libutilityFunctionObjects.so");
        writeControl    writeTime;
        UName U;
        meanVelCorr false; //based on AIJ paper Sec 3.1.4.4
        gust false;
        windDirections
        (
            dir1
            dir2
            dir3
            dir4
        );
        WeibullParameters
        (
            (2.1 4.63 1.70)
            (2.0 5.06 1.80)
            (2.4 5.40 1.87)
            (2.6 5.65 1.92)
        );
        referenceHeight
        (
            100
            200
            300
            100
        );
        referenceVelocity
        (
            10
            10
            10
            10
        );
        //city of London
        scaleFactor table
        (
            (100 0.96)
            (120 1)
            (160 1.07)
            (200 1.13)
            (250 1.19)
            (300 1.24)
            (450 1.37)
            (600 1.48)
        );
        criteriaThresholds
        (
            2.5
            4
            6
            8
        );
        criteriaExceedance // 5%
        (
           5
           5
           5
           5
        );
    };
    \endverbatim

SourceFiles
    pedestrianComfort.C

\*---------------------------------------------------------------------------*/

#ifndef pedestrianComfort_H
#define pedestrianComfort_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "fields/volFields/volFields.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
            Class pedestrianComfortFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class pedestrianComfort
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        // velocity field name
        word UName_;

        // k field name (optional) -- default is k
        word kName_;

        // mean velocity correction based on 2*tke
        bool meanVelCorr_;

        // use gust equivalent mean velocity instead of
        // standard velocity
        bool gustCorr_;
        // peak factor
        scalar gustG_;
        // gust equivalent mean factor
        scalar gustT_;

        // list of wind directions (regions name)
        List<word> windDirections_;

        // list of solution instance names
        List<word> instanceNames_;

        // Number of wind directions
        label numberOfWindDirections_;

        // Weibull parameters for each wind direction
        // p [%] probability wind approaching from specific dir
        // c scale factor [m/s]
        // k shape factor
        List<Vector<scalar>> WeibullParameters_;

        // reference height for each wind direction
        List<scalar> hRef_;

        // reference velocity for each wind direction
        List<scalar> URef_;

        // reference height scale factor
        // ref height [m] vs scale factor
        autoPtr<Function1<scalar>> scaleFactor_;

        // list of velocities from Lawson criteria
        // or safety criteria used to evaluate
        // exceedance probabilities i.e. for London City
        // Lawson Criterial (2.5 4 6 8 ) m/s
        List<scalar> criteriaThresholds_;

        // criteria probability of exceedance,
        // i.e. Lawson 5%, safety 0.022 %
        List<scalar> criteriaExceedance_;

        // List of pointers to total probability
        // of exceedance for each threshold velocity
        List<volScalarField*> totProbExceedPtr_;

        // Pedestrian Comfort field
        // with marker for each comfort category.
        autoPtr<volScalarField> pComfort_;

private:
        // Calculate total probability exceedance
        // for each velocity threshold specified
        void calculateTotalProbabilityExceedance();

        // Compare total probabilities exceedance
        // against the wind comfort criteria and
        // created Pedestrian comfort field with
        // markers for each comfort category
        void assessComfort();

public:

    //- Runtime type information
    TypeName("pedestrianComfort")

    // Constructors

        pedestrianComfort
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

    virtual ~pedestrianComfort();

    // Member Functions

        //- Execute, currently does nothing
        bool execute() override;

        //- Sample and write
        bool write() override;

        //- Read the parameters
        bool read(const dictionary&) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


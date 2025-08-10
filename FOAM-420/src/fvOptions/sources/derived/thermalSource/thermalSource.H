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
    (c) 2013-2022 Esi Ltd.

Class
    Foam::fv::thermalSource

Group
    grpFvOptionsSources

Description
    Simple volumetric heat source/sink model. Operates in three modes:
    1. Uniform fixed heat transfer, Q_t [W]
    2. Fixed heat transfer coefficient (HTC), h [W/m3/K]
    3. Fixed total heat transfer, with scaled HTC. Heat transfer coefficient
       is calculated such that the total heat transfer Q_t [W] is satisfied.

    This is a very basic heat source designed to stop excessive heat
    collecting in flow dead spots. Effectiveness heat exchanger should be used
    where more accuracy is required.

    The total heat transfer rate is given by:
    \f[
        Q_t = \dot{h} \sum((T_{wm} - T_{air}) dV)
    \f]

    where:
    \vartable
        Q_t       | total heat source [W]
        \dot{h}   | uniform heat transfer coefficient [W/m3/K]
        T_{wm}    | mean water-side temperature [K]
        dV        | cell volume
    \endvartable

    The discretisation is abstracted as:
    \f[
        Q_t = Q_{expl} + Q_{impl}
            = htc_/VDash_ (T_{wm} - T_{air})
            + htcScaled he - fvm::Sp(htcScaled, he)
    \f]

    with
    \f[
        htcScaled = htc_ if he = "T"
        htcScaled = htc_/Cpv if he = "h"/"e"
    \f]

    \heading Source usage

    Example usage:
    \verbatim
    thermalSource1
    {
        type            thermalSource;
        active          yes;

        thermalSourceCoeffs
        {
            selectionMode    cellZone;
            cellZone         porosity;
            volumeMode       absolute;  // specific
            heatTransferMode scaledHTC; // uniformFlux / fixedHTC
            fieldName        T; // h / e

            Q                1000;
            Tinf             336;
            htc              1;
            dTmin            1;//[K] Minimum mean difference in temperature
                                 //Optional, sets limit for scaled HTC
        }
    }
    \endverbatim


Note
 - volumeMode:
    - "specific" all coefficients operate as normal
    - "absolute" htc [W/K]
    - Both htc and Q must be specified in the correct mode/units
 - heatTransferMode "uniformFlux" will result in an explicit source.
 - heatTransferMode "fixedHTC" will result in a semi-implicit source.
 - heatTransferMode "scaledHTC" will result in a semi-implicit source with the
   heat transfer coefficient scaled to provide the target heat flux.

SourceFiles
    thermalSource.C

\*---------------------------------------------------------------------------*/

#ifndef thermalSource_H
#define thermalSource_H

#include "cellSetOption/cellSetOption.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class thermalSource Declaration
\*---------------------------------------------------------------------------*/

class thermalSource
:
    public cellSetOption
{

public:

    // Public data

        //- Enumeration for volume types
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Word list of volume mode type names
        static const NamedEnum<volumeModeType, 2> volumeModeTypeNames_;

        //- Enumeration for heat transfer types
        enum heatTransferModeType
        {
            hmUniformFlux,
            hmFixedHTC,
            hmScaledHTC,
            hmScaledTinf
        };

        //- Word list of volume mode type names
        static const NamedEnum<heatTransferModeType, 4> heatTransferModeTypeNames_;

protected:

    // Protected data

        //- Volume mode
        volumeModeType volumeMode_;

        //- Volume normalisation
        scalar VDash_;

        //- Heat transfer mode
        heatTransferModeType heatTransferMode_;

        //- Heat flux [W] or [W/m3]
        autoPtr<Function1<scalar>> QF1_;
        scalar Q_;

        //- Heat transfer coefficient [W/K] or[W/m3/K]
        autoPtr<Function1<scalar>> htcF1_;
        scalar htc_;

        //- Mean bulk temperature [K]
        autoPtr<Function1<scalar>> TinfF1_;
        scalar Tinf_;

        //- Minimum aggregate temperature differential / limiter for auto-htc
        scalar DTmin_;

        //- Name of temperature field
        word TName_;

        //- Reporting y/n
        bool reportOn_;


    // Protected functions

        //- Scale heat transfer coefficient based on he dimensions
        tmp<scalarField> htcScaled(const volScalarField&) const;

        //- Set the local field data
        void updateRate(const fvMatrix<scalar>& eqn);


public:

    //- Runtime type information
    TypeName("thermalSource");


    // Constructors

        //- Construct from components
        thermalSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~thermalSource()
    {}


    // Member Functions

        // Evaluation

            //- Add contribution to equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );

            //- Add contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldI
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

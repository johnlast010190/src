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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::fv::tabulatedEffectivenessHeatTransfer

Group
    grpFvOptionsSources

Description
    Tabulated effectiveness heat transfer model.

    Heat transfer coefficient estimated by:

    \f[
        htc = \frac{q_{target}}{dT V_{core}}
    \f]

    Where \f$ q_{target} \f$ is given by:

    \f[
        q_{target} = effectiveness C_{min} (T_{in,1} - T_{in,2})
        C_{min} = min \left(Cp_1 \dot{m}_1, Cp_2 \dot{m}_2 \right)
    \f]

    \heading Example usage

    \verbatim
    coolerToAir
    {
        type            tabulatedEffectivenessHeatTransfer;
        active          yes;

        tabulatedEffectivenessHeatTransferCoeffs
        {
            interpolationMethod cellVolumeWeight;
            nbrRegion       air;
            master          true;

            fields          (h);
            outOfBounds     clamp;
            file            "dataTable";
            nbrModel        airToWater;
            semiImplicit    no;

            inletPatch      inlet_water;
            inletPatchNbr   inlet_air;
        }
    }
    \endverbatim

    The effectiveness table is described in terms of the primary (=master)
    and secondary mass flow rates.  For example, the table:

    \verbatim
                           secondary MFR
                       |  0.1   0.2   0.3
                  -----+-----------------
                  0.02 |   A     B     C
     primary MFR  0.04 |   D     E     F
                  0.06 |   G     H     I
    \endverbatim


    Is specified by the following:

    \verbatim
        (
            0.02
            (
                (0.1    A)
                (0.2    B)
                (0.3    C)
            ),
            0.04
            (
                (0.1    D)
                (0.2    E)
                (0.3    F)
            ),
            0.06
            (
                (0.1    G)
                (0.2    H)
                (0.3    I)
            )
        );
    \endverbatim

SourceFiles
    tabulatedEffectivenessHeatTransfer.C

SeeAlso
    interRegionHeatTransferModel.H

\*---------------------------------------------------------------------------*/

#ifndef tabulatedEffectivenessHeatTransfer_H
#define tabulatedEffectivenessHeatTransfer_H

#include "sources/interRegion/interRegionHeatTransfer/interRegionHeatTransferModel/interRegionHeatTransferModel.H"
#include "memory/autoPtr/autoPtr.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class tabulatedEffectivenessHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class tabulatedEffectivenessHeatTransfer
:
    public interRegionHeatTransferModel
{
private:

    // Private data

        //- Name of velocity field; default = U
        word UName_;

        //- Name of neighbour velocity field; default = U
        word UNbrName_;

        //- Name of density field; default = rho
        word rhoName_;

        //- Name of neighbour density field; default = rho
        word rhoNbrName_;

        //- Name of temperature field; default = T
        word TName_;

        //- Name of neighbour temperature field; default = T
        word TNbrName_;

        //- Pointer to 2-D look-up table of effectiveness f(mDot1, mDot2)
        autoPtr<interpolation2DTable<scalar>> effectivenessTable_;

        //- Heat exchanger core volume
        scalar Vcore_;

        //- Inlet area [kg/s]
        scalar mDot_;

        //- Neighbour region inlet area [kg/s]
        scalar mDotNbr_;

        //- Fluid Cp value
        scalar Cp_;

        //- Neighbour fluid Cp value
        scalar CpNbr_;

        //- Inlet temperature [K]
        scalar Tin_;

        //- Neighbour region inlet temperature [K]
        scalar TinNbr_;

        //- Inlet patch name
        word inletPatch_;

        //- Neighbour region inlet patch name
        word inletPatchNbr_;

        //- Model target energy source value
        scalar Qtarget_;


    // Private Member functions

        //- Effectiveness table helper
        const interpolation2DTable<Foam::scalar>& effectivenessTable();

        //- Thermophysical properties helper
        const basicThermo& thermo(const fvMesh& mesh) const;

        //- Initialise geometry
        void initialiseGeometry();

        //- Calculate and return mass flow rate
        scalar massFlowRate() const;
        scalar massFlowRateNbr() const;

        //- Calculate and return inlet temperature
        scalar inletTemperature() const;
        scalar inletTemperatureNbr() const;

        //- Calculate and return inlet Cp
        scalar inletCp() const;
        scalar inletCpNbr() const;


public:

    //- Runtime type information
    TypeName("tabulatedEffectivenessHeatTransfer");


    // Constructors

        //- Construct from dictionary
        tabulatedEffectivenessHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~tabulatedEffectivenessHeatTransfer();


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc(scalar deltaT);

        //- Set the target heat transfer
        scalar calcTargetEnergy();
        virtual scalar getTargetEnergy();


        // I-O

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

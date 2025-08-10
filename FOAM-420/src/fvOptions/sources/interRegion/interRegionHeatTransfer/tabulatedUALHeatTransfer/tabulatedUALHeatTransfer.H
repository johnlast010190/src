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
    Foam::fv::tabulatedUALHeatTransfer

Group
    grpFvOptionsSources

Description
    Tabulated UAL heat transfer model.

    The heat flux is calculated based on the mass transfer coefficient UAL.
    A 2-D table of UAL as functions of the primary and secondary mass flow
    rates is required.

    The exchanger geometry can be specified using:
        - user must provide either one of below
           - inlet area of each region
           - mass flow rate of each region
           - inlet patch name for each region from which mDot is calculated
        - core volume automatically calculated by the mesh-to-mesh
          interpolation volume if not provided by user
        - Cp can be specified

    Heat transfer coefficient calculated by:

    \f[
        htc = \frac{UAL}{V_{core}}
    \f]

    \heading Example usage

    \verbatim
    coolerToAir
    {
        type            tabulatedUALHeatTransfer;
        active          yes;

        tabulatedUALHeatTransferCoeffs
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

    The NTU table is described in terms of the primary (=master)
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
    tabulatedUALHeatTransfer.C

SeeAlso
    interRegionHeatTransferModel.H

\*---------------------------------------------------------------------------*/

#ifndef tabulatedUALHeatTransfer_H
#define tabulatedUALHeatTransfer_H

#include "sources/interRegion/interRegionHeatTransfer/interRegionHeatTransferModel/interRegionHeatTransferModel.H"
#include "memory/autoPtr/autoPtr.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class tabulatedUALHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class tabulatedUALHeatTransfer
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

        //- Pointer to 2-D look-up table of UAL f(mDot1, mDot2)
        autoPtr<interpolation2DTable<scalar>> ualTable_;

        //- Inlet area [m2]
        scalar Ain_;

        //- Neighbour region inlet area [m2]
        scalar AinNbr_;

        //- Heat exchanger core volume
        scalar Vcore_;

        //- Inlet area [kg/s]
        scalar mDot_;

        //- Neighbour region inlet area [kg/s]
        scalar mDotNbr_;

        //- Inlet patch name
        word inletPatch_;

        //- Neighbour region inlet patch name
        word inletPatchNbr_;


    // Private Member functions

        //- UAL table helper
        const interpolation2DTable<Foam::scalar>& ualTable();

        //- Initialise geometry
        void initialiseGeometry();

        //- Calculate and return mass flow rate
        tmp<scalarField> massFlowRate() const;
        tmp<scalarField> massFlowRateNbr() const;


public:

    //- Runtime type information
    TypeName("tabulatedUALHeatTransfer");


    // Constructors

        //- Construct from dictionary
        tabulatedUALHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~tabulatedUALHeatTransfer();


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc(scalar deltaT);


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

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
    (c) 2015 OpenCFD Ltd.
    (c) 2018 Esi Ltd.

Class
    Foam::fv::directionalPressureGradientExplicitSource

Group
    grpFvOptionsSources

Description
    Creates an explicit  pressure gradient source in such a way to deflect the
    flow towards an specific direction (flowDir). Alternatively add an extra
    pressure drop in the flowDir direction using a model.

    \heading Source usage
    Example usage:
    \verbatim
    airDeflection
    {
        type            directionalPressureGradientExplicitSource;
        active          true;

        directionalPressureGradientExplicitSourceCoeffs
        {
            selectionMode   cellZone;
            cellZone        cZone;

            fields      (U);            // Name of the field
            flowDir     (1 1 0);        // Desired flow direction
            faceZone    f0Zone;         // Face zone upstream cell zone
            relaxationFactor    0.3;    // Relaxation factor for flow
                                        // deflection (default 0.3)

            //Pressure drop model [Pa]
            model       volumetricFlowRate;//pressureDrop, DarcyForchheimer;

            //DarcyForchheimer model
            // deltaP = (D*mu + 0.5*rho*magUn)*magUn*length_

            D           5e7;
            I           0;
            length      1e-3;

            //Pressure drop Function1 interface
            //time dependent if model == pressureDrop
            //flow rate dependent if model == volumetricFlowRate
            pressureDrop    constant 40;
        }
    }
    \endverbatim

    NOTE: In order to obtain the upwind velocities this function loops over
    the slaves cells of the faceZone specified in the dictionary, on the other
    hand, the cellZone to which this source term is applied should be composed
    of the master cells and they should be 'downwind' the faceZone.

SourceFiles
    directionalPressureGradientExplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef directionalPressureGradientExplicitSource_H
#define directionalPressureGradientExplicitSource_H

#include "cellSetOption/cellSetOption.H"
#include "memory/autoPtr/autoPtr.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "cfdTools/general/fvOptions/fvOption.H"
#include "primitives/functions/Function1/Function1/Function1.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
          Class directionalPressureGradientExplicitSource Declaration
\*---------------------------------------------------------------------------*/

class directionalPressureGradientExplicitSource
:
    public cellSetOption
{
public:

    //  Public enumeration

        //- Modes of pressure drop
        enum pressureDropModel
        {
            pVolumetricFlowRate,
            pPressureDrop,
            pDarcyForchheimer
        };


private:

    // Private data

        static const NamedEnum<pressureDropModel, 3> PressureDropModelNames_;

         //- Pressure drop model
        pressureDropModel model_;

        //- Pressure gradient before correction
        vectorField gradP0_;

        //- Change in pressure gradient
        vectorField dGradP_;

        //- Pressure drop due to porous media
        vectorField gradPporous_;

        //- Flow direction
        vector flowDir_;

        //- Matrix 1/A coefficients field pointer
        autoPtr<volScalarField> invAPtr_;

        //- Darcy pressure loss coefficient
        scalar D_;

        //- Inertia pressure lost coefficient
        scalar I_;

        //- Porous media length
        scalar length_;

        //- Constant pressure drop
        autoPtr<Function1<scalar>> pressureDrop_;

        //- Name of the faceZone at the heat exchange inlet
        word faceZoneName_;

        //- Id for the face zone
        label zoneID_;

         //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch ID per face
        labelList facePatchId_;

        //- Relaxation factor
        scalar relaxationFactor_;

        //- Cells faces mapping
        labelList cellFaceMap_;


    // Private Member Functions

        //- Init
        void init();

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const vectorField& gradP) const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn);

        //- Disallow default bitwise copy construct
        directionalPressureGradientExplicitSource
        (
            const directionalPressureGradientExplicitSource&
        );

        //- Disallow default bitwise assignment
        void operator=(const directionalPressureGradientExplicitSource&);


public:

    //- Runtime type information
    TypeName("directionalPressureGradientExplicitSource");


    // Constructors

        //- Construct from explicit source name and mesh
        directionalPressureGradientExplicitSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions

        // Evaluate

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add implicit contribution to block Momentum equation
            virtual void addSup
            (
                fvBlockMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to block Momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvBlockMatrix<vector>& eqn,
                const label fieldi
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvBlockMatrix<vector>& eqn,
                const label fieldI
            );


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

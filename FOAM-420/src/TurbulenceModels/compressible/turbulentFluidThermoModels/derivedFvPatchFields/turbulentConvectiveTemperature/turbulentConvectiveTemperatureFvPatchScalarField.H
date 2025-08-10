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
    (c) 2010-2022 Esi Ltd.
    (c) 1991-2009 OpenCFD Ltd.

Class
    Foam::turbulentConvectiveTemperatureFvPatchScalarField

Description
    Foam::turbulentConvectiveTemperatureFvPatchScalarField

    Temperature boundary conditions for resistive wall heat transfer
    Includes radiation if present.

    alphaWall is defined as Q = alphaWall(Twall - Tinf) [W/m2]
    lambda is thermal conductivity

    @verbatim
        wall
        {
            type            turbulentConvectiveTemperature;
            Tinf            uniform 300;
            //alphaWall     uniform 0.01;

            //layers
            //{
            //  outside_metal
            //  {
            //      thickness uniform 0.003;
            //      lambda    uniform 16;
            //  }
            //
            //  foam
            //  {
            //      thickness uniform 0.05;
            //      lambda    uniform 0.03;
            //  }
            //
            //  inside_metal
            //  {
            //      thickness uniform 0.002;
            //      lambda    uniform 16;
            //  }
            //}
        }
    @endverbatim

SourceFiles
    turbulentConvectiveTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentConvectiveTemperatureFvPatchScalarFields_H
#define turbulentConvectiveTemperatureFvPatchScalarFields_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "turbulentFluidThermoModels/derivedFvPatchFields/boundaryKappa/boundaryKappa.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
      Class turbulentConvectiveTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentConvectiveTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public boundaryKappa
{
    // Private data

        //- Tinf
        scalarField Tinf_;

        //- alphaWall (assumed to be infinite if it doesn't exist)
        autoPtr<scalarField> alphaWall_;

        //- Old optional heat flux
        autoPtr<scalarField> QoptOld_;

        //- Relaxation value
        scalar relax_;

        //- material layers
        //- layer names
        List<word> layerNames_;

        //- layer thickness
        PtrList<scalarField> t_;

        //- layer conductivity
        PtrList<scalarField> lambda_;

        //- total conductivity
        mutable autoPtr<scalarField> at_;

        //- additional heat flux
        autoPtr<scalarField> qadd_;

    // Private Member Functions

        void calculateTotalConductivity() const;

public:

    //- Runtime type information
    TypeName("compressible::turbulentConvectiveTemperature");


    // Constructors

        //- Construct from patch and internal field
        turbulentConvectiveTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentConvectiveTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // turbulentConvectiveTemperatureFvPatchScalarField
        // onto a new patch
        turbulentConvectiveTemperatureFvPatchScalarField
        (
            const turbulentConvectiveTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentConvectiveTemperatureFvPatchScalarField
        (
            const turbulentConvectiveTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentConvectiveTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentConvectiveTemperatureFvPatchScalarField
        (
            const turbulentConvectiveTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentConvectiveTemperatureFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        // Access

            //- Return Tinf
            const scalarField& Tinf() const
            {
                return Tinf_;
            }

            //- Return reference to Tinf to allow adjustment
            scalarField& Tinf()
            {
                return Tinf_;
            }

            //- Return alphaWall
            const scalarField& alphaWall() const
            {
                return alphaWall_();
            }

            //- Return reference to alphaWall to allow adjustment
            scalarField& alphaWall();

            //- Return layer names
            const List<word>& layerNames() const
            {
                return layerNames_;
            }

            //- Return reference to layer names to allow adjustment
            List<word>& layerNames();

            //- Return thickness
            const PtrList<scalarField>& thickness() const
            {
                return t_;
            }

            //- Return reference to thickness to allow adjustment
            PtrList<scalarField>& thickness();

            //- Return conductivity
            const PtrList<scalarField>& lambda() const
            {
                return lambda_;
            }

            //- Return reference to conductivity to allow adjustment
            PtrList<scalarField>& lambda();

            //- Return total conductivity
            const scalarField& conductivity() const;

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

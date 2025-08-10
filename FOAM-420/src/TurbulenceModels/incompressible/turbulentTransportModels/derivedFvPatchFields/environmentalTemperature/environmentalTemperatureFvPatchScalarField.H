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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2009 OpenCFD Ltd.

Class
    Foam::environmentalTemperatureFvPatchScalarField

Description
    Foam::environmentalTemperatureFvPatchScalarField

    Temperature boundary conditions outside environment.

    Heat transfer is defined as:

    Qoutside = alphaOutside*(Tconv - Twall) + eps*sigma*(Trad^4 - Twall^4)
        + Qrad                                                    [W/m2]
    or
    Qoutside = Q //if Q is specified all else ignored             [W/m2]

    Qinside = alphaInside*(Twall - T0)                            [W/m2]

    No internal radiation for now. ***combine with layer BC later***

    @verbatim
        option1
        {
            q uniform 10; //fixed heat flux W/m2
        }

        option2
        {
            Twall 300; //fixed uniform temperature
        }

        option3
        {
            Tenv 300; //fixed uniform convection temperature
            alphaConv 0.1; //fixed uniform convection coefficient

            //alphaName lambda; //name of mesh side thermal difusivity
            //CpName Cp; //name of mesh side heat capacitance

            //Trad 300; //optional fixed environemtnal radiation temp
            //emissivity 1;//optional emissivity of this and environment
            //sigma 5.6704e-8;//Stefan Boltzmann constant

            //directed radiation optional
            directedRadiation
            (
                sun
                {
                    intensity 100; //W/m2
                    direction (0 0 -1);
                }

                skyEast
                {
                    intensity 50; //W/m2
                    direction (-1 0 -1);
                }
            );

        }


    @endverbatim


SourceFiles
    environmentalTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef environmentalTemperatureFvPatchScalarFields_H
#define environmentalTemperatureFvPatchScalarFields_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
      Class environmentalTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class environmentalTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- environmental radiation temperature, K
        autoPtr<scalar> Trad_;

        //- emmissivity
        scalar eps_;

        //- Boltzman constant
        scalar sigma_;

        //- environmental convection/conduction temperature, K
        autoPtr<scalar> Tconv_;

        //- outside heat transfer coefficient
        autoPtr<scalar> alphaOutside_;

        //- inside heat transfer coefficient field name
        word alphaName_;
        word CpName_;

        //- environemtnal radiation settings and heat flux
        PtrList<entry> rSourcesPtr_;
        autoPtr<scalarField> Qenv_;

        //- fixed heat flux value [W/m2]
        autoPtr<scalarField> QPtr_;


    // Private Member Functions

        //- calculate environmental radiation flux
        void calculateEnvRadiation();


public:

    //- Runtime type information
    TypeName("incompressible::environmentalTemperature");


    // Constructors

        //- Construct from patch and internal field
        environmentalTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        environmentalTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // environmentalTemperatureFvPatchScalarField
        // onto a new patch
        environmentalTemperatureFvPatchScalarField
        (
            const environmentalTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        environmentalTemperatureFvPatchScalarField
        (
            const environmentalTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new environmentalTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        environmentalTemperatureFvPatchScalarField
        (
            const environmentalTemperatureFvPatchScalarField&,
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
                new environmentalTemperatureFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        // Access


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

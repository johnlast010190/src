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
    (c) 2010-2016 Esi Ltd.
    (c) 2015 OpenFOAM Foundation

Class
    Foam::radiation::domSolar

Description

    A solar and environmental radiation model based on a simple finite volume
    flux convection scheme.

    The domSolar model can be used in conjuntion with fvDOM
    radiation model but only using a single band spectrum.

    The model is aplicable to non-interacting media only.

    @verbatim
    domSolarCoeffs
    {
        sources
        (
            sun
            {
                radiantInputType userDefined;
                intensity constant 10;
                direction consant (0 -1 -1);
            }
            sky
            {
                radiantInputType userDefined;
                intensity 15; direction (0 0 -1); solve false; initTol 1e-4;
                localSourceParches ( skylight ); // overwrite sourcePatches if specified
            }
        );
        // solve and initTol are optional, default values displayed
        // initTol = absolute tolerance of maximum residual for initial
        // field generation with explicit wave solver
        // solve = whether to use higher order implicit solver to refine
        // env radiation. Especially useful for solar or other highly
        // directional radiation


        sourcePatches
        (
            "window.*" "*.window" skylight
        );
    }


SourceFiles
    domSolar.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModeldomSolar_H
#define radiationModeldomSolar_H

#include "radiationModels/radiationModel/radiationModel.H"
#include "fvMesh/singleCellFvMesh/singleCellFvMesh.H"
#include "fields/volFields/volFields.H"
#include "submodels/solarCalculator/solarCalculator.H"
#include "radiationModels/domSolar/radiantField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class domSolar Declaration
\*---------------------------------------------------------------------------*/

class domSolar
:
    public radiationModel,
    public PtrList<radiantField>
{
protected:

    // Static data

private:

    // Private data


        //- Radiation fields / [W/m2]
        PtrList<radiantField> IenvPtr_;

        //- Incident environmental radiative heat flux [W/m2]
        volScalarField Qenv_;

        //- patch source coefficients
        FieldField<Field, scalar> patchSourceCoeffs_;

        //- patch source coefficients for local sources
        PtrList<FieldField<Field, scalar>> patchLocalSourceCoeffs_;

        //- role switch between supplemental and primary radiation model
        Switch primary_;


    // Private Member Functions

        //- Initialise
        void initialise(const dictionary&);

        //- Disallow default bitwise copy construct
        domSolar(const domSolar&);

        //- Disallow default bitwise assignment
        void operator=(const domSolar&);


public:

    //- Runtime type information
    TypeName("domSolar");


    // Constructors

        //- Construct from volScalarField
        domSolar(const volScalarField& T);

        //- Construct from dictionary and volScalarField
        domSolar(const dictionary& dict, const volScalarField& T);

        //- Constructor from local components. Does not create a radiationModel.
        //  radWallFieldName is the solar heat field name
        domSolar
        (
            const dictionary& dict,
            const volScalarField& T,
            const word radWallFieldName
        );


    //- Destructor
    virtual ~domSolar();


    // Member functions

        // Edit

            //- Solve
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh>> Ru() const;


    // Access

            const radiantField& operator[](const label i) const
            {
                return IenvPtr_[i];
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

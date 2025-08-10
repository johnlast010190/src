/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM® : Professional Open-source CFD
|   o   O   o    |
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------

License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM® <http://www.openfoam.org/>.

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
    © 2018 ESI Ltd.

Author
    2018. Nikolaos Magoulas (Esi Ltd.). All rights reserved.

InClass
    Foam::fv::artificialViscosityUSource

Description
    Artificial viscosity source for  momentum equation. A sensor field is
    needed to provide an input field, which used to calculate an artificial
    viscosity field, which is then used to add a diffusion term in the
    velocity equation.

    \heading Source usage

    Example usage:
    \verbatim
    artVisc
    {
        type    artificialViscosityUSource;
        active  yes;
        artificialViscosityUSourceCoeffs
        {
            selectionMode   all;
            fields          (U);
        }
        viscosityCoeff  1e-7;
        writeViscosity  true;
        maxViscosity    2e+2;

        sensor
        {
            type        maxSingularValue;
            fieldName   U;
        }
    }
    \endverbatim

SourceFiles
    artificialViscosityUSource.C

\*---------------------------------------------------------------------------*/

#ifndef artificialViscosityUSource_H
#define artificialViscosityUSource_H

#include "cellSetOption/cellSetOption.H"
#include "cfdTools/general/sensor/sensor/sensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{


/*---------------------------------------------------------------------------*\
                   Class artificialViscosityUSource Declaration
\*---------------------------------------------------------------------------*/

class artificialViscosityUSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Sensor function to provide the input field for the artificial
        //  viscosity
        autoPtr<sensor<vector>> sensor_;

        //- Multiplication factor for the sensor input field
        scalar viscosityCoeff_;

        //- If enabled writes the viscosity field
        Switch writeViscosity_;

        //- Clip the artificial viscosity value
        scalar maxViscosity_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        artificialViscosityUSource(const artificialViscosityUSource&);

        //- Disallow default bitwise assignment
        void operator=(const artificialViscosityUSource&);


public:

    //- Runtime type information
    TypeName("artificialViscosityUSource");


    // Constructors

        //- Construct from components
        artificialViscosityUSource
        (
            const word& name,
            const word& modelType,
            const dictionary& optionDict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~artificialViscosityUSource()
    {}


    // Member Functions

        // Application of artificial viscosity

            //- Add diffusion source term in equation
            virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);

            //- Add diffusion source term in equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvBlockMatrix<vector>& eqn,
                const label fieldI
            );

        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

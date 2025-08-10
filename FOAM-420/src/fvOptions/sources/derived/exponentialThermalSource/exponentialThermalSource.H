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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::fv::exponentialThermalSource

Description
    Exponential thermal source
    Q = Cm_*(1-exp(-magU/Ce))*(T0 - T) //W/m3


    \heading Source usage

    Example usage:
    \verbatim
    exponentialThermalSourceCoeffs
    {
        Cm 0.291e6;
    Ce 1.369;
        UName U;
        T0 353;
    eName h; // or T for ico
    }
    \endverbatim


SourceFiles
    exponentialThermalSource.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialThermalSource_H
#define exponentialThermalSource_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class exponentialThermalSource Declaration
\*---------------------------------------------------------------------------*/

class exponentialThermalSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Energy field name, default = h
        word eName_;

        //- Velocity field name, default = U
        word UName_;

        //- radiator temperature
        scalar T0_;

        //- multiplicative constant
        scalar Cm_;

        //- exponent constant
        scalar Ce_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        exponentialThermalSource(const exponentialThermalSource&);

        //- Disallow default bitwise assignment
        void operator=(const exponentialThermalSource&);


public:

    //- Runtime type information
    TypeName("exponentialThermalSource");


    // Constructors

        //- Construct from components
        exponentialThermalSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~exponentialThermalSource()
    {}


    // Member Functions

        //- Initialise after all fvOptions created
        virtual bool initialise() override;


        // Add explicit and implicit contributions

            //- Vector
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            ) override;

            //- Add contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldI
            ) override;

        // I-O

            //- Write data
            virtual void writeData(Ostream&) const override;

            //- Read dictionary
            virtual bool read(const dictionary& dict) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2021 Esi Ltd.

Class
    Foam::fv::scalarFieldLimiterSource

Description
    scalarField limiter source, used to suppress unphysical scalarField by
    either adding resistance to the scalarField matrix and/or clipping
    the scalarField. Usage is controled through the 'scalarFieldClip'
    and 'scalarFieldDamp' switches.

    \heading Source usage

    Example usage:
    \verbatim
    scalarFieldLimiterSourceCoeffs
    {
        limit 1;

        //optional (active settings are defaults)
        scalarFieldDamp false; //controls injection of resistance term
        scalarFieldClip true; //controls clipping of excess scalarField
        verbose false; //controls verbosity of limiter output

        fieldName alpha.air; // name of scalarField field to limit (no default)
        //fieldNames (alpha.air alpha.water);
    }
    \endverbatim

SourceFiles
    scalarFieldLimiterSource.C

\*---------------------------------------------------------------------------*/

#ifndef scalarFieldLimiterSource_H
#define scalarFieldLimiterSource_H

#include "cellSetOption/cellSetOption.H"
#include "primitives/bools/Switch/Switch.H"
#include "primitives/functions/Polynomial/Polynomial.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{


/*---------------------------------------------------------------------------*\
                   Class scalarFieldLimiterSource Declaration
\*---------------------------------------------------------------------------*/

class scalarFieldLimiterSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- scalarField at which to start damping
        scalar limit_;

        //- Switches

            //- switch for clipping excess scalarField
            Switch scalarFieldClip_;

            //- switch for clipping excess scalarField
            Switch scalarFieldDamp_;

            //- switch controlling output verbosity
            Switch verbose_;

            //- switch to control output of diagnostic fields
            Switch outputDiagnostics_;

    // Protected Member Functions

        //- damping calculator
        void addDamping(fvMatrix<scalar>& eqn);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        scalarFieldLimiterSource(const scalarFieldLimiterSource&);

        //- Disallow default bitwise assignment
        void operator=(const scalarFieldLimiterSource&);


public:

    //- Runtime type information
    TypeName("scalarFieldLimiterSource");


    // Constructors

        //- Construct from components
        scalarFieldLimiterSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~scalarFieldLimiterSource()
    {}


    // Member Functions

        // Clip values

            //- Scalar
            virtual void correct
            (
                volScalarField& vf
            );

            //- Constrain scalar matrix
            virtual void constrain(fvMatrix<scalar>& eqn, const label fieldI);

        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

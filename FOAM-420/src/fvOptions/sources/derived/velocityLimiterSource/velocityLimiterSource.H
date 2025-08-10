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
    (c) 2022 Esi Ltd.

Class
    Foam::fv::velocityLimiterSource

Description
    Velocity limiter source, used to suppress unphysical velocity by
    either adding resistance to the velocity matrix and/or clipping
    the velocity field. Usage is controled through the 'velocityClip'
    and 'velocityDamp' switches.

    \heading Source usage

    Example usage:
    \verbatim
    velocityLimiterSourceCoeffs
    {
        Ulimit 100;

        //optional (active settings are defaults)
        velocityDamp false; //controls injection of resistance term
        velocityClip true; //controls clipping of excess velocity
        verbose false; //controls verbosity of limiter output

        UName U; // default, name of velocity field to limit
        //UNames (U.air U.water);
    }
    \endverbatim

SourceFiles
    velocityLimiterSource.C

\*---------------------------------------------------------------------------*/

#ifndef velocityLimiterSource_H
#define velocityLimiterSource_H

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
                   Class velocityLimiterSource Declaration
\*---------------------------------------------------------------------------*/

class velocityLimiterSource
:
    public cellSetOption
{

    // Private data

        //- velocity at which to start damping
        scalar Ulimit_;

        //- Switches

            //- switch for clipping excess velocity
            Switch velocityClip_;

            //- switch for clipping excess velocity
            Switch velocityDamp_;

            //- switch controlling output verbosity
            Switch verbose_;

            //- switch to control output of diagnostic fields
            Switch outputDiagnostics_;

    // Private Member Functions

        //- damping calculator
        void addDamping(fvMatrix<vector>& eqn);

        //- damping calculator
        void addDamping(fvBlockMatrix<vector>& eqn);

        //- printing damping information
        void printDampingInformation
        (
            const label&,
            const volVectorField&,
            const scalarField&
        ) const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        velocityLimiterSource(const velocityLimiterSource&);

        //- Disallow default bitwise assignment
        void operator=(const velocityLimiterSource&);


public:

    //- Runtime type information
    TypeName("velocityLimiterSource");


    // Constructors

        //- Construct from components
        velocityLimiterSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~velocityLimiterSource()
    {}


    // Member Functions

        // Clip values

            //- Vector
            virtual void correct
            (
                volVectorField& U
            );

            //- Constrain vector matrix
            virtual void constrain(fvMatrix<vector>& eqn, const label fieldI);

            //- Constrain block vector matrix
            virtual void constrain
            (
                fvBlockMatrix<vector>& eqn,
                const label fieldI
            );

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

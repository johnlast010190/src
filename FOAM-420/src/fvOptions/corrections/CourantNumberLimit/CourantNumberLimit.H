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
    (c) 2017 Esi Ltd.
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

Class
    Foam::fv::CourantNumberLimit

Group
    grpFvOptionsConstraints

Description
    Constraint to clip velocity extrema in transient calculations

    Velocity is limited based on Courant number. No source or sinks are
    applied, excessive velocity is simply scaled to comply with the
    constraint CFL.

    Constraint described by:

        CourantNumberLimitCoeffs
        {
            cflMax        100;

            // Optional: name of velocity field (default: U)
            //UName       U;

            //Optional flux name
            //phiName phi;

            //Optional phase field name
            //phaseName none;

            //Optional add scaled clipped velocity back into field
            //residualScale 0;

        }


SourceFiles
    CourantNumberLimit.C

\*---------------------------------------------------------------------------*/

#ifndef CourantNumberLimit_H
#define CourantNumberLimit_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class CourantNumberLimit Declaration
\*---------------------------------------------------------------------------*/

class CourantNumberLimit
:
    public cellSetOption
{

protected:

    // Protected data

        //- Maximum Courant number
        scalar cflMax_;

        //- scale clipped portion of velocity, useful for diagnostics
        scalar residualScale_;

        //- Flux field used to calculate Co number
        word phiName_;

        //- multiphase field for interface CFL limiting
        word phaseName_;

        //- density name
        word rhoName_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        CourantNumberLimit(const CourantNumberLimit&);

        //- Disallow default bitwise assignment
        void operator=(const CourantNumberLimit&);


public:

    //- Runtime type information
    TypeName("CourantNumberLimit");


    // Constructors

        //- Construct from components
        CourantNumberLimit
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~CourantNumberLimit()
    {}


    // Member Functions

        // Clip values
            //- Vector
            virtual void correct
            (
                volVectorField& U
            );


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

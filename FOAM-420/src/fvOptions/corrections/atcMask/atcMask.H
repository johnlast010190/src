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
    Foam::fv::atcMask

Description
    Masks the ATC field according to a sensor field.

Usage
    Example usage:
    \verbatim
    mask
    {
        type            atcMask;
        active          yes;

        atcMaskCoeffs
        {
            selectionMode   all;
            fieldName       atcMask_drag
            max             1e10;
        }

        sensor
        {
            type        maxSingularValue;
            fieldName   U;
        }
    }
    \endverbatim

SourceFiles
    atcMask.C

\*---------------------------------------------------------------------------*/

#ifndef atcMask_H
#define atcMask_H

#include "cellSetOption/cellSetOption.H"
#include "cfdTools/general/sensor/sensor/sensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class atcMask Declaration
\*---------------------------------------------------------------------------*/

class atcMask
:
    public cellSetOption
{

protected:

    // Protected data

        //- Name of field to be masked
        word fieldName_;

        //- Sensor function to be used as input for the masking
        autoPtr<sensor<vector>> sensor_;

        //- Minimum value of sensor field for which masking doesn't an effect
        scalar min_;

        //- Maximum value of sensor field after which masking is fully applied
        scalar max_;

        //- Write mask
        Switch writeMask_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        atcMask(const atcMask&);

        //- Disallow default bitwise assignment
        void operator=(const atcMask&);


public:

    //- Runtime type information
    TypeName("atcMask");


    // Constructors

        //- Construct from components
        atcMask
        (
            const word& name,
            const word& modelType,
            const dictionary& optionDict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~atcMask()
    {}


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Correct the energy field
        virtual void correct(volScalarField& atcMask);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

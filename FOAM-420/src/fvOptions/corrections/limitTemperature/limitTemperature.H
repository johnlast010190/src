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
    (c) 2012-2017 OpenFOAM Foundation

Class
    Foam::fv::limitTemperature

Group
    grpFvOptionsCorrections

Description
    Correction for temperature to apply limits between minimum and maximum
    values.

Usage
    Example usage:
    \verbatim
    limitT
    {
        type            limitTemperature;
        active          yes;
        field           T;
        selectionMode   all;
        min             200;
        max             500;
    }
    \endverbatim

SourceFiles
    limitTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef limitTemperature_H
#define limitTemperature_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class limitTemperature Declaration
\*---------------------------------------------------------------------------*/

class limitTemperature
:
    public cellSetOption
{

protected:

    // Protected data

        //- Energy eqn. field name, default = h
        word TName_;

        //- Minimum temperature limit [K]
        scalar Tmin_;

        //- Maximum temperature limit [K]
        scalar Tmax_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        limitTemperature(const limitTemperature&);

        //- Disallow default bitwise assignment
        void operator=(const limitTemperature&);


public:

    //- Runtime type information
    TypeName("limitTemperature");


    // Constructors

        //- Construct from components
        limitTemperature
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~limitTemperature()
    {}


    // Member Functions

        //- Initialise after all fvOptions created
        virtual bool initialise() override;

        //- Read dictionary
        virtual bool read(const dictionary& dict) override;

        //- Correct the energy field
        virtual void correct(volScalarField& he) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

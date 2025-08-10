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
    (c) 2016-2017 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

Class
    Foam::fv::EquationFlux

Group
    grpFvOptionsPostProcessing

Description
    Write out the boundary flux density for an equation.

Usage
    For example to write the boundary heat fluxes:
    \verbatim
    EquationFlux
    {
        type            scalarEquationFlux;
        active          yes;
        equation        h;
        writeField      yes; // Defaults to 'no' if outputFieldName not specified, otherwise 'yes'
        outputFieldName wallHeatFlux; // Defaults to <equation>EquationFlux
    }
    \endverbatim

See also
    Foam::fvOption

SourceFiles
    EquationFlux.C
    EquationFluxes.C

\*---------------------------------------------------------------------------*/

#ifndef EquationFlux_H
#define EquationFlux_H

#include "cellSetOption/cellSetOption.H"
#include "db/functionObjects/writeFile/writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class EquationFlux Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class EquationFlux
:
    public option,
    public functionObjects::writeFile
{
    // Private member data

        //- Equation name
        word equation_;

        //- Output field name
        word outputFieldName_;

        //- Whether to write the field to disk
        Switch writeField_;

        //- Whether to store the field in object registry
        Switch storeField_;

        //- Whether header written yet
        bool headerWritten;

protected:

    // Protected member functions

        void writeFileHeader(Ostream& os, const fvMatrix<Type>& eqn) const;

public:

    //- Runtime type information
    TypeName("EquationFlux");


    // Constructors

        //- Construct from components
        EquationFlux
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        //- Grab flux from assembled equation
        virtual void constrain(fvMatrix<Type>& eqn, const label fieldi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EquationFlux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "vtk/format/foamVtkLegacyAsciiFormatter.H"
#include "vtk/output/foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::legacyAsciiFormatter::legacyName_ = "ASCII";

const Foam::vtk::outputOptions
Foam::vtk::legacyAsciiFormatter::opts_(formatType::LEGACY_ASCII);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::legacyAsciiFormatter::legacyAsciiFormatter
(
    std::ostream& os
)
:
    asciiFormatter(os)
{}


Foam::vtk::legacyAsciiFormatter::legacyAsciiFormatter
(
    std::ostream& os,
    unsigned precision
)
:
    asciiFormatter(os, precision)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::vtk::outputOptions&
Foam::vtk::legacyAsciiFormatter::opts() const
{
    return opts_;
}


const char* Foam::vtk::legacyAsciiFormatter::name() const
{
    return legacyName_;
}


const char* Foam::vtk::legacyAsciiFormatter::encoding() const
{
    return legacyName_;
}


// ************************************************************************* //

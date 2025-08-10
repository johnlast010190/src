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
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "vtk/format/foamVtkBase64Formatter.H"
#include "vtk/output/foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::base64Formatter::name_ = "binary";

const Foam::vtk::outputOptions
Foam::vtk::base64Formatter::opts_(formatType::INLINE_BASE64);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::base64Formatter::base64Formatter(std::ostream& os)
:
    foamVtkBase64Layer(os)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::base64Formatter::~base64Formatter()
{
    if (base64Layer::close())
    {
        os().put('\n');
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::vtk::outputOptions&
Foam::vtk::base64Formatter::opts() const
{
    return opts_;
}


const char* Foam::vtk::base64Formatter::name() const
{
    return name_;
}


void Foam::vtk::base64Formatter::flush()
{
    if (base64Layer::close())
    {
        os().put('\n');
    }
}


// ************************************************************************* //

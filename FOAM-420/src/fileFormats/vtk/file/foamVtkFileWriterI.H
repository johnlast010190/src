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
    (c) 2018 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline std::ofstream& Foam::vtk::fileWriter::os()
{
    return os_;
}


inline Foam::vtk::formatter& Foam::vtk::fileWriter::format()
{
    return *format_;
}


inline bool Foam::vtk::fileWriter::isState(outputState test) const
{
    return (test == state_);
}


inline bool Foam::vtk::fileWriter::notState(outputState test) const
{
    return (test != state_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::vtk::fileTag Foam::vtk::fileWriter::contentType() const
{
    return contentType_;
}


inline Foam::vtk::outputOptions Foam::vtk::fileWriter::opts() const
{
    return opts_;
}


inline Foam::word Foam::vtk::fileWriter::ext() const
{
    return opts_.ext(contentType_);
}


inline bool Foam::vtk::fileWriter::legacy() const
{
    return opts_.legacy();
}


inline bool Foam::vtk::fileWriter::parallel() const
{
    return parallel_;
}


inline const Foam::word& Foam::vtk::fileWriter::state() const
{
    return stateNames[state_];
}


inline const Foam::fileName& Foam::vtk::fileWriter::output() const
{
    return outputFile_;
}


inline Foam::label Foam::vtk::fileWriter::nCellData() const
{
    return nCellData_;
}


inline Foam::label Foam::vtk::fileWriter::nPointData() const
{
    return nPointData_;
}


// ************************************************************************* //

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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

if (writeVacf)
{
    OFstream vacfFile(runTime.path()/"vacf");

    if (!vacf.writeAveraged(vacfFile))
    {
        FatalErrorInFunction
            << "Failed writing to "
            << vacfFile.name()
            << abort(FatalError);
    }
}

Info<< "Diffusion coefficient = "
    << vacf.integral() << endl;

if (writePacf)
{
    OFstream pacfFile(runTime.path()/"pacf");

    if (!pacf.writeAveraged(pacfFile))
    {
        FatalErrorInFunction
            << "Failed writing to "
            << pacfFile.name()
            << abort(FatalError);
    }
}

Info<< "Viscosity = "
    << pacf.integral()/averageTemperature/moleculeCloud::kb/meshVolume
    << endl;

if (writeHFacf)
{
    OFstream hfacfFile
    (
        runTime.path()/ + "hfacf"
    );

    if (!hfacf.writeAveraged(hfacfFile))
    {
        FatalErrorInFunction
            << "Failed writing to "
            << hfacfFile.name()
            << abort(FatalError);
    }
}

Info<< "Thermal conductivity = "
    << hfacf.integral()
        /averageTemperature
        /averageTemperature
        /moleculeCloud::kb
        / meshVolume
    << endl;
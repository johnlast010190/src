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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

Info<< nl << "Creating autocorrelation functions." << endl;

IOdictionary mdTransportProperitesDict
(
    IOobject
    (
        "mdTransportProperitesDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    )
);

const dictionary& autocorrelationFunctionDict
(
    mdTransportProperitesDict.subDict("autocorrelationFunctions")
);

//- Velocity autocorrelation function

Info<< tab << "velocty" << endl;

const dictionary& velocityACFDict
(
    autocorrelationFunctionDict.subDict("velocity")
);

correlationFunction<vector> vacf
(
    mesh,
    velocityACFDict,
    molecules.size()
);

bool writeVacf(Switch(velocityACFDict.lookup("writeFile")));

//- Pressure autocorrelation function

Info<< tab << "pressure" << endl;

const dictionary& pressureACFDict
(
    autocorrelationFunctionDict.subDict("pressure")
);

correlationFunction<vector> pacf
(
    mesh,
    pressureACFDict,
    1
);

bool writePacf(Switch(pressureACFDict.lookup("writeFile")));

//- Heat flux autocorrelation function

Info<< tab << "heat flux" << endl;

const dictionary& heatFluxACFDict
(
    autocorrelationFunctionDict.subDict("heatFlux")
);

correlationFunction<vector> hfacf
(
    mesh,
    heatFluxACFDict,
    1
);

bool writeHFacf(Switch(heatFluxACFDict.lookup("writeFile")));

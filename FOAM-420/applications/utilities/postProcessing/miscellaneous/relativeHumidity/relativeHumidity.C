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
    (c) 2010-2011, Esi Ltd

Application
    ptot

Description
    For each time: calculate the total pressure.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "rhoThermo/rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "include/createMesh.H"

    dimensionedScalar Mvap("Mwv", dimensionSet(1, 0, 0, 0, 1), 18.02e-3);
    dimensionedScalar Mair("Mwv", dimensionSet(1, 0, 0, 0, 1), 28.96e-3);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Pheader
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject Theader
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject wheader
        (
            "w",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );


        // Check p and U exist
        if
        (
            Pheader.typeHeaderOk<volScalarField>(true)
            && Theader.typeHeaderOk<volScalarField>(true)
            && wheader.typeHeaderOk<volScalarField>(true)
        )
        {
            mesh.readUpdate();

            Info<< "    Reading p" << endl;
            volScalarField p(Pheader, mesh);
            Info<< "    Reading T" << endl;
            volScalarField T(Theader, mesh);
            Info<< "    Reading w" << endl;
            volScalarField w(wheader, mesh);

            if (p.dimensions() == dimPressure)
            {
                rhoThermo& thermo =
                    refCast<rhoThermo>
                    (
                        basicThermo::lookupOrCreate(mesh.thisDb())
                    );
                p += thermo.pRef();
            }
            else
            {
                IOdictionary transportProperties
                (
                    IOobject
                    {
                        "transportProperties",
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    }
                );
                dimensionedScalar rho
                (
                    "rho",
                    transportProperties.lookup("rho")
                );
                dimensionedScalar pRef
                (
                    "pRef",
                    transportProperties.lookup("pRef")
                );
                p *= rho;
                p += pRef;
            }

            Info<< "    Calculating relative humidity" << endl;

            volScalarField Psat
            (
                IOobject
                (
                    "Psat",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                mesh,
                dimensionedScalar("Psat", dimPressure, 0.0)
            );

            Psat.primitiveFieldRef() =
                2337*Foam::exp(6879*(1/293.15 - 1/T.primitiveField())
              - 5.031*Foam::log(T.primitiveField()/293.15));

            forAll(Psat.boundaryField(), patchi)
            {
                Psat.boundaryFieldRef()[patchi] =
                    2337*Foam::exp
                    (
                        6879*(1/293.15 - 1/T.boundaryField()[patchi]
                    )
                  - 5.031*Foam::log(T.boundaryField()[patchi]/293.15));
            }

            volScalarField relativeHumidity
            (
                IOobject
                (
                    "relativeHumidity",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                w*p*Mair/(Psat*(Mvap + w*(Mair - Mvap)))
            );
            relativeHumidity.write();

        }
        else
        {
            Info<< "    No p, T or w" << endl;
        }

        Info<< endl;
    }

    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //

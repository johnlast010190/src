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

Description
    Compare two scalar or vector fields, write the difference as a new field.
    Also write the overall difference to the screen.

\*---------------------------------------------------------------------------*/
#include "cfdTools/general/include/fvCFD.H"
#include "meshes/polyMesh/polyPatches/derived/baseWall/wall.H"


int main(int argc, char *argv[])
{

	argList:: addNote("compareFields compares two scalar fields or vector fields, write \
the difference into a new field for visualization. It also print the overall \
relative difference on screen.");
	argList::noParallel();

    argList::addOption
    (
        "field",
        "p",
        "Specify the field name to be compared"
    );

    argList::addOption
    (
        "baseDir",
        "Specify the directory for the base case"
    );

    argList::addOption
    (
        "isVector",
        "false",
        "Specify whether the two fields are vectors"
    );

    argList::addOption
    (
        "timeName",
        "Specify the time name for the field to compare"
    );

    DynamicList<word> argNames;
	DynamicList<word> argValues;

	if (argc>1)
	{
		for (label i=1;i<argc-1;i+=2)
		{
			argNames.append(argv[i]);
			argValues.append(argv[i+1]);
		}
	}

	#include "include/setRootCase.H"
    #include "include/createTime.H"

	word fieldName="p";
	word baseDir="./base";

	bool isVector=false;
	word timeName=runTime.timeName();

	if (argNames.size()>=1)
	{
		for (label i=0;i<argNames.size();i++)
		{
			word key=argNames[i];
			word value=argValues[i];

			if (key=="-field")
			{
				fieldName=value;
				continue;
			}

			if (key=="-baseDir")
			{
				baseDir=value;
				continue;
			}

			if (key=="-isVector" && value=="true")
			{
				isVector=true;
				continue;
			}

			if (key=="-timeName")
			{
				timeName=value;
				continue;
			}
		}
	}

	scalar ftime=atof(timeName.c_str());
	runTime.setTime(ftime,0);

	#include "include/createMesh.H"

	scalar relerror=0;
	word diffName=fieldName+"_diff";

	word baseField=fieldName+"_base";

	cp(baseDir+"/"+timeName+"/"+fieldName, "./"+timeName+"/"+baseField);

	if (!isVector)
	{
		volScalarField sc
		(
			IOobject
			(
				fieldName,
				timeName,
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);

		volScalarField sc0
		(
			IOobject
			(
				baseField,
				timeName,
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);


		volScalarField scDiff
		(
			IOobject
			(
				diffName,
				timeName,
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			sc*0.0
		);

		scDiff=sc-sc0;
		scalar fldSum=0;
		scalar volSum=0;

		forAll(mesh.V(),i)
		{
			fldSum+=fabs(scDiff.primitiveField()[i])*mesh.V()[i];
			volSum+=fabs(sc0.primitiveField()[i])*mesh.V()[i];
		}

		relerror=fldSum/volSum;

		Info<< "Time now = " << runTime.timeName() << endl;
		runTime.writeNow();
	}

	else
	{
		volVectorField vel
		(
			IOobject
			(
				fieldName,
				timeName,
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);

		volVectorField vel0
		(
			IOobject
			(
				baseField,
				timeName,
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);

		volVectorField velDiff
		(
			IOobject
			(
				diffName,
				timeName,
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			vel*0.0
		);

		velDiff=vel-vel0;

		scalar fldSum=0;
		scalar volSum=0;

		forAll(mesh.V(),i)
		{
			fldSum+=mag(velDiff.primitiveField()[i])*mesh.V()[i];
			volSum+=mag(vel0.primitiveField()[i])*mesh.V()[i];
		}

		relerror=fldSum/volSum;

		Info<< "Time now = " << runTime.timeName() << endl;
		runTime.writeNow();

	}


	Info<<"Over-all relative difference of the two fields "<<fieldName<<" and "
	<<baseField<<" is:"<<relerror<<endl;

	Info<<"runTime: "<<runTime.timeName()<<endl;
	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;

    Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //

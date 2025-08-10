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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2017 CSIR
    (c) 2020 Esi Ltd.

Class
    Foam::functionObjects::wallHeatFlux

Group
    grpForcesFunctionObjects

Description
    Calculates and writes the heat flux at wall patches as the
    volScalarField field 'wallHeatFlux', and calculates the min/max and
    heat flux at patches.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

Usage
    Example of function object specification:
    \verbatim
    wallHeatFlux1
    {
        type        wallHeatFlux;
        ...
        patches     (".*Wall");
        kappaMethod lookup;
        kappa       k;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallHeatFlux    | yes        |
        patches  | list of patches to process | no         | all wall patches
        kappaMethod  | Thermal conductivity method        | no  | fluidThermo
        kappa        | Name of thermal conductivity field | no  | none
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    wallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallHeatFlux_H
#define functionObjects_wallHeatFlux_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "fields/volFields/volFieldsFwd.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "fvMesh/fvMesh.H"
#include "turbulentFluidThermoModels/derivedFvPatchFields/boundaryKappa/boundaryKappa.H"
#include "turbulentTransportModels/turbulentTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFlux
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

//        const fvMesh& mesh_;

        //- Optional list of patches to process
        labelHashSet patchSet_;

        //- Optional phase (group) name
        word phaseName_;

        PtrList<boundaryKappa> boundaryKappa_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the heat-flux
        void calcHeatFlux
        (
            const volScalarField& T,
            volScalarField& wallHeatFlux
        );

        //- Calculate the heat-flux (incompressible special-case)
        void calcHeatFluxIco
        (
            const incompressible::turbulenceModel& turb,
            const volScalarField& T,
            volScalarField& wallHeatFlux
        );


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        wallHeatFlux(const wallHeatFlux&);

        //- Disallow default bitwise assignment
        void operator=(const wallHeatFlux&);


public:

    //- Runtime type information
    TypeName("wallHeatFlux");


    // Constructors

        //- Construct from Time and dictionary
        wallHeatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~wallHeatFlux();


    // Member Functions

        //- Read the wallHeatFlux data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

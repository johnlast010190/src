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
    (c) 2022 Esi Ltd.

Authors
    2019. Nikolaos Magoulas (Esi Ltd.). All rights reserved.

Class
    Foam::objectiveFunctionObject

Description
    Abstract base class for function objects of objectives.

SourceFiles
    objectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef OBJECTIVEFUNCTIONOBJECT_H
#define OBJECTIVEFUNCTIONOBJECT_H

#include "cfdTools/general/include/fvCFD.H"
#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class objectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class objectiveFunctionObject
:
    public functionObjects::fvMeshFunctionObject,
    public functionObjects::writeFile
{

protected:

    // Protected member data

        //- Name of velocity field
        word UName_;

        //- Name of pressure field
        word pName_;

        //- Name of pressure field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Name of species field
        word psiName_;

        //- Name of temperature field
        word TName_;

        //- Value of objective
        scalar objectiveValue_;

        //- List of boundaries this objective applies to (if any)
        boolList objectivePatch_;

        //- Name of the base geometry
        word baseGeometryName_;

        //- list of names of MRF's to optimize
        // if list has zero size then all the zones will be optimized
        wordList zoneNames_;

        //- Pointer to the inletPatches
        autoPtr<boolList> inletPatchPtr_;

        //- Write objective values to file
        bool logToFile_;

        //- Pointer to file in which obj values are written
        autoPtr<OFstream> objFilePtr_;


    // Protected member functions

        // Access

            //- Return const reference to the velocity field
            inline const volVectorField& U() const
            {
                return lookupObject<volVectorField>(UName_);
            }

            //- Return const reference to the pressure field
            inline const volScalarField& p() const
            {
                return lookupObject<volScalarField>(pName_);
            }

            //- Return const reference to the phi field
            inline const surfaceScalarField& phi() const
            {
                return lookupObject<surfaceScalarField>(phiName_);
            }

            //- Return const reference to the species field
            inline const volScalarField& psi() const
            {
                return lookupObject<volScalarField>(psiName_);
            }

            //- Return const reference to the temperature field
            inline const volScalarField& temperature() const
            {
                return lookupObject<volScalarField>(TName_);
            }

            //- Return the real pressure field
            tmp<volScalarField> P() const;

            //- Return the density field
            tmp<volScalarField> rho() const;

            //- Return the laminar viscosity field
            tmp<volScalarField> mu() const;

            //- Return the effective viscosity field
            tmp<volScalarField> muEff() const;

            //- Return the devRhoReff field
            tmp<volSymmTensorField> devRhoReff() const;


        // Process

            //- Select objective patches
            boolList getObjectivePatches(const dictionary&) const;

            //- Calculate total area of objective inlet patches
            scalar objectivePatchInletArea() const;

            //- Calculate total area of objective outlet patches
            scalar objectivePatchOutletArea() const;

            //- Calculate the flow rate of the inlets
            scalar inletFlowRate() const;

            //- Calculate the power loss
            scalar calculatePowerLoss() const;


        // Write

            //- Write the header of the output file in postProcessing
            virtual void writeFileHeader(Ostream&);

            //- write objective values to log
            void logToFile();

            //- Create log files in postProcessing folder
            void createFiles(bool);


private:

    // Private Member Functions

        //- Select inlet outlet patches (used for subdomain solution)
        void setInletPatches(const dictionary&);

        //- Disallow default bitwise copy construct and assignment
        objectiveFunctionObject(const objectiveFunctionObject&);
        void operator=(const objectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("objectiveFunctionObject");


    // Constructors

        //- Construct from components for composite
        objectiveFunctionObject
        (
            const Time&,
            dictionary&,
            bool useAdjointFileFormat = false
        );


    //- Destructor
    virtual ~objectiveFunctionObject();


    // Member Functions

        // Access

            //- Return objective value
            inline scalar objectiveValue() const
            {
                return objectiveValue_;
            }

            //- Return activity switch for specified boundary index
            inline bool objectivePatch(label i) const
            {
                return objectivePatch_[i];
            };

            //- Return true if compressible case
            bool compressible() const;

            //- Read the field average data
            virtual bool read(const dictionary&) = 0;

            //- Execute the averaging
            virtual bool execute() = 0;

            //- Track and write
            virtual bool write() = 0;

            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            virtual void writeOutputToFile();

            virtual bool writeOutputToAdjointFile() const;

        // Modify

            //- Set the name of the base geometry
            void setBaseGeometryName(const word& geometryName);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

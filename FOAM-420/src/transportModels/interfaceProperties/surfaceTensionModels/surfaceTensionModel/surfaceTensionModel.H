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
    (c) 2017 OpenFOAM Foundation

Class
    Foam::surfaceTensionModel

Description
    Abstract base-class for surface tension models which return the surface
    tension coefficient field.

Usage
    Example of the surface tension specification:
    \verbatim
        sigma
        {
            type                <surface tension model type>;
            <coefficient name>  <coefficient value>;
            .
            .
            .
        }
    \endverbatim
    For simplicity and backward-compatibility the constant value format is
    also supported, e.g.
    \verbatim
        sigma           0.07;
    \endverbatim

SourceFiles
    surfaceTensionModel.C
    newSurfaceTensionModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTensionModel_H
#define surfaceTensionModel_H

#include "db/regIOobject/regIOobject.H"
#include "dimensionedTypes/dimensionedTypes.H"
#include "fields/volFields/volFieldsFwd.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class surfaceTensionModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionModel
:
    public regIOobject
{
protected:

    // Protected member data

        //- Reference to mesh
        const fvMesh& mesh_;


    // Protected member functions

        static const dictionary& sigmaDict(const dictionary& dict)
        {
            return dict.subDict("sigma");
        }


public:

    //- Runtime type information
    TypeName("surfaceTensionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceTensionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Static data members

        //- Surface tension coefficient dimensions
        static const dimensionSet dimSigma;


    // Constructors

        // Construct from mesh
        surfaceTensionModel(const fvMesh& mesh);


    //- Destructor
    virtual ~surfaceTensionModel();


    // Selectors

        static autoPtr<surfaceTensionModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Surface tension coefficient
        virtual tmp<volScalarField> sigma() const = 0;

        //- Update surface tension coefficient from given dictionary
        virtual bool readDict(const dictionary& dict) = 0;

        //- Write in dictionary format
        virtual bool writeData(Ostream& os) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

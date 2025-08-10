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
    (c) 2010-2016 Esi Ltd.
    (c) 2012-2015 OpenFOAM Foundation

Class
    Foam::porosityModel

Description
    Top level model for porosity models

SourceFiles
    porosityModel.C
    porosityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef porosityModel_H
#define porosityModel_H

#include "fvMesh/fvMesh.H"
#include "db/dictionary/dictionary.H"
#include "fvMatrices/fvMatricesFwd.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "coordinate/systems/coordinateSystem.H"
#include "dimensionedTypes/dimensionedVector/dimensionedVector.H"
#include "primitives/strings/keyType/keyType.H"
#include "finiteVolume/fvc/fvc.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class porosityModel Declaration
\*---------------------------------------------------------------------------*/

class porosityModel
:
    public regIOobject
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        porosityModel(const porosityModel&);

        //- Disallow default bitwise assignment
        void operator=(const porosityModel&);


protected:

    // Protected data

        //- Porosity name
        word name_;

        //- Reference to the solution registry
        const objectRegistry& obr_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Dictionary used for model construction
        const dictionary dict_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Porosity active flag
        bool active_;

        //- Use relative Velocity
        bool porousRelVelocity_;

        //- Name(s) of cell-zone
        keyType zoneName_;

        //- Cell zone IDs
        labelList cellZoneIDs_;

        //- Reference Frame Ptr
        coordinateFrame* coorFramePtr_;


    // Protected Member Functions


        //- Transform the model data wrt mesh changes
        virtual void calcTransformModelData() = 0;

        //- Adjust negative resistance values to be multiplier of max value
        void adjustNegativeResistance(dimensionedVector& resist);

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const = 0;

        virtual void correct(fvVectorMatrix& UEqn) const = 0;
        virtual void correct(fvBlockMatrix<vector>&) const
        {
            NotImplemented;
        };

        virtual void correct
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const = 0;

        virtual void correct
        (
            fvBlockMatrix<vector>&,
            const volScalarField&,
            const volScalarField&
        ) const
        {
            NotImplemented;
        };

        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const = 0;


        virtual void adjointCorrect
        (
            fvVectorMatrix& UaEqn,
            const volVectorField& Uprimal
        ) const = 0;


        virtual void adjointCorrect
        (
            fvVectorMatrix& UaEqn,
            const volScalarField& rho,
            const volScalarField& mu,
            const volVectorField& Uprimal
        ) const = 0;



        virtual void adjointCorrect
        (
            fvBlockMatrix<vector>&,
            const volVectorField&
        ) const
        {
            NotImplemented;
        }


        virtual void adjointCorrect
        (
            fvBlockMatrix<vector>&,
            const volScalarField&,
            const volScalarField&,
            const volVectorField&
        ) const
        {
            NotImplemented;
        }

        //- Local coordinate system
        inline const coordinateSystem& csys() const;

        inline bool validParentFrame() const;

        //- Return label index
        label fieldIndex(const label index) const;

public:

    //- Runtime type information
    TypeName("porosityModel");

    //- Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        porosityModel,
        mesh,
        (
            const word& modelName,
            const word& name,
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cellZoneName
        ),
        (modelName, name, obr, mesh, dict, cellZoneName)
    );

    //- Constructor
    porosityModel
    (
        const word& name,
        const word& modelType,
        const objectRegistry& obr,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName = word::null
    );

    //- Return pointer to new porosityModel object created on the freestore
    //  from an Istream
    class iNew
    {
        const objectRegistry& obr_;
        const fvMesh& mesh_;
        const word& name_;

    public:

        iNew
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            const word& name
        )
        :
            obr_(obr),
            mesh_(mesh),
            name_(name)
        {}

        autoPtr<porosityModel> operator()(Istream& is) const
        {
            const dictionary dict(is);

            return autoPtr<porosityModel>
            (
                porosityModel::New
                (
                    name_,
                    obr_,
                    mesh_,
                    dict
                )
            );
        }
    };

    //- Selector
    static autoPtr<porosityModel> New
    (
        const word& name,
        const objectRegistry& obr,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName = word::null
    );

    //- Destructor
    virtual ~porosityModel();


    // Member Functions

        //- Return const access to the porosity model name
        inline const word& name() const;

        //- Return const access to the porosity active flag
        inline bool active() const;

        //- Return const access to the cell zone IDs
        inline const labelList& cellZoneIDs() const;

        //- Transform the model data wrt mesh changes
        virtual void transformModelData();

        //- Return the force over the cell zone(s)
        virtual tmp<vectorField> force
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu
        );

        //- Add resistance
        virtual void addResistance(fvVectorMatrix& UEqn);

        virtual void addResistance(fvBlockMatrix<vector>&);

        virtual void addResistance
        (
            fvBlockMatrix<vector>&,
            const volScalarField&,
            const volScalarField&
        );

        //- Add resistance
        virtual void addResistance
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        );

        //- Add resistance
        virtual void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC
        );

        //- Add adjoint resistance
        virtual void addAdjointResistance
        (
            fvVectorMatrix& UaEqn,
            const volVectorField& Uprimal
        );

        //- Add adjoint resistance
        virtual void addAdjointResistance
        (
            fvVectorMatrix& UaEqn,
            const volScalarField& rho,
            const volScalarField& mu,
            const volVectorField& Uprimal
        );

        //- Add resistance
        virtual void addAdjointResistance
        (
            fvBlockMatrix<vector>&,
            const volVectorField&
        );

        //- Add resistance
        virtual void addAdjointResistance
        (
            fvBlockMatrix<vector>&,
            const volScalarField&,
            const volScalarField&,
            const volVectorField&
        );

    // Access

        //- cellZone boundary faces
        labelList zoneBoundaryFaces() const;


    // I-O

        //- Write
        virtual bool writeData(Ostream& os) const;

        //- Inherit read from regIOobject
        using regIOobject::read;

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }

        //- Read porosity dictionary
        virtual bool read(const dictionary& dict);
};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<porosityModel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "porosityModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

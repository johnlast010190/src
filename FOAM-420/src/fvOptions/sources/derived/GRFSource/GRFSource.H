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
    (c) 2018 Esi Ltd.
    (c) 2012-2013 OpenFOAM Foundation

Class
    Foam::fv::GRFSource

Description

    Generalized Reference Frame (GRF) source

    Adds a rotational flux based source plus Coriolis. Relative/absolute fluxes
    are not used. The pressure equation is solved in the absolute frame.
    This class contains the information of all the GRFs used and combines the
    various fluxes accordingly.

    \heading Source usage

    Example usage:
    \verbatim

    \endverbatim


SourceFiles
    GRFSource.C

\*---------------------------------------------------------------------------*/
#ifndef GRFSource_H
#define GRFSource_H

#include "cellSetOption/cellSetOption.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "sources/derived/GRFSource/GRF/GRF.H"
#include "sources/derived/frameSources/frameSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class GRFSource Declaration
\*---------------------------------------------------------------------------*/

class GRFSource
:
    public option
{

protected:

    // Protected data

        //- Generalized Reference Frames
        List<autoPtr<GRF>> GRFPtrList_;

        //- Density field name, default = rho
        word rhoName_;

        word UName_;

        //- fvm mask field
        volScalarField fvmMask_;

        mutable surfaceScalarField phiOmega_;

        mutable surfaceScalarField relativePhi_;

        const Switch consistentDiscretization_;

        //- relax GRF patch cells in block solver
        scalar relax_;

        //- First time read
        bool isInitialised_;

    // Protected Member Functions
        void addLocalFvSchemes();

        //- set mask
        void setFvmMask();

        //- set rot flux
        void setFields() const;

        void relaxGRFBoundary(fvBlockMatrix<vector>&) const;

        //- Add rotational convection for scalar
        template<class RhoFieldType>
        void addRotation
        (
            const RhoFieldType& rhof,
            fvMatrix<scalar>& eqn
        ) const;

        //- Add rotational convection for vectors
        template<class RhoFieldType>
        void addRotation
        (
            const RhoFieldType& rhof,
            fvMatrix<vector>& eqn
        ) const;


        //- Add volumetric centrifugal force
        template<class RhoFieldType>
        void addRotation
        (
            const RhoFieldType& rho,
            fvBlockMatrix<vector>& eqn
        ) const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GRFSource(const GRFSource&);

        //- Disallow default bitwise assignment
        void operator=(const GRFSource&);


public:

    //- Runtime type information
    TypeName("GRFSource");


    // Constructors

        //- Construct from components
        GRFSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~GRFSource()
    {}


    // Member Functions

        //- Initialise after all fvOptions created
        virtual bool initialise() override;

        //- Return rotation induced flux field for the entire domain
        virtual tmp<surfaceScalarField> phiOmega() const;

        // Add explicit and implicit contributions

            //- Scalar
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;

            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;

            //- Vector
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            ) override;

            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            ) override;

            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;

            //- Add the frame acceleration to the matrix
            virtual void addAcceleration
            (
                const volVectorField& U,
                volVectorField& ddtU
            ) const override
            {}

            //- Add the frame acceleration contribution to the momentum
            virtual void addAcceleration
            (
                fvVectorMatrix& UEqn,
                bool rhs = false
            ) const override
            {}

            //- Add the frame acceleration contribution to the momentum
            virtual void addAcceleration
            (
                const volScalarField& rho,
                fvVectorMatrix& UEqn,
                bool rhs = false
            ) const override
            {}

            //- Add the frame acceleration contribution to the block matrix
            virtual void addAcceleration
            (
                fvBlockMatrix<vector>& eqn
            ) const override;

            //- Add the frame acceleration contribution to the block matrix
            virtual void addAcceleration
            (
                const volScalarField& rho,
                fvBlockMatrix<vector>& eqn
            ) const override;

            //- Make the cell centre velocity relative
            virtual void makeRelative(volVectorField& U) const override;

            //- Make the given absolute flux relative
            virtual void makeRelative(surfaceScalarField& phi) const override
            {}

            //- Make the given absolute boundary flux relative
            virtual void makeRelative
            (
                FieldField<fvsPatchField, scalar>& phi
            ) const override
            {}

            //- Make the given absolute patch flux relative
            virtual void makeRelative
            (
                Field<scalar>& phi,
                const label patchi
            ) const override
            {}

            //- Make the given absolute mass-flux relative
            virtual void makeRelative
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const override
            {}

            //- Make the given relative flux absolute
            virtual void makeAbsolute(surfaceScalarField& phi) const override
            {}

            //- Make the given relative mass-flux absolute
            virtual void makeAbsolute
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const override
            {}

            //- Make the cell centre velocity absolute
            virtual void makeAbsolute(volVectorField& U) const override
            {}


        // Frame manipulations

            //- Filter-out the MRF region contribution from the given field
            // setting the corresponding values to zero
            virtual void zero
            (
                surfaceScalarField& phi
            ) const override;


        // I-O

            //- Read dictionary
            virtual bool read(const dictionary& dict) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GRFSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

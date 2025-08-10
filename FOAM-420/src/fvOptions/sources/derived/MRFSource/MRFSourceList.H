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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

Class
    Foam::MRFSourceList

Description
    List container for MRF sources

SourceFiles
    MRFSourceList.C

\*---------------------------------------------------------------------------*/

#ifndef MRFSourceList_H
#define MRFSourceList_H

#include "fvMesh/fvMesh.H"
#include "db/dictionary/dictionary.H"
#include "fvMatrices/fvMatricesFwd.H"
#include "MRFSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// Forward declaration of friend functions and operators
class MRFSourceList;
Ostream& operator<<(Ostream& os, const MRFSourceList& models);

/*---------------------------------------------------------------------------*\
                         Class MRFSourceList Declaration
\*---------------------------------------------------------------------------*/

class MRFSourceList
:
    public PtrList<MRFSource>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        MRFSourceList(const MRFSourceList&);

        //- Disallow default bitwise assignment
        void operator=(const MRFSourceList&);


protected:

    // Protected data

        //- Reference to the database
        const objectRegistry& obr_;


        void initialise();
public:

    //- Constructor
    MRFSourceList
    (
        const objectRegistry& obr,
        const dictionary& dict
    );

    //- Destructor
    ~MRFSourceList();


    // Member Functions

        //- Return active status
        bool active(const bool warn = false) const;

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Return the MRF with a given name
        const MRFSource& getFromName(const word& name) const;

        //- Add the frame acceleration
        void addAcceleration
        (
            const volVectorField& U,
            volVectorField& ddtU
        ) const;

        //- Add the frame acceleration contribution to the momentum equation
        void addAcceleration
        (
            fvVectorMatrix& UEqn,
            bool rhs = false
        ) const;

        //- Add the frame acceleration contribution to the momentum equation
        void addAcceleration
        (
            const volScalarField& rho,
            fvVectorMatrix& UEqn,
            bool rhs = false
        ) const;

        //- Return the frame acceleration
        tmp<volVectorField> DDt
        (
            const volVectorField& U
        ) const;

        //- Return the frame acceleration
        tmp<volVectorField> DDt
        (
            const volScalarField& rho,
            const volVectorField& U
        ) const;

        //- Make the given absolute velocity relative within the MRF region
        void makeRelative(volVectorField& U) const;

        //- Make the given absolute flux relative within the MRF region
        void makeRelative(surfaceScalarField& phi) const;

        //- Return the given absolute flux relative within the MRF region
        tmp<surfaceScalarField> relative
        (
            const tmp<surfaceScalarField>& phi
        ) const;

        //- Return the given absolute boundary flux relative within
        //  the MRF region
        tmp<FieldField<fvsPatchField, scalar>> relative
        (
            const tmp<FieldField<fvsPatchField, scalar>>& tphi
        ) const;

        //- Return the given absolute patch flux relative within
        //  the MRF region
        tmp<Field<scalar>> relative
        (
            const tmp<Field<scalar>>& tphi,
            const label patchi
        ) const;

        //- Make the given absolute mass-flux relative within the MRF region
        void makeRelative
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        //- Make the given relative velocity absolute within the MRF region
        void makeAbsolute(volVectorField& U) const;

        //- Make the given relative flux absolute within the MRF region
        void makeAbsolute(surfaceScalarField& phi) const;

        //- Return the given relative flux absolute within the MRF region
        tmp<surfaceScalarField> absolute
        (
            const tmp<surfaceScalarField>& phi
        ) const;

        //- Make the given relative mass-flux absolute within the MRF region
        void makeAbsolute
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        //- Filter-out the MRF region contribution from the given field
        // setting the corresponding values to zero
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> zeroFilter
        (
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tphi
        ) const;

        // I-O

            //- Read dictionary
            bool read(const dictionary& dict);

            //- Write data to Ostream
            bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const MRFSourceList& models
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MRFSourceListTemplates.C"
#endif

#endif

// ************************************************************************* //

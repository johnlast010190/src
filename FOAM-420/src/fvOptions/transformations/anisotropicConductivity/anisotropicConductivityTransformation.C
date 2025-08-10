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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*----------------------------------------------------------------------------*/

#include "anisotropicConductivityTransformation.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "basicThermo/basicThermo.H"
#include "solidThermo/solidThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "coordinate/systems/cartesianCS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(anisotropicConductivityTransformation, 0);

        addToRunTimeSelectionTable
        (
            option,
            anisotropicConductivityTransformation,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::anisotropicConductivityTransformation
::anisotropicConductivityTransformation
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    csysPtr_(),
    coorFramePtr_(nullptr),
    localCartesian_()
{
    if (coeffs_.found("referenceFrame"))
    {
        word frameName = coeffs_.lookup<word>("referenceFrame");
        coorFramePtr_ = &coordinateFrame::New(mesh_, frameName);
    }
    else if (coeffs_.found("coordinateSystem"))
    {
        csysPtr_ = coordinateSystem::New(mesh_, coeffs_, coordinateSystem::typeName_());
    }
    else
    {
        csysPtr_.reset(new coordSystem::cartesian());
    }

    fieldNames_.setSize(1, "Anialpha");
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::coordinateSystem&
Foam::fv::anisotropicConductivityTransformation::csys() const
{
    if (coorFramePtr_)
    {
        return coorFramePtr_->coorSys();
    }
    return *csysPtr_;
}

Foam::tmp<Foam::volTensorField> Foam::fv::anisotropicConductivityTransformation
::updateTransformationTensor() const
{
    tmp<volTensorField> tft
    (
        new volTensorField
        (
            IOobject
            (
                "lcc",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor
            (
                "lcc",
                dimless,
                tensor::I
            ),
            calculatedFvPatchField<tensor>::typeName
        )
    );

    forAll(tft->boundaryField(), pi)
    {
        tft->boundaryFieldRef()[pi] = tensor::I;
    }


    // set internal fields values
    forAll(cells_, cellI)
    {
        label cI = cells_[cellI];
        tft->operator[](cellI) = (csys().R(mesh_.C()[cI]));
    }

    // set boundary values

    // first make a list of all non-coupled boundary element
    List<DynamicList<label>> cellSetBoundaryFaces(mesh_.boundaryMesh().size());

    // loop through all set cells and add any boundary faces to List
    forAll(cells_, i)
    {
        const labelList& cfaces(mesh_.cells()[cells_[i]]);

        forAll(cfaces, cfI)
        {
            label gFaceI = cfaces[cfI];
            if (gFaceI >= mesh_.nInternalFaces())
            {
                label patchI = mesh_.boundaryMesh().whichPatch(gFaceI);

                if (!mesh_.boundaryMesh()[patchI].coupled())
                {
                    cellSetBoundaryFaces[patchI].append
                    (
                        gFaceI - mesh_.boundaryMesh()[patchI].start()
                    );
                }
            }
        }
    }

    // now calculate and assign local coordinate tensors
    volTensorField::Boundary& tftbf = tft->boundaryFieldRef();
    forAll(cellSetBoundaryFaces, si)
    {
        cellSetBoundaryFaces[si].shrink();

        if (!tftbf[si].coupled() && (tftbf[si].size()))
        {
            forAll(cellSetBoundaryFaces[si], sfi)
            {
                label faceI = cellSetBoundaryFaces[si][sfi];
                tftbf[si][faceI]
                    = (csys().R(mesh_.boundary()[si].Cf()[faceI]));
            }
        }
    }

    return tft;
}

Foam::symmTensor Foam::fv::anisotropicConductivityTransformation
::transformPrincipal
(
    const tensor& tt,
    const vector& st
) const
{
    return symmTensor
    (
        tt.xx()*st.x()*tt.xx()
      + tt.xy()*st.y()*tt.xy()
      + tt.xz()*st.z()*tt.xz(),

        tt.xx()*st.x()*tt.yx()
      + tt.xy()*st.y()*tt.yy()
      + tt.xz()*st.z()*tt.yz(),

        tt.xx()*st.x()*tt.zx()
      + tt.xy()*st.y()*tt.zy()
      + tt.xz()*st.z()*tt.zz(),

        tt.yx()*st.x()*tt.yx()
      + tt.yy()*st.y()*tt.yy()
      + tt.yz()*st.z()*tt.yz(),

        tt.yx()*st.x()*tt.zx()
      + tt.yy()*st.y()*tt.zy()
      + tt.yz()*st.z()*tt.zz(),

        tt.zx()*st.x()*tt.zx()
      + tt.zy()*st.y()*tt.zy()
      + tt.zz()*st.z()*tt.zz()
    );
}

void Foam::fv::anisotropicConductivityTransformation::transformField
(
    symmTensorField& fld,
    const tensorField& tt,
    const vectorField& st
) const
{
    forAll(fld, i)
    {
        fld[i] = transformPrincipal(tt[i], st[i]);
    }
}

void Foam::fv::anisotropicConductivityTransformation::transformField
(
    volSymmTensorField& fld,
    const volTensorField& tt,
    const volVectorField& st
) const
{
    transformField
    (
        fld.ref(),
        tt.internalField(),
        st.internalField()
    );

    forAll(fld.boundaryField(), pI)
    {
        transformField
        (
            fld.boundaryFieldRef()[pI],
            tt.boundaryField()[pI],
            st.boundaryField()[pI]
        );
    }
}


bool Foam::fv::anisotropicConductivityTransformation::alwaysApply() const
{
    return true;
}

void Foam::fv::anisotropicConductivityTransformation::correct
(
    volSymmTensorField& kbcp
)
{
    // modify transformation tensor if necessary
    if (mesh_.changing() || !localCartesian_.valid())
    {
        localCartesian_.reset(updateTransformationTensor().ptr());
    }

    // first lookup thermo
    const solidThermo& thermo =
        obr_.lookupObject<solidThermo>(basicThermo::dictName);

    //fail if isotropic
    if (thermo.isotropic())
    {
        FatalErrorInFunction
            << "Anisotropic conductivity transformation called for an "
            << "isotropic medium."
            << exit(FatalError);
    }

    tmp<volVectorField> KappaByCp(thermo.Kappa()/thermo.Cp());

    // transform diffusivity
    transformField
    (
        kbcp,
        localCartesian_(),
        KappaByCp()
    );

    if (debug)
    {
        Info<< "anisotropicConductivityTransformation : "
            << "dumping converted Kxx, Kyy, Kzz"
            << endl;
        {
            volVectorField Kxx
            (
                IOobject
                (
                    "Kxx",
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh_,
                dimless
            );
            Kxx.primitiveFieldRef() = transform
            (
                localCartesian_->internalField(),
                vectorField
                (
                    localCartesian_->internalField().size(),
                    point(1, 0, 0)
                )
            );
            forAll(Kxx.boundaryField(), patchI)
            {
                Kxx.boundaryFieldRef()[patchI] = transform
                (
                    localCartesian_->boundaryField()[patchI],
                    vectorField
                    (
                        localCartesian_->boundaryField()[patchI].size(),
                        point(1, 0, 0)
                    )
                );
            }
            Kxx.write();
        }
        {
            volVectorField Kyy
            (
                IOobject
                (
                    "Kyy",
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh_,
                dimless
            );
            Kyy.primitiveFieldRef() = transform
            (
                localCartesian_->internalField(),
                vectorField
                (
                    localCartesian_->internalField().size(),
                    point(0, 1, 0)
                )
            );
            forAll(Kyy.boundaryField(), patchI)
            {
                Kyy.boundaryFieldRef()[patchI] = transform
                (
                    localCartesian_->boundaryField()[patchI],
                    vectorField
                    (
                        localCartesian_->boundaryField()[patchI].size(),
                        point(0, 1, 0)
                    )
                );
            }
            Kyy.write();
        }
        {
            volVectorField Kzz
            (
                IOobject
                (
                    "Kzz",
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh_,
                dimless
            );
            Kzz.primitiveFieldRef() = transform
            (
                localCartesian_->internalField(),
                vectorField
                (
                    localCartesian_->internalField().size(),
                    point(0, 0, 1)
                )
            );
            forAll(Kzz.boundaryField(), patchI)
            {
                Kzz.boundaryFieldRef()[patchI] = transform
                (
                    localCartesian_->boundaryField()[patchI],
                    vectorField
                    (
                        localCartesian_->boundaryField()[patchI].size(),
                        point(0, 0, 1)
                    )
                );
            }
            Kzz.write();
        }
    }


}

void Foam::fv::anisotropicConductivityTransformation::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::anisotropicConductivityTransformation::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        if (coeffs_.found("referenceFrame"))
        {
            word frameName = coeffs_.lookup<word>("referenceFrame");
            coorFramePtr_ = &coordinateFrame::New(mesh_, frameName);
        }
        else
        {
            csysPtr_ = coordinateSystem::New(mesh_, coeffs_, coordinateSystem::typeName_());
        }

        localCartesian_.clear();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

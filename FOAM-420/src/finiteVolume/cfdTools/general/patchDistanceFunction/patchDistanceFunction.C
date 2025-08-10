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
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "cfdTools/general/patchDistanceFunction/patchDistanceFunction.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

defineTypeNameAndDebug(patchDistanceFunction, 0);

template<>
const char* NamedEnum<patchDistanceFunction::distTypes, 7>::names[] =
{
    "x",
    "y",
    "z",
    "wallDistance",
    "point",
    "vector",
    "radius"
};


const NamedEnum <patchDistanceFunction::distTypes, 7>
patchDistanceFunction::distTypeNames_(0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistanceFunction::patchDistanceFunction(const fvPatch& p)
:
    patch_(p),
    distanceType_(dtZ),
    dictCoeffs_(dictionary()),
    origin_(Zero),
    axis_(vector(0, 0, 1)),
    value_(),
    dir_(),
    delta_(),
    coorFramePtr_(nullptr)
{}


Foam::patchDistanceFunction::patchDistanceFunction
(
    const fvPatch& p,
    const dictionary& dict
)
:
    patch_(p),
    distanceType_(distTypeNames_.readOrDefault(dict.lookup("distanceType"), dtX)),
    dictCoeffs_
    (
        dict.optionalSubDict(word(distTypeNames_[distanceType_]) + "Coeffs")
    ),
    origin_(dictCoeffs_.lookupOrDefault<vector>("origin", Zero)),
    axis_
    (
        normalised
        (
            dictCoeffs_.lookupOrDefault<vector>("axis", vector(0, 0, 1))
        )
    ),
    value_(),
    dir_(),
    delta_(),
    coorFramePtr_(loadFrame(dict))
{
    if
    (
        coorFramePtr_
     && (dictCoeffs_.found("origin") || dictCoeffs_.found("axis"))
    )
    {
        WarningInFunction
            << "origin and axis in  \"" << patchDistanceFunction::type() << "\""
            << " on patch name: " << "\"" << patch_.name() << "\""
            << " isn't used. Frame name: \"" << coorFramePtr_->name()
            << "\" origin is used instead." << endl;
    }
}


Foam::patchDistanceFunction::patchDistanceFunction
(
    const patchDistanceFunction& pdf
)
:
    tmp<patchDistanceFunction>::refCount(),
    patch_(pdf.patch_),
    distanceType_(pdf.distanceType_),
    dictCoeffs_(pdf.dictCoeffs_),
    origin_(pdf.origin_),
    axis_(pdf.axis_),
    value_(),
    dir_(),
    delta_(),
    coorFramePtr_(pdf.coorFramePtr_)
{
    if
    (
        coorFramePtr_
     && (dictCoeffs_.found("origin") || dictCoeffs_.found("axis"))
    )
    {
        WarningInFunction
            << "origin and axis in  \"" << patchDistanceFunction::type() << "\""
            << " on patch name: " << "\"" << patch_.name() << "\""
            << " isn't used. Frame name: \"" << coorFramePtr_->name()
            << "\" origin is used instead." << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::coordinateFrame* Foam::patchDistanceFunction::loadFrame
(
    const dictionary& dict
)
{
    if
    (
        dictCoeffs_.found("referenceFrame")
     && dictCoeffs_.lookupOrAddDefault<Switch>("definedInFrame", false)
    )
    {
        return
            coordinateFrame::lookupNew
            (
                patch_.boundaryMesh().mesh(),
                dictCoeffs_
            );
    }

    return nullptr;
}


Foam::vector Foam::patchDistanceFunction::localToGlobalDir(vector dir) const
{
    if (coorFramePtr_)
    {
        dir = coorFramePtr_->coorSys().globalVector(dir);
    }
    return dir;
}


const Foam::scalarField& Foam::patchDistanceFunction::value() const
{
    if (!value_.valid())
    {
        value_.reset(new scalarField(patch_.size(), Zero));
        vectorField Cf(patch_.Cf());
        if (coorFramePtr_)
        {
            origin_ = coorFramePtr_->coorSys().origin();
            axis_ = coorFramePtr_->axis();
            Cf = coorFramePtr_->coorSys().localPosition(Cf);
        }
        else
        {
            Cf -= origin_;
        }
        switch(distanceType_)
        {
            case dtX:
            {
                value_() = Cf.component(0)();
                break;
            }
            case dtY:
            {
                value_() = Cf.component(1)();
                break;
            }
            case dtZ:
            {
                value_() = Cf.component(2)();
                break;
            }
            case dtWall:
            {
                wallDist d(patch_.boundaryMesh().mesh());
                value_() =
                    d.y().boundaryField()[patch_.index()].patchInternalField();
                break;
            }
            case dtPoint:
            {
                value_() = mag(Cf);
                break;
            }
            case dtVector:
            {
                value_() = (axis_ & Cf);
                break;
            }
            case dtRadius:
            {
                value_() = mag(Cf - axis_*(axis_ & Cf));
                break;
            }
            default:
            {
                FatalError
                    << distTypeNames_[distanceType_]
                    << " is not a valid distanceType." << endl
                    << "Valid formats are: x, y, z, wallDistance, "
                    << "point, radius and vector."
                    << exit(FatalError);
            }
        }
    }

    return value_();
}


const Foam::vectorField& Foam::patchDistanceFunction::dir() const
{
    if (!dir_.valid())
    {
        if (coorFramePtr_)
        {
            origin_ = coorFramePtr_->coorSys().origin();
            axis_ = coorFramePtr_->axis();
        }
        switch(distanceType_)
        {
            case dtX:
            {
                const vector xDirection(localToGlobalDir(vector(1, 0, 0)));
                dir_.reset(new vectorField(patch_.size(), xDirection));
                break;
            }
            case dtY:
            {
                const vector yDirection(localToGlobalDir(vector(0, 1, 0)));
                dir_.reset(new vectorField(patch_.size(), yDirection));
                break;
            }
            case dtZ:
            {
                const vector zDirection(localToGlobalDir(vector(0, 0, 1)));
                dir_.reset(new vectorField(patch_.size(), zDirection));
                break;
            }
            case dtWall:
            {
                wallDist d(patch_.boundaryMesh().mesh());
                dir_.reset
                (
                    d.n().boundaryField()
                    [
                        patch_.index()
                    ].patchInternalField().ptr()
                );
                break;
            }
            case dtPoint:
            {
                dir_.reset(new vectorField(patch_.Cf() - origin_));
                dir_() /= mag(dir_());
                break;
            }
            case dtVector:
            {
                dir_.reset(new vectorField(patch_.size(), axis_));
                break;
            }
            case dtRadius:
            {
                dir_.reset(new vectorField(patch_.Cf() - origin_));
                dir_() -= axis_*(axis_ & dir_());
                dir_() /= max(mag(dir_()), SMALL);
                break;
            }
            default:
            {
                FatalError
                    << distTypeNames_[distanceType_]
                    << " is not a valid distanceType." << endl
                    << "Valid formats are: x, y, z, wallDistance, "
                    << "point, radius and vector."
                    << exit(FatalError);
            }
        }
    }

    return dir_();
}


const Foam::scalarField& Foam::patchDistanceFunction::delta() const
{
    // Since direction is considering frame rotation the delta should be
    // computed consistently with the frame rotation
    if (!delta_.valid())
    {
        const fvMesh& mesh(patch_.boundaryMesh().mesh());
        const faceList& meshFaces(mesh.faces());
        const pointField& meshPoints(mesh.points());
        delta_.reset(new scalarField(patch_.size(), Zero));

        const vectorField& direction = dir();

        forAll(patch_, fi)
        {
            label gFaceI = patch_.patch().start() + fi;
            pointField facePoints = meshFaces[gFaceI].points(meshPoints);

            scalarList pointDotDir(facePoints & direction[fi]);

            delta_->operator[](fi) = 0.5*(max(pointDotDir) - min(pointDotDir));
        }
    }

    return delta_();
}


bool Foam::patchDistanceFunction::update()
{
    const fvMesh& mesh(patch_.boundaryMesh().mesh());
    bool updated = false;

    // The frames could have access if there were changed or not
    if
    (
        mesh.moving()
     || mesh.changing()
     || (coorFramePtr_ && coorFramePtr_->anyDynamic())
    )
    {
        value_.clear();
        dir_.clear();
        delta_.clear();
    }

    if (!value_.valid())
    {
        updated = true;
    }

    return updated;
}


void Foam::patchDistanceFunction::write(Ostream& os, bool writeFrame) const
{
    const word distType(distTypeNames_[distanceType_]);
    os.writeEntry("distanceType", distType);

    // In sub-dict case there is no problem to write the frame since
    // there will not be any duplicates
    if (dictCoeffs_.dictName() == word(distType + "Coeffs"))
    {
        os.indent();
        os << dictCoeffs_.dictName();
        os << dictCoeffs_;
    }
    else
    {
        if (!coorFramePtr_)
        {
            if (origin_ != vector::zero)
            {
                os.writeEntry("origin", origin_);
            }
            if (axis_ != vector(0, 0, 1))
            {
                os.writeEntry("axis", axis_);
            }
        }
        else if (writeFrame)
        {
            os.writeEntry("referenceFrame", coorFramePtr_->name());
            os.writeEntry("definedInFrame", Switch(true));
        }
    }
}


// ************************************************************************* //

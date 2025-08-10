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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Class
    volumetricNurbs

Description
    Creates a Nurbs based morphing method

\*---------------------------------------------------------------------------*/

#include "vNurbsDeformation/volumetricNurbs/volumetricNurbs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

point volumetricNurbs::getRationalPoint
(
    const scalar u,
    const scalar v,
    const scalar w
) const
{
    // calculate denom W

    scalar W(0);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

        scalar weight = weights[index];

                W += N_u*N_v*N_w*weight;
            }
        }
    }

    // Calculate point

    point result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

                result += (N_u*N_v*N_w*weights[index]*box.getPoint(icp, jcp, kcp));
            }
        }
    }

    result /= W;

    return result;
}

vector volumetricNurbs::getRationalUDer
(
    const scalar u,
    const scalar v,
    const scalar w
) const
{
    // calculate denom W, W'

    scalar W(0), Wprime(0);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);
    const scalar N_uP = uBasis.derivative(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

        scalar weight = weights[index];

                W += N_u*N_v*N_w*weight;
                Wprime += N_uP*N_v*N_w*weight;
            }
        }
    }

    // Calculate Uder

    vector result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);
    const scalar N_uP = uBasis.derivative(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

                const scalar& weight = weights[index];

                vector value = (N_uP*N_v*N_w*weight*box.getPoint(icp, jcp, kcp) * W - N_u*N_v*N_w*weight*box.getPoint(icp,jcp,kcp)*Wprime);

                result += value;
            }
        }
    }

    result /= (W*W);

    return result;
}

vector volumetricNurbs::getRationalVDer
(
    const scalar u,
    const scalar v,
    const scalar w
) const
{
    // calculate denom W, W'

    scalar W(0), Wprime(0);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);
        const scalar N_vP = vBasis.derivative(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

        scalar weight = weights[index];

                W += N_u*N_v*N_w*weight;
                Wprime += N_u*N_vP*N_w*weight;
            }
        }
    }

    // Calculate Vder

    vector result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);
        const scalar N_vP = vBasis.derivative(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

                const scalar& weight = weights[index];

                vector value = (N_u*N_vP*N_w*weight*box.getPoint(icp, jcp, kcp) * W - N_u*N_v*N_w*weight*box.getPoint(icp,jcp,kcp)*Wprime);

                result += value;
            }
        }
    }

    result /= (W*W);

    return result;
}

vector volumetricNurbs::getRationalWDer
(
    const scalar u,
    const scalar v,
    const scalar w
) const
{
    // calculate denom W, W'

    scalar W(0), Wprime(0);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);
        const scalar N_wP = wBasis.derivative(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

        scalar weight = weights[index];

                W += N_u*N_v*N_w*weight;
                Wprime += N_u*N_v*N_wP*weight;
            }
        }
    }

    // Calculate Wder

    vector result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);
        const scalar N_wP = wBasis.derivative(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

                const scalar& weight = weights[index];

                vector value = (N_u*N_v*N_wP*weight*box.getPoint(icp, jcp, kcp) * W - N_u*N_v*N_w*weight*box.getPoint(icp,jcp,kcp)*Wprime);

                result += value;
            }
        }
    }

    result /= (W*W);

    return result;
}

vector volumetricNurbs::getRationalDisplacement
(
    const scalar u,
    const scalar v,
    const scalar w
) const
{
    // calculate denom W

    scalar W(0);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

        scalar weight = weights[index];

                W += N_u*N_v*N_w*weight;
            }
        }
    }

    // Calculate point

    vector result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

                result += N_u*N_v*N_w*weights[index]*box.getDisplacement(icp, jcp, kcp);
            }
        }
    }

    result /= W;

    return result;
}

scalarField volumetricNurbs::getRationalBasesList
(
    const scalar u,
    const scalar v,
    const scalar w
) const
{

    // calculate denom W

    scalar W(0);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

        scalar w = weights[index];

                W += N_u*N_v*N_w*w;
            }
        }
    }

    // Calculate basesList

    scalarField result(weights.size(), 0.0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

                result[index] = N_u*N_v*N_w*weights[index]/W;
            }
        }
    }

    return result;
}

void volumetricNurbs::checkWeights()
{
    if (!checkRational) return;

    checkRational = false;
    isRational = false;

    scalar valueToCompare = weights[0];

    label size = weights.size();

    for (label i = 1; i < size; ++i)
    {
        scalar valueToCheck = weights[i];

        // in-equality means that it is rational

        scalar diff = mag(valueToCheck-valueToCompare);

        if (diff > 1e-12)
        {
            isRational = true;

            return;
        }
    }
}

// ------------------------------------------------------------------------- //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

volumetricNurbs::volumetricNurbs
(
    const dictionary& dict
) : box
    (
        dict
    ),
    uBasis
    (
        readLabel(dict.lookup("nUCPs")),
        readLabel(dict.lookup("UDegree"))
    ),
    vBasis
    (
        readLabel(dict.lookup("nVCPs")),
        readLabel(dict.lookup("VDegree"))
    ),
    wBasis
    (
        readLabel(dict.lookup("nWCPs")),
        readLabel(dict.lookup("WDegree"))
    ),
    weights
    (
        dict.lookupOrDefault("Weights", scalarField(uBasis.nbPoints()*vBasis.nbPoints()*wBasis.nbPoints(), 1.0))
    )
{
    checkRational = true;
}

//- Construct as copy
volumetricNurbs::volumetricNurbs
(
     const volumetricNurbs& vNurbs
) : box(vNurbs.box),
    uBasis(vNurbs.uBasis),
    vBasis(vNurbs.vBasis),
    wBasis(vNurbs.wBasis),
    weights(vNurbs.weights)
{
    checkRational = true;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const scalarField& volumetricNurbs::getWeights() const
{
    return weights;
}

scalarField volumetricNurbs::getBasesList
(
    const scalar u,
    const scalar v,
    const scalar w
)
{
    this->checkWeights();

    if (isRational)
        return this->getRationalBasesList(u, v, w);

    scalarField result(box.nbXCPs()*box.nbYCPs()*box.nbZCPs(), 0.0);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

            label index = icp + jcp*box.nbXCPs() + kcp*box.nbXCPs()*box.nbYCPs();

        result[index] = N_u*N_v*N_w;
        }
        }
    }

    return result;
}

point volumetricNurbs::getPoint
(
        const scalar u,
           const scalar v,
        const scalar w
)
{
    this->checkWeights();

    if (isRational)
        return this->getRationalPoint(u, v, w);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    point result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

        result += N_u*N_v*N_w*box.getPoint(icp, jcp, kcp);
        }
    }
    }

    return result;
}

vector volumetricNurbs::getDisplacement
(
         const scalar u,
        const scalar v,
        const scalar w
)
{
    this->checkWeights();

    if (isRational)
        return this->getRationalDisplacement(u, v, w);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    vector result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

        result += N_u*N_v*N_w*box.getDisplacement(icp, jcp, kcp);
        }
    }
    }

    return result;
}

vector volumetricNurbs::getUDer
(
         const scalar u,
        const scalar v,
        const scalar w
)
{
    this->checkWeights();

    if (isRational)
        return this->getRationalUDer(u, v, w);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    vector result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.derivative(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

        result += N_u*N_v*N_w*box.getPoint(icp, jcp, kcp);
        }
    }
    }

    return result;
}


vector volumetricNurbs::getVDer
(
         const scalar u,
        const scalar v,
        const scalar w
)
{
    this->checkWeights();

    if (isRational)
        return this->getRationalVDer(u, v, w);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    vector result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.derivative(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.value(kcp, wBasis.degree(), w);

        result += N_u*N_v*N_w*box.getPoint(icp, jcp, kcp);
        }
    }
    }

    return result;
}


vector volumetricNurbs::getWDer
(
         const scalar u,
        const scalar v,
        const scalar w
)
{
    this->checkWeights();

    if (isRational)
        return this->getRationalWDer(u, v, w);

    const label uSpan = uBasis.span(u);
    const label vSpan = vBasis.span(v);
    const label wSpan = wBasis.span(w);

    vector result(0, 0, 0);

    for (label i = 0; i <= uBasis.degree(); i++)
    {
        label icp = uSpan - uBasis.degree() + i;
    const scalar N_u = uBasis.value(icp, uBasis.degree(), u);

    for (label j = 0; j <= vBasis.degree(); j++)
    {
        label jcp = vSpan - vBasis.degree() + j;
        const scalar N_v = vBasis.value(jcp, vBasis.degree(), v);

        for (label k = 0; k <= wBasis.degree(); k++)
        {
            label kcp = wSpan - wBasis.degree() + k;
        const scalar N_w = wBasis.derivative(kcp, wBasis.degree(), w);

        result += N_u*N_v*N_w*box.getPoint(icp, jcp, kcp);
        }
    }
    }

    return result;
}


inline bool volumetricNurbs::contains(const point& pt) const
{
    point minpt = box.getMinPt();
    point maxpt = box.getMaxPt();

    return
    (
        pt.x() >= minpt.x() && pt.x() <= maxpt.x()
     && pt.y() >= minpt.y() && pt.y() <= maxpt.y()
     && pt.z() >= minpt.z() && pt.z() <= maxpt.z()
    );
}


point volumetricNurbs::invert
(
    const point& p
)
{
    scalar u(-1), v(-1), w(-1);

    const polyMesh* cBoxMesh = box.getPolyMesh();

    if (cBoxMesh == nullptr)
        FatalErrorInFunction
            << "The control box must be first converted to a polyMesh object"
            << exit(FatalError);


    const cellList& cells = cBoxMesh->cells();
    const faceList& faces = cBoxMesh->faces();
    const pointField& points = cBoxMesh->points();

    if (!contains(p)) return point(-1, -1, -1); //points

    forAll(cells, cI)
    {
        // check for bndbox

        const cell& c = cells[cI];

        pointField cellPoints = c.points(faces, points);

        if (!contains(p)) continue; //cellPoints

        if
        (
            cBoxMesh->pointInCell
            (
                p,
                cI,
                Foam::polyMesh::FACE_PLANES
            )
        )
        {
            label ic, jc, kc;

            kc = std::floor(scalar(cI) / scalar((box.nbXCPs()-1)*(box.nbYCPs()-1)));
            jc = std::floor(scalar(cI - kc*(box.nbXCPs()-1)*(box.nbYCPs()-1)) / scalar(box.nbXCPs()-1));
            ic = cI - (jc*(box.nbXCPs()-1) + kc*(box.nbXCPs()-1)*(box.nbYCPs()-1));

            u = scalar(ic+1) / scalar(box.nbXCPs()+1);
            v = scalar(jc+1) / scalar(box.nbYCPs()+1);
            w = scalar(kc+1) / scalar(box.nbZCPs()+1);

            break;
        }
    }

    if (u < 0 || u > 1)
    {
        return point(-1, -1, -1);
    }

    if (v < 0 || v > 1)
    {
        return point(-1, -1, -1);
    }

    if (w < 0 || w > 1)
    {
        return point(-1, -1, -1);
    }


    point p_star = this->getPoint(u, v, w);

    scalar dist = mag(p-p_star);

    scalar eps = 1e-12;

    Foam::vector vRepeat(-1, -1, -1);

    label calcCount(0);

    while (dist > eps)
    {
        // if statement to catch repeating exceptions
        // may be unnecessary
        if (calcCount > 2000)
        {
            if (vRepeat.x() < 0)
            {
                vRepeat = Foam::vector(u, v, w);

                goto next;
            }

            Foam::vector vToCheck(u, v, w);

            if (mag(vToCheck-vRepeat) < 1e-12)
                return point(-1, -1, -1);
        }

        next:

        // dv : distance vector
        // Bu, Bv, Bw : directional derivatives

        vector dv = p-p_star;

        vector Bu = this->getUDer(u, v, w);
        vector Bv = this->getVDer(u, v, w);
        vector Bw = this->getWDer(u, v, w);

        tensor M
        (
            Bu.x(), Bv.x(), Bw.x(),
            Bu.y(), Bv.y(), Bw.y(),
            Bu.z(), Bv.z(), Bw.z()
        );

        if (det(M) == 0)
        {
            FatalErrorInFunction
                << "Directional derivative matrix is singular. "
                << "Check for degenerate regions due to duplicate CPs."
                << exit(FatalError);
        }

        // Invert derivative matrix
        // calculate update step for N - R
        tensor Minv = inv(M);

        vector delta = Minv&dv;

        // calculate relaxation
        // start with OR and reduce
        scalar step(1.2), unew, vnew, wnew;

        label n(0);

        while (n < 100)
        {
            unew = u+step*delta.x();
            vnew = v+step*delta.y();
            wnew = w+step*delta.z();

            if (unew < 0) unew = 0.0;
            if (unew > 1) unew = 1.0;
            if (vnew < 0) vnew = 0.0;
            if (vnew > 1) vnew = 1.0;
            if (wnew < 0) wnew = 0.0;
            if (wnew > 1) wnew = 1.0;

            if (dist > mag(p-this->getPoint(unew,vnew,wnew)))
            {
                break;
            }

            step = step / 2.0;
            n++;
        }

        if (n > 99) return point(-1, -1, -1);

        u = unew;
        v = vnew;
        w = wnew;

        // bound u v w to 0-1 for robustness
        if (u < 0) u = 0.0;
        if (u > 1) u = 1.0;
        if (v < 0) v = 0.0;
        if (v > 1) v = 1.0;
        if (w < 0) w = 0.0;
        if (w > 1) w = 1.0;

        p_star = this->getPoint(u, v, w);

        dist = mag(p-p_star);

        if ((dist > eps) && (mag(delta)==0))
        {
            WarningInFunction
                << "Newton not converging. 1st derivative probably zero."
                << " Change initiallization of u v w." << endl;
        }

        calcCount++;
    }

    point result(u, v, w);

    return result;
}

controlBox& volumetricNurbs::getBox()
{
    return box;
}

const controlBox& volumetricNurbs::getBox() const
{
    return box;
}

void volumetricNurbs::assignUknots(const scalarField& uKnots)
{
    uBasis.changeKnots(uKnots);
}

void volumetricNurbs::assignVknots(const scalarField& vKnots)
{
    vBasis.changeKnots(vKnots);
}

void volumetricNurbs::assignWknots(const scalarField& wKnots)
{
    wBasis.changeKnots(wKnots);
}

void volumetricNurbs::setWeight
(
    const label& i,
    const label& j,
    const label& k,
    const scalar& w
)
{
    label index = i + j*box.nbXCPs() + k*box.nbXCPs()*box.nbYCPs();
    weights[index] = w;

    checkRational = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

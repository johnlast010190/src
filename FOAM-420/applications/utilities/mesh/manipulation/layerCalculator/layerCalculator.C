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
    (c) 2011-2015 OpenFOAM Foundation

Application
    layerCalculator

Group
    grpMeshManipulationUtilities

Description
    A layer mesh profile can be characterised by specifying any three of
    the following five quantities: first cell height, final cell height,
    total layer height, stretching and number of layers. This calculator
    inputs any three of these five variables and calculates the remaining two

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "numLayers",
        "label",
        "Number of layers to grow"
    );

    argList::addOption
    (
        "expansionRatio",
        "scalar",
        "Expansion ratio of layers"
    );

    argList::addOption
    (
        "firstLayerThickness",
        "scalar",
        "First layer thickness (m)"
    );

    argList::addOption
    (
        "finalLayerThickness",
        "scalar",
        "Final layer thickness (m)"
    );

    argList::addOption
    (
        "totalLayerThickness",
        "scalar",
        "Total layer thickness (m)"
    );
    #include "include/setRootCase.H"

    label numLayers = args.optionLookupOrDefault("numLayers",-1);
    scalar expansionRatio =
        args.optionLookupOrDefault("expansionRatio",-1.0);
    scalar finalLayerThickness =
        args.optionLookupOrDefault("finalLayerThickness",-1.0);
    scalar totalLayerThickness =
        args.optionLookupOrDefault("totalLayerThickness",-1.0);
    scalar firstlayerThickness =
        args.optionLookupOrDefault("firstLayerThickness",-1.0);

    //Calaculate layer profiles based on input
    if (numLayers > 0)
    {
        if (expansionRatio > 0)
        {
            if (finalLayerThickness > 0)
            {
                Info<<"Calculating based on number of layers,expansion ratio"
                    <<" and final layer thickness"<<endl;

                firstlayerThickness = finalLayerThickness
                    / Foam::pow(expansionRatio, numLayers-1);
                if (expansionRatio == 1)
                {
                    totalLayerThickness = firstlayerThickness*numLayers;
                }
                else
                {
                    totalLayerThickness = firstlayerThickness*
                    (
                        (scalar(1.)-Foam::pow(expansionRatio,numLayers))
                        / (scalar(1.)-expansionRatio)
                     );
                }
            }
            else if (firstlayerThickness > 0)
            {
                Info<<"Calculating based on number of layers,expansion ratio"
                    <<" and first layer thickness"<<endl;

                finalLayerThickness = firstlayerThickness
                    * Foam::pow(expansionRatio, numLayers -1);
                if (expansionRatio == 1)
                {
                    totalLayerThickness = firstlayerThickness*numLayers;
                }
                else
                {
                    totalLayerThickness = firstlayerThickness*
                    (
                        (scalar(1)-Foam::pow(expansionRatio,numLayers))
                        / (scalar(1)-expansionRatio)
                     );
                }
            }
            else if (totalLayerThickness > 0)
            {
                Info<<"Calculating based on number of layers,expansion ratio"
                    <<" and total layer thickness"<<endl;

                if (expansionRatio == 1)
                {
                   finalLayerThickness  = totalLayerThickness / numLayers;
                   firstlayerThickness = finalLayerThickness;
                }
                else
                {
                    firstlayerThickness = totalLayerThickness
                        * (1 - expansionRatio)
                        / (1 - Foam::pow(expansionRatio, numLayers));
                    finalLayerThickness = firstlayerThickness
                        *Foam::pow(expansionRatio, numLayers-1);
                }
            }
        }
        else if (finalLayerThickness > 0)
        {
            if (firstlayerThickness > 0)
            {
                Info<<"Calculating based on number of layers, "
                    <<" first layer thickness and final layer thickness" <<endl;

                if (numLayers == 1)
                {
                    expansionRatio = 1.;
                }
                else
                {
                    scalar ratio = finalLayerThickness/firstlayerThickness;
                    expansionRatio = Foam::pow(ratio, scalar(1./(numLayers-1)));
                }

                if (expansionRatio == 1)
                {
                    totalLayerThickness = firstlayerThickness*numLayers;
                }
                else
                {
                    totalLayerThickness = firstlayerThickness*
                    (
                        (scalar(1)-Foam::pow(expansionRatio,numLayers))
                        / (scalar(1)-expansionRatio)
                    );
                }
            }
            else if (totalLayerThickness > 0)
            {
                Info<<"Calculating based on number of layers, "
                    <<" final layer thickness and total layer thickness" <<endl;

                scalar xn = 0.0;
                scalar x = -1.0;

                scalar a1 = finalLayerThickness;
                scalar a2 = totalLayerThickness;
                label repeat = 0;

                scalar tol = 0.001;

                while (!(repeat++ >= 100 || mag(x - xn) < tol))
                {
                    x = xn;
                    scalar f1 = a1 * Foam::pow(x, numLayers) + a2*(1 - x) - a1;
                    scalar f2 = a1 * numLayers
                        * Foam::pow(x, numLayers -1) - a2;
                    xn = x - f1/(f2+SMALL);
                }

                if (xn < 1.0 + tol && xn > 1.0 - tol)
                {
                    xn = 10.0;
                    x = -1.0;
                    repeat = 0;
                    while (!(repeat++ >= 100 || mag(x - xn) < tol))
                    {
                        x = xn;
                        scalar f1 = a1 * Foam::pow(x, numLayers)
                            + a2*(1 - x) - a1;
                        scalar f2 = a1 * numLayers
                            * Foam::pow(x, numLayers -1) - a2;
                        xn = x - f1/(f2+SMALL);
                    }
                }

                //now calculate inverse
                expansionRatio = 1. / xn;
                firstlayerThickness = finalLayerThickness
                    / Foam::pow(expansionRatio,numLayers-1);
            }
        }
        else if (firstlayerThickness > 0)
        {
            if (totalLayerThickness > 0)
            {
                Info<<"Calculating based on number of layers, "
                    <<" first layer thickness and total layer thickness" <<endl;

                scalar xn = 10.0;
                scalar x = -1.0;

                scalar a1 = firstlayerThickness;
                scalar a2 = totalLayerThickness;
                label repeat = 0;

                scalar tol = 0.001;

                while (!(repeat++ >= 100 || mag(x - xn) < tol))
                {
                    x = xn;
                    scalar f1 = a1 * Foam::pow(x, numLayers) + a2*(1 - x) - a1;
                    scalar f2 = a1 * numLayers
                        * Foam::pow(x, numLayers -1) - a2;
                    xn = x - f1/(f2+SMALL);
                }
                if (xn < 1.0 + tol && xn > 1.0 - tol)
                {
                    xn = 0.0;
                    x = -1.0;
                    repeat = 0;
                    while (!(repeat++ >= 100 || mag(x - xn) < tol))
                    {
                        x = xn;
                        scalar f1 = a1 * Foam::pow(x, numLayers)
                            + a2*(1 - x) - a1;
                        scalar f2 = a1 * numLayers
                            * Foam::pow(x, numLayers -1) - a2;
                        xn = x - f1/(f2+SMALL);
                    }
                }
                expansionRatio = xn;
                finalLayerThickness = a1 * Foam::pow(xn, numLayers - 1);
            }
        }
    }
    else if (expansionRatio > 0)
    {
        if (finalLayerThickness > 0)
        {
            if (firstlayerThickness > 0)
            {
                Info<<"Calculating based on expansion ratio, "
                    <<" first layer thickness and final layer thickness" <<endl;

                if (expansionRatio == 1)
                {
                    if (finalLayerThickness != firstlayerThickness)
                    {
                        Info<<"Expansion ratio set to 1 but final layer thickness"
                            <<" not equal to first layer thickness "<<endl;
                    }
                    else
                    {
                        Info<<"Unique solution is not possible"<<endl;
                    }
                }
                else
                {
                    scalar ratio = finalLayerThickness / firstlayerThickness;
                    numLayers =
                        (Foam::log(ratio) / Foam::log(expansionRatio)) + 1;
                    if (expansionRatio == 1)
                    {
                        totalLayerThickness = firstlayerThickness*numLayers;
                    }
                    else
                    {
                        totalLayerThickness = firstlayerThickness*
                        (
                            (scalar(1)-Foam::pow(expansionRatio,numLayers))
                            / (scalar(1)-expansionRatio)
                        );
                    }
                }
            }
            else if (totalLayerThickness > 0)
            {
                Info<<"Calculating based on expansion ratio, "
                    <<" final layer thickness and total layer thickness" <<endl;

                if (expansionRatio == 1)
                {
                    if (finalLayerThickness > totalLayerThickness)
                    {
                        Info<<"Final layer thickness is greater than the "
                            <<"total layer thickness"<<endl;
                    }
                    else
                    {
                        firstlayerThickness = finalLayerThickness;
                        numLayers = totalLayerThickness / firstlayerThickness;
                    }
                }
                else
                {
                    scalar invStr = 1. / expansionRatio;
                    scalar flt = finalLayerThickness;
                    scalar mlt = totalLayerThickness;
                    scalar ratio = 1.  - (mlt * (1. - invStr) / flt);

                    if (ratio < SMALL)
                    {
                        Info<< "Cannot find solution for  expansionRatio, "
                             << "finalLayerThickness and totalLayerThickness"
                             <<endl;
                    }
                    else
                    {
                        numLayers = Foam::log(ratio)
                            / Foam::log(invStr);
                        firstlayerThickness = finalLayerThickness
                            / Foam::pow(expansionRatio, numLayers-1);
                    }
                }
            }
        }
        else if (firstlayerThickness > 0)
        {
            if (totalLayerThickness > 0)
            {
                Info<<"Calculating based on expansion ratio, "
                    <<" first layer thickness and total layer thickness" <<endl;

                if (expansionRatio == 1)
                {
                    if (firstlayerThickness > totalLayerThickness)
                    {
                        Info<<"First layer thickness is greater than the "
                            <<"total layer thickness"<<endl;
                    }
                    else
                    {
                        numLayers = totalLayerThickness/firstlayerThickness;
                        finalLayerThickness = firstlayerThickness;
                    }
                }
                else
                {
                    scalar ratio = 1.
                        - (totalLayerThickness * (1. - expansionRatio)
                           / firstlayerThickness);

                    if (ratio < SMALL)
                    {
                        Info<< "Cannot find solution for  expansionRatio, "
                             << "firstlayerThickness and totalLayerThickness"
                             <<endl;
                    }
                    else
                    {
                        numLayers = Foam::log(ratio)/Foam::log(expansionRatio);
                        finalLayerThickness = firstlayerThickness
                            * Foam::pow(expansionRatio, numLayers -1);
                    }
                }
            }
        }
    }
    else if (finalLayerThickness > 0)
    {
        if (firstlayerThickness > 0)
        {
            if (totalLayerThickness > 0)
            {
                Info<<"Calculating based on first layer thickness, "
                    <<" final layer thickness and total layer thickness" <<endl;

                scalar flt = finalLayerThickness;
                scalar mlt = totalLayerThickness;
                expansionRatio = (firstlayerThickness - mlt)/(flt - mlt);

                if (expansionRatio < SMALL)
                {
                    Info<< "Cannot find solution for finalLayerThickness, "
                         << "firstlayerThickness and totalLayerThickness"
                         <<endl;
                }
                else
                {
                    scalar ratio = flt*expansionRatio/firstlayerThickness;
                    numLayers = Foam::log(ratio)/Foam::log(expansionRatio);
                }
            }
        }
    }

    //Write out layer profiles
    Info<<"Number of Layers : "<< numLayers <<endl;
    Info<<"ExpansionRatio : "<< expansionRatio <<endl;
    Info<<"First Layer Thickness (m) : "<< firstlayerThickness <<endl;
    Info<<"Final Layer Thickness (m ) : "<< finalLayerThickness <<endl;
    Info<<"Total Layer Thickness (m) : "<< totalLayerThickness <<endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

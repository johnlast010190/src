word fluxScheme("Kurganov");
if (mesh.schemes().dict().readIfPresent("fluxScheme", fluxScheme))
{
    if ((fluxScheme == "Tadmor") || (fluxScheme == "Kurganov"))
    {
        Info<< "fluxScheme: " << fluxScheme << endl;
    }
    else
    {
        FatalErrorInFunction
            << "fluxScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: Tadmor, Kurganov"
            << abort(FatalError);
    }
}

scalar eMin = mesh.solution().dict().lookupOrDefault<scalar>("eMin", -GREAT);
scalar rhoMin = mesh.solution().dict().lookupOrDefault<scalar>("rhoMin", -GREAT);

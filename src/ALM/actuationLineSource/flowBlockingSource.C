/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "flowBlockingSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(flowBlockingSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        flowBlockingSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flowBlockingSource::flowBlockingSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),

    //- Set the pointer to runTime
    runTime_(mesh_.time()),

    // Set the current simulation time.
    time(runTime_.timeName()),
    //- f field of the volume of cell
    f
    (
        IOobject
        (
            "f",
            time,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("f",dimless,0.0)
    ),
    // - bodyForce vector 
    bodyForce(
       IOobject
       (
           "bodyForce",
           time,
           mesh_,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedVector("bodyForce",dimForce/dimVolume/dimDensity,vector::zero)
    ),
    // - sp field
       sp(
       IOobject
       (
           "sp",
           time,
           mesh_,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       mesh_,
       dimensionedScalar("sp",dimless/dimTime,0.0)
    )
{
    // assign -1 value to cell centres so min() max() 
    // function works in the case of multiple geometries
    forAll(mesh_.C(),cellI)
    {
        f[cellI] =-1.0;
    }

    // - Calculate maximum & minumum edge
    const edgeList& edges       = mesh.edges();    
    const pointField& pp        = mesh.points();
    scalarField eLengths(edges.size(),0.0);
    forAll(edges,edgei)
    {
        eLengths[edgei] = edges[edgei].mag(pp);
    }    
    maxEdge_ = Foam::max(eLengths);
    minEdge_ = Foam::min(eLengths);
    characteristicLength_ = minEdge_;
    Info<<"Characteristic Length = "<<characteristicLength_<<endl;
    Info<<"Max edge = "<<maxEdge_<<endl;
    Info<<"Min edge = "<<minEdge_<<endl;

    // - Read general variables
    read(dict);
    // - Read specific variables
    forAll(geometryNames_,i)
    {
        word geometryType = word(coeffs_.subDict(geometryNames_[i]).lookup("type"));

        if (geometryType == "circularCylinder")
        {
            vector axisStart = vector(coeffs_.subDict(geometryNames_[i]).lookup("axisStart"));
            vector axisEnd   = vector(coeffs_.subDict(geometryNames_[i]).lookup("axisEnd"));

            scalar radius = scalar(readScalar(coeffs_.subDict(geometryNames_[i]).lookup("radius")));
            List<scalar> temp(int(2),0.);
            List<List<scalar> > surfaceNonDimProfile(2,temp);
            surfaceNonDimProfile[0][0] = 0.;
            surfaceNonDimProfile[0][1] = radius;
            surfaceNonDimProfile[1][0] = 1.;
            surfaceNonDimProfile[1][1] = radius;
            
            spAxisymmetric(axisStart,axisEnd,surfaceNonDimProfile);
            
        }        
        else if (geometryType == "sphere")
        {
            vector centre = vector(coeffs_.subDict(geometryNames_[i]).lookup("centre"));
            scalar radius = scalar(readScalar(coeffs_.subDict(geometryNames_[i]).lookup("radius")));

            spSphere(centre,radius);
        }        
        else if (geometryType == "cube")
        {
            vector centre = vector(coeffs_.subDict(geometryNames_[i]).lookup("centre"));
            vector sides  = vector(coeffs_.subDict(geometryNames_[i]).lookup("sides"));

            spCube(centre,sides[0],sides[1],sides[2]);
        }        
        else if (geometryType == "nacelle")
        {
            vector axisStart   = vector(coeffs_.subDict(geometryNames_[i]).lookup("axisStart"));
            vector axisEnd     = vector(coeffs_.subDict(geometryNames_[i]).lookup("axisEnd"));
            scalar radius      = scalar(readScalar(coeffs_.subDict(geometryNames_[i]).lookup("radius")));
            scalar halfLength  = scalar(readScalar(coeffs_.subDict(geometryNames_[i]).lookup("halfLength")));
            scalar axisLength  = mag(axisEnd-axisStart);
            vector axisDirUnit = (axisEnd-axisStart)/axisLength;
            vector centre = axisStart + halfLength*axisDirUnit;

            //spSphere(centre,radius);
            spEllipse(centre,axisEnd,radius);

            List<scalar> temp(int(2),0.);
            List<List<scalar> > surfaceNonDimProfile(2,temp);

            surfaceNonDimProfile[0][0] = 0.;
            surfaceNonDimProfile[0][1] = radius;
            surfaceNonDimProfile[1][0] = 1.;
            surfaceNonDimProfile[1][1] = radius;
            
            spAxisymmetric(axisStart,centre,surfaceNonDimProfile);
        }        
        else if (geometryType == "arbitraryAxisymmetricSurface")
        {
            vector axisStart = vector(coeffs_.subDict(geometryNames_[i]).lookup("axisStart"));
            vector axisEnd   = vector(coeffs_.subDict(geometryNames_[i]).lookup("axisEnd"));
            List<List<scalar> > surfaceNonDimProfile;
            coeffs_.subDict(geometryNames_[i]).lookup("nonDimDistribution")>>surfaceNonDimProfile;

            spAxisymmetric(axisStart,axisEnd,surfaceNonDimProfile);
        }        
        else
        {
            FatalErrorIn("Foam::fv::flowBlockingSource::read()")
            <<"Invalid geometryType type "<<geometryType
            <<"Valid geometryType types are :\n\tcircularCylinder\n\tsquareCylinder\n\tsphere\n\tcube\n\tnacelle\n\tarbitraryAxisymmetricSurface\n"
            <<exit(FatalIOError);
        }
    }
    f.write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::flowBlockingSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
        const volVectorField& U = eqn.psi();

        bodyForce *= 0.;
        sp *= 0.;

        forAll(mesh_.cells(),cellI)
        {
            sp[cellI] += f[cellI]*characteristicVelocity_/characteristicLength_/mesh_.V()[cellI];
            bodyForce[cellI] += U[cellI]*sp[cellI]*mesh_.V()[cellI];
        }
        eqn -= fvm::Sp(sp,U);
        
}

void Foam::fv::flowBlockingSource::spSphere
(
    vector centre,
    scalar radius
)
{
    const faceList   & ff = mesh_.faces();
    const pointField & pp = mesh_.points();
    scalar outsideVol  = 0.;
    scalar insideVol   = 0.;
    scalar boundaryVol = 0.;
    scalar totalVol    = 0.;

    
    forAll(mesh_.C(),cellI)
    {
        totalVol += mesh_.V()[cellI];

        vector point = mesh_.C()[cellI];
        if (mag(point-centre) < 0.5*radius)
        {
            f[cellI] = 1.;
            insideVol+=mesh_.V()[cellI];
        }
        else if (mag(point-centre)> 1.5*radius)
        {
            f[cellI] = Foam::max(0.,f[cellI]);
            outsideVol+=mesh_.V()[cellI];
        }
        else 
        {
            int in  = 0;
            int out = 0;
            const cell & cc = mesh_.cells()[cellI];
            labelList pLabels(cc.labels(ff));
            pointField pLocal(pLabels.size(),vector::zero);
            forAll(pLabels,pLabelsi)
            {
                vector vertex = pp[pLabels[pLabelsi]];    
                
                if ( (mag(vertex-centre)<radius) )   in++;
                else                                  out++;
            }//end forAll vertices

            if (out == 0 ) 
            {
                f[cellI] = 1.;
                insideVol+=mesh_.V()[cellI];
            }
            else if (in == 0 ) 
            {
                f[cellI] = Foam::max(0.,f[cellI]);
                outsideVol+=mesh_.V()[cellI];
            }
            else 
            {
                f[cellI] = Foam::max(double(in)/double(in+out),f[cellI]);
                boundaryVol+=mesh_.V()[cellI];
            }
        }//end if (boundary)
    }// end forAll mesh cells

    reduce(totalVol,sumOp<scalar>());
    reduce(insideVol,sumOp<scalar>());
    reduce(outsideVol,sumOp<scalar>());
    reduce(boundaryVol,sumOp<scalar>());

    Info<<" In spSphere found points:\n"
    <<"\tInside: "    <<(insideVol)/(totalVol)*100.<<"% \n"
    <<"\tOutside: "    <<(outsideVol)/(totalVol)*100.<<"% \n"
    <<"\tBoundary: "<<(boundaryVol)/(totalVol)*100.<<"% \n";


}

void Foam::fv::flowBlockingSource::spEllipse
(
    vector axisStart,
    vector axisEnd,
    scalar radius
)
{
    scalar nPoints = 100;
    List<scalar> temp(int(2),0.);
    List<List<scalar> > surfaceNonDimProfile(nPoints,temp);
    for (int i=0;i<nPoints;i++)
    {
        scalar x = double(i)/double(nPoints-1);
        surfaceNonDimProfile[i][0] = x;
        surfaceNonDimProfile[i][1] = Foam::sqrt(1-x*x)*radius;
    }

    spAxisymmetric(axisStart,axisEnd,surfaceNonDimProfile);
}
void Foam::fv::flowBlockingSource::spAxisymmetric
(
    vector axisStart,
    vector axisEnd,
    List<List<scalar> >  surfaceNonDimProfile
)
{
    // - Assign axis info, volume info, mesh info
    vector axisDir     = axisEnd-axisStart;
    scalar axisLength  = mag(axisDir);
    vector axisDirUnit = axisDir/axisLength;

    scalar outsideVol  = 0;
    scalar insideVol   = 0;
    scalar boundaryVol = 0;
    scalar totalVol    = 0;

    const faceList   & ff = mesh_.faces();
    const pointField & pp = mesh_.points();

    forAll(mesh_.C(),cellI)
    {
        totalVol += mesh_.V()[cellI];

        vector projection    = axisStart+((mesh_.C()[cellI]-axisStart)&axisDirUnit)*axisDirUnit;
        scalar distFromStart = mag(projection-axisStart);
        scalar distFromEnd   = mag(projection-axisEnd);
        if ( (distFromStart>axisLength) || (distFromEnd>axisLength) ) // case where cell is above or below the axis
        {
            f[cellI] = Foam::max(0.,f[cellI]);    
            outsideVol+=mesh_.V()[cellI];
        }
        else
        {
            scalar s = mag(projection-axisStart)/axisLength;
            scalar r = radiusAtS(s,surfaceNonDimProfile);
            if ( mag (mesh_.C()[cellI] -projection)> 1.5*r )
            {
                f[cellI] = Foam::max(0.,f[cellI]);
                outsideVol+=mesh_.V()[cellI];
            }
            else if (mag(mesh_.C()[cellI]-projection)< 0.5*r )
            {
                f[cellI] = 1.;
                insideVol+=mesh_.V()[cellI];
            //    Info<<mesh_.C()[cellI]<<" "<<f[cellI]<<endl;
            }
            else
            {
                int in  = 0;
                int out = 0;
                const cell & cc = mesh_.cells()[cellI];
                labelList pLabels(cc.labels(ff));
                pointField pLocal(pLabels.size(),vector::zero);
                forAll(pLabels,pLabelsi)
                {
                    vector vertex = pp[pLabels[pLabelsi]];    
                    vector projectionOfVertex = axisStart+((vertex-axisStart)&axisDirUnit)*axisDirUnit;

                    scalar vertexDistFromStart = mag(projectionOfVertex-axisStart);
                    scalar vertexDistFromEnd   = mag(projectionOfVertex-axisEnd);
                    
        
                    if ( (vertexDistFromStart>axisLength) || (vertexDistFromEnd>axisLength) )  
                    {
                        out++;
                        continue;
                    }
                    else if ( mag(vertex-projectionOfVertex) < r )                            
                    {
                        in++;
                        continue;
                    }
                    else                         
                    {
                        out++;
                        continue;
                    }

                }
                if (out == 0 ) 
                {
                    f[cellI] = 1.;
                    insideVol+=mesh_.V()[cellI];
                //    Info<<mesh_.C()[cellI]<<" "<<f[cellI]<<endl;
                }
                else if (in == 0 ) 
                {    
                    f[cellI] = Foam::max(0.,f[cellI]);
                    outsideVol+=mesh_.V()[cellI];
                }
                else 
                {
                    f[cellI] = Foam::max(double(in)/double(in+out),f[cellI]);
                    boundaryVol+=mesh_.V()[cellI];
                //    Info<<mesh_.C()[cellI]<<" "<<f[cellI]<<endl;    
                }
            }
        }
    }

    
    reduce(totalVol,sumOp<scalar>());
    reduce(insideVol,sumOp<scalar>());
    reduce(outsideVol,sumOp<scalar>());
    reduce(boundaryVol,sumOp<scalar>());

    Info<<" In spAxisymmetric found points:\n"
    <<"\tInside: "    <<(insideVol)/(totalVol)*100.<<"% \n"
    <<"\tOutside: "    <<(outsideVol)/(totalVol)*100.<<"% \n"
    <<"\tBoundary: "<<(boundaryVol)/(totalVol)*100.<<"% \n";


}

void Foam::fv::flowBlockingSource::spCube
(
    vector centre,
    scalar sideX,
    scalar sideY,
    scalar sideZ
)
{
// side becomes half side so that I can use the centre
    sideX = sideX*0.5+minEdge_/1000.;
    sideY = sideY*0.5+minEdge_/1000.;
    sideZ = sideZ*0.5+minEdge_/1000.;

    scalar xup05   = centre[0] + 0.5*sideX;
    scalar yup05   = centre[1] + 0.5*sideY;
    scalar zup05   = centre[2] + 0.5*sideZ;

    scalar xdown05 = centre[0] - 0.5*sideX;
    scalar ydown05 = centre[1] - 0.5*sideY;
    scalar zdown05 = centre[2] - 0.5*sideZ;

    scalar xup15   = centre[0] + 1.5*sideX;
    scalar yup15   = centre[1] + 1.5*sideY;
    scalar zup15   = centre[2] + 1.5*sideZ;

    scalar xdown15 = centre[0] - 1.5*sideX;
    scalar ydown15 = centre[1] - 1.5*sideY;
    scalar zdown15 = centre[2] - 1.5*sideZ;

    scalar xup     = centre[0] + sideX;
    scalar yup     = centre[1] + sideY;
    scalar zup     = centre[2] + sideZ;

    scalar xdown   = centre[0] - sideX;
    scalar ydown   = centre[1] - sideY;
    scalar zdown   = centre[2] - sideZ;

    scalar outsideVol  = 0.;
    scalar insideVol   = 0.;
    scalar boundaryVol = 0.;
    scalar totalVol    = 0.;

    const faceList   & ff = mesh_.faces();
    const pointField & pp = mesh_.points();
    forAll(mesh_.C(),cellI)
    {
        totalVol += mesh_.V()[cellI];
        
        vector point  = mesh_.C()[cellI];
        if ( (point[0] < xup05 ) && (point[0] > xdown05 ) &&
             (point[1] < yup05 ) && (point[1] > ydown05 ) &&
             (point[2] < zup05 ) && (point[2] > zdown05 ) )
        {
            f[cellI] = 1.;
            insideVol+=mesh_.V()[cellI];
            //Info<<mesh_.C()[cellI]<<" f="<<f[cellI]<<endl;
        }
        else if ( ( (point[0] > xup15 ) || (point[0] < xdown15 ) ) ||
                  ( (point[1] > yup15 ) || (point[1] < ydown15 ) ) ||
                  ( (point[2] > zup15 ) || (point[2] < zdown15 ) ) )
        {
            f[cellI] = Foam::max(0.,f[cellI]);
            outsideVol+=mesh_.V()[cellI];
        }
        else
        {

            int in  = 0;
            int out = 0;
            const cell & cc = mesh_.cells()[cellI];
            labelList pLabels(cc.labels(ff));
            pointField pLocal(pLabels.size(),vector::zero);
            forAll(pLabels,pLabelsi)
            {
                vector vertex = pp[pLabels[pLabelsi]];    
                
                if ( (vertex[0] < xup ) && (vertex[0] > xdown ) &&
                     (vertex[1] < yup ) && (vertex[1] > ydown ) &&
                     (vertex[2] < zup ) && (vertex[2] > zdown ) )
                {    
                    in++;
                }
                else
                {
                    out++;
                }
            } //for All vertices

            if ( out == 0 ) // all in 
            {
                f[cellI] = 1.;
                insideVol+=mesh_.V()[cellI];
            //    Info<<mesh_.C()[cellI]<<" f="<<f[cellI]<<endl;
            }
            else if ( in == 0 ) // all out
            {
                f[cellI] = Foam::max(0.,f[cellI]);
                outsideVol+=mesh_.V()[cellI];
            }
            else 
            {
                f[cellI] = Foam::max(double(in)/double(in+out),f[cellI]);
                boundaryVol+=mesh_.V()[cellI];
            //    Info<<mesh_.C()[cellI]<<" f="<<f[cellI]<<endl;
            }    
        }
    } //forAll cells

    reduce(totalVol,sumOp<scalar>());
    reduce(insideVol,sumOp<scalar>());
    reduce(outsideVol,sumOp<scalar>());
    reduce(boundaryVol,sumOp<scalar>());

    Info<<" In spCube found points:\n"
    <<"\tInside: "    <<(insideVol)/(totalVol)*100.<<"% \n"
    <<"\tOutside: "    <<(outsideVol)/(totalVol)*100.<<"% \n"
    <<"\tBoundary: "<<(boundaryVol)/(totalVol)*100.<<"% \n";


}


scalar Foam::fv::flowBlockingSource::radiusAtS
(
    scalar s,
    List<List<scalar> > surfaceNonDimProfile
    
)
{

    DynamicList<scalar> xOld;
    DynamicList<scalar> yOld;
    forAll(surfaceNonDimProfile,i)
    {
    xOld.append(surfaceNonDimProfile[i][0]);
    yOld.append(surfaceNonDimProfile[i][1]);
    }
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(s - xOld[i]);
        if(diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (s < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
        }
        return
            yOld[indexM]
          + (
                (yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM])
            )
           *(s - xOld[indexM]);
    }
    else if (s > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
        }
        return
            yOld[indexM]
          + (
                (yOld[indexP] - yOld[indexM])/(xOld[indexP] - xOld[indexM])
            )
           *(s - xOld[indexM]);
    }
    else if (s == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}


void Foam::fv::flowBlockingSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::flowBlockingSource::read(const dictionary& dict)
{
   if (option::read(dict)) 
   {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);
        coeffs_.lookup("geometryNames")>>geometryNames_;
        coeffs_.lookup("characteristicVelocity")>>characteristicVelocity_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

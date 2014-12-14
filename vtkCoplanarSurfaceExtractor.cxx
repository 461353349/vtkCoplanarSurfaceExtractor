
//////boolean intersection operation on coplanar regions of two polydata meshes
//////idea and implementation by Roman Grothausmann

/////vtkDistancePolyDataFilter and vtkPartialVolumeModeller taken as examples


////the filter's own header
#include "vtkCoplanarSurfaceExtractor.h"

#include <vtkSmartPointer.h>

////vtk-filter headers
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>

////headers needed for this specific filter
#include <vtkMath.h>
#include <vtkTriangleFilter.h> 
#include <vtkPolyDataNormals.h> 
#include <vtkFloatArray.h>
#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkLine.h>
#include <vtkCleanPolyData.h>
#include <vtkMergePoints.h>
#include <vtkTransform.h>
#include <vtkDelaunay2D.h>
#include <vtkAppendPolyData.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkConvexHull2D.h>
#include <vtkGeometryFilter.h>


vtkStandardNewMacro(vtkCoplanarSurfaceExtractor);




vtkCoplanarSurfaceExtractor::vtkCoplanarSurfaceExtractor(){

    this->SetNumberOfInputPorts(2);

    this->DistanceTolerance= 0.001;
    this->FaceOrientationTolerance= 0.001;
    this->LineOrientationTolerance= 0.001;

    this->MeshMode= VTK_USE_DELAUNAY2D;
    }



int vtkCoplanarSurfaceExtractor::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector){

    // get the info objects
    vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
    vtkInformation *inInfo1 = inputVector[1]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkPolyData *input0 = vtkPolyData::SafeDownCast(inInfo0->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *input1 = vtkPolyData::SafeDownCast(inInfo1->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


    vtkSmartPointer< vtkTriangleFilter> triangulate0= vtkSmartPointer<vtkTriangleFilter>::New(); //converts vtkPolyLine to vtkLine
    triangulate0->SetInputData(input0);
    triangulate0->PassLinesOff();//Very important because cells of type vtkLine in the output confuses vtkPolyDataNormals with ComputeCellNormalsOn!!! If it is off, then the input lines will be ignored and the output will have no lines. //causes this error in paraview or dumpXML: vtkXMLPolyDataReader (0x2370350): Cannot read cell data array "Normals" from PointData in piece 0.  The data array in the element may be too short.

    vtkSmartPointer<vtkTriangleFilter> triangulate1= vtkSmartPointer<vtkTriangleFilter>::New(); //converts vtkPolyLine to vtkLine
    triangulate1->SetInputData(input1);
    triangulate1->PassLinesOff();

    vtkSmartPointer<vtkPolyDataNormals> PDnormals0= vtkSmartPointer<vtkPolyDataNormals>::New();
    PDnormals0->SetInputConnection(triangulate0->GetOutputPort());
    PDnormals0->ComputePointNormalsOff(); 
    PDnormals0->ComputeCellNormalsOn();
    PDnormals0->Update();

    vtkSmartPointer<vtkPolyDataNormals> PDnormals1= vtkSmartPointer<vtkPolyDataNormals>::New();
    PDnormals1->SetInputConnection(triangulate1->GetOutputPort());
    PDnormals1->ComputePointNormalsOff(); 
    PDnormals1->ComputeCellNormalsOn();
    PDnormals1->Update();

    vtkPolyData *mesh0= PDnormals0->GetOutput(); 
    vtkPolyData *mesh1= PDnormals1->GetOutput();

    vtkIdType N0= mesh0->GetNumberOfCells();
    vtkDebugMacro(<< "Mesh 1 contains " << N0 << " cells and " << mesh0->GetNumberOfPoints() << " points."); 
    vtkIdType N1= mesh1->GetNumberOfCells();
    vtkDebugMacro(<< "Mesh 2 contains " << N1 << " cells and " << mesh1->GetNumberOfPoints() << " points."); 

    double x[3];
    double y[3];

    vtkDataArray *normals0= mesh0->GetCellData()->GetNormals();
    vtkDataArray *normals1= mesh1->GetCellData()->GetNormals();

    if(normals0){
        vtkDebugMacro(<< "There are " << normals0->GetNumberOfTuples() << " normals in normals0.");
        }
    else {
        vtkErrorMacro(<< "No normals in normals0! Aborting.");
        return VTK_ERROR;
        }
        
    if(normals1){
        vtkDebugMacro(<< "There are " << normals1->GetNumberOfTuples() << " normals in normals1.");
        }
    else {
        vtkErrorMacro(<< "No normals in normals1! Aborting.");
        return VTK_ERROR;
        }
        

    vtkSmartPointer<vtkAppendPolyData> appendPD= vtkSmartPointer<vtkAppendPolyData>::New();
    ////inital input for vtkAppendPolyData needed in the case of disjunct sets!!!
    vtkSmartPointer<vtkPolyData> empty_mesh= vtkSmartPointer<vtkPolyData>::New();
    appendPD->AddInputData(empty_mesh);

    vtkIdType counter= 0;
    int abortExecute=0;
    vtkIdType progressInterval = N0*N1/1000 + 1;

    for (vtkIdType i= 0; i < N0 && !abortExecute; i++){

	vtkDebugMacro(<< "i: " << i);
        
        vtkCell* cell0= mesh0->GetCell(i);
        if (cell0->GetCellType() != VTK_TRIANGLE){
	    vtkDebugMacro(<< "Cell " << i << " is no triangle (" << vtkCellTypes::GetClassNameFromTypeId(cell0->GetCellType()) << ").");
            continue;
            }

        for (vtkIdType j= 0; j < N1; j++){

            vtkDebugMacro(<< "j: " << j);
        
            vtkCell* cell1= mesh1->GetCell(j);
            if (cell1->GetCellType() != VTK_TRIANGLE){
                vtkDebugMacro(<< "Cell " << j << " is no triangle (" << vtkCellTypes::GetClassNameFromTypeId(cell1->GetCellType()) << ").");
                continue;
                }


	    //////cells are coplanar, if the containing planes are identical, i.e. are parallel and have the same displacement from origin
	    ////check if cells are nearly parallel using face-tollerance (this->FaceOrientationTolerance)
	    ////check if containing planes have similar displacement from origin
	    ////if faces are nearly coplanar:
            if (coplanar_check(normals0->GetTuple(i), normals1->GetTuple(j), cell0->GetPoints()->GetPoint(0), cell1->GetPoints()->GetPoint(0), this->FaceOrientationTolerance, this->DistanceTolerance)){
                ////find intersections of edges -> add to out point list
                vtkIdType m0= cell0->GetNumberOfEdges();
                vtkIdType m1= cell1->GetNumberOfEdges();

                vtkCell *edge0; 
                vtkCell *edge1; 

                double l00[3];
                double l01[3];
                double l10[3];
                double l11[3];
                double p0[3];
                double p1[3];
                double t1, t2;

		vtkSmartPointer<vtkPoints> outPoints= vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> outCells= vtkSmartPointer<vtkCellArray>::New();

		vtkSmartPointer<vtkPolyData> point_pd= vtkSmartPointer<vtkPolyData>::New();
		point_pd->SetPoints(outPoints);


                ////merge dublicate points in point-set before Delaunay2D
                ///vtkCleanPolyData removes points not used by any cell, use vtkMergePoints instead

		vtkSmartPointer<vtkMergePoints> mergePoints =  vtkSmartPointer<vtkMergePoints>::New();
		mergePoints->SetDataSet(point_pd);
		const int divisions= 0;//the smaller, the faster
		mergePoints->SetDivisions(divisions, divisions, divisions);//docs do not explain what this is for, big speedup!!!
		mergePoints->InitPointInsertion(point_pd->GetPoints(), cell0->GetBounds());


                vtkDebugMacro(<< "Checking intersections... ");
                for (vtkIdType k= 0; k < m0; k++){
                    edge0= cell0->GetEdge(k);
                    for (vtkIdType l= 0; l < m1; l++){
                        edge1= cell1->GetEdge(l);
                        
                        edge0->GetPoints()->GetPoint(0,l00);
                        edge0->GetPoints()->GetPoint(1,l01);
                        edge1->GetPoints()->GetPoint(0,l10);
                        edge1->GetPoints()->GetPoint(1,l11);

			////check if edges are not (nearly) parallel, because then the closest point is not well defined and adding it makes no sense; using line-tollerance (this->LineOrientationTolerance)
                        double b=angle_in_deg_unori(l00, l01, l10, l11);
                        if (b <= this->LineOrientationTolerance){
                            //if (b == 0){
                            vtkDebugMacro(<< "Angle: " << b << " No intersection considered!");
                            continue;
                            }

                        ////get closest points between lines
                        double d= vtkLine::DistanceBetweenLineSegments(l00, l01, l10, l11, p0, p1, t1, t2);//return value is the shortest distance squared between the two line-segments.
                        vtkDebugMacro(<< "Distance: " << d);

			int added;
                        if(point_in_both_cells(p0, cell0, cell1, this->DistanceTolerance)){
                            vtkIdType id;
                            added= mergePoints->InsertUniquePoint(p0,id);
			    if(added)
                                vtkDebugMacro(<< "Added point: " << p0[0] << ", " << p0[1] << ", " << p0[2]);
                            }
			////only add one point -> therefore else if!
                        else if(point_in_both_cells(p1, cell0, cell1, this->DistanceTolerance)){
                            vtkIdType id;
                            added= mergePoints->InsertUniquePoint(p1,id);
			    if(added)
				vtkDebugMacro(<< "Added point: " << p1[0] << ", " << p1[1] << ", " << p1[2]);
                            }//if(point_in_both_cells
                        }//l
                    }//k


                ////if we have at least 3 points: create polygon-cell -> add to out cell list
                if(point_pd->GetPoints()->GetNumberOfPoints() >= 3){
                    
                    //////points have to lie in xy-plane to consistently perform a correct delaunay or convex-hull
                    ////calc N from points
                    
                    double n[3];
                    double angle;
                    double N[3];
                    double x0[3];
                    double x1[3];
                    double x2[3];

		    ////Ids 0-2 given, as point_pd->GetPoints()->GetNumberOfPoints() >= 3
                    point_pd->GetPoints()->GetPoint(0, x0);
                    point_pd->GetPoints()->GetPoint(1, x1);
                    point_pd->GetPoints()->GetPoint(2, x2);

                    z_normal_of_3points(x0, x1, x2, n, angle, N);
		    vtkDebugMacro(<< "N: " << N[0] << ", " << N[1] << ", " << N[2] << "; angle: " << angle);
  

		    if(this->MeshMode== VTK_USE_DELAUNAY2D){

                        ////rotate flat polygon to lie in xy-plane
                        vtkSmartPointer<vtkTransform> xfd0= vtkSmartPointer<vtkTransform>::New();
                        xfd0->RotateWXYZ(angle, N); //angle expected in deg!!!!!!!


                        vtkSmartPointer<vtkDelaunay2D> delaunay2D= vtkSmartPointer<vtkDelaunay2D>::New();
                        delaunay2D->SetInputData(point_pd);

                        ////since vtkDelaunay2D executes in xy-plane a transform of the points (that reside within a single plane) is needed
                        delaunay2D->SetTransform(xfd0);
                        ////Graphics/vtkDelaunay2D.h: VTK_DELAUNAY_XY_PLANE 0, VTK_SET_TRANSFORM_PLANE 1, VTK_BEST_FITTING_PLANE 2
                        delaunay2D->SetProjectionPlaneMode(VTK_DELAUNAY_XY_PLANE); //use VTK_DELAUNAY_XY_PLANE as the data is transformed to lie in the xy-plane
                        delaunay2D->SetTolerance(0.0);
                        delaunay2D->SetAlpha(0);
                        delaunay2D->BoundingTriangulationOff();
                        delaunay2D->Update();
		      
                        vtkSmartPointer<vtkPolyData> out_mesh= vtkSmartPointer<vtkPolyData>::New(); 
                        out_mesh= delaunay2D->GetOutput();

                        appendPD->AddInputData(out_mesh);

                        vtkDebugMacro(<< "Created delaunay-2d mesh with: " << out_mesh->GetNumberOfPoints() << " points and " << out_mesh->GetNumberOfCells() << " cells.");
                        }
		    else{//!useDelaunay2D

                        ////rotate flat polygon to lie in xy-plane, translate first in a way that one point lies in the origin
                        double z0[3];
                        point_pd->GetPoints()->GetPoint(0, z0);

                        /////first translate, then rotate
                        vtkSmartPointer<vtkTransform> xfd0= vtkSmartPointer<vtkTransform>::New();
                        xfd0->RotateWXYZ(angle, N); //angle expected in deg!!!!!!!
                        xfd0->Translate(-z0[0], -z0[1], -z0[2]);

                        vtkSmartPointer<vtkTransformPolyDataFilter> xform0= vtkSmartPointer<vtkTransformPolyDataFilter>::New();
                        xform0->SetInputData(point_pd);
                        xform0->SetTransform(xfd0);
                        xform0->Update();

                        vtkSmartPointer<vtkConvexHull2D> hull2D = vtkSmartPointer<vtkConvexHull2D>::New();
                        hull2D->SetInputData(xform0->GetOutput());
			hull2D->SetMinHullSizeInWorld(0.0);
			hull2D->OutlineOff();
                        hull2D->Update();

                        ////first rotate back, then translate back
                        vtkSmartPointer<vtkTransform> xfd1= vtkSmartPointer<vtkTransform>::New();
                        xfd1->Translate(z0);
                        xfd1->RotateWXYZ(-angle, N); //angle expected in deg!!!!!!!
		     
                        vtkSmartPointer<vtkTransformPolyDataFilter> xform1= vtkSmartPointer<vtkTransformPolyDataFilter>::New();
                        xform1->SetInputData(hull2D->GetOutput());
                        xform1->SetTransform(xfd1);
                        xform1->Update();

                        vtkSmartPointer<vtkPolyData> out_mesh= vtkSmartPointer<vtkPolyData>::New(); 
                        out_mesh= xform1->GetOutput();

                        appendPD->AddInputData(out_mesh);

                        vtkDebugMacro(<< "Created vtkConvexHull2D polygon with: " << out_mesh->GetNumberOfPoints() << " points and " << out_mesh->GetNumberOfCells() << " cells.");
                        }//useDelaunay2D
                    }
                else if(point_pd->GetPoints()->GetNumberOfPoints() > 0){
                    vtkDebugMacro(<< "Only " << point_pd->GetNumberOfPoints() << " points found. This can happen!");
                    }

                }//coplanar cells
            //// progress | abort
	    counter++;
	    if(!(counter%progressInterval)){
                //printf("\r#c0: %9lld, #c1: %9lld; %6.2f%%", i, j, (counter)*100.0/(N0*N1));
                this->UpdateProgress(static_cast<double>(counter)/N0/N1);
		abortExecute = this->GetAbortExecute();
                }

            }//cell1
        }//cell0

    vtkSmartPointer<vtkGeometryFilter> usg2poly2= vtkSmartPointer<vtkGeometryFilter>::New();
    usg2poly2->SetInputConnection(appendPD->GetOutputPort());

    ////since each polygon is added individually, dublicate points and edges exist, vtkCleanPolyData can be used to clean this up as the mesh should not contain any unused points now
    vtkSmartPointer<vtkCleanPolyData> clean= vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputConnection(usg2poly2->GetOutputPort());
    clean->Update();

    output->ShallowCopy(clean->GetOutput());

    return 1;
    }


int vtkCoplanarSurfaceExtractor::FillInputPortInformation(int, vtkInformation *info){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
    }

void vtkCoplanarSurfaceExtractor::PrintSelf(ostream& os, vtkIndent indent){
    this->Superclass::PrintSelf(os,indent);

    os << indent << "Mesh Mode: " << (this->MeshMode == VTK_USE_DELAUNAY2D ?
        "Using vtkDelaunay2D.\n" : "Using vtkConvexHull2D.\n");
    os << indent << "Distance Tolerance: " << this->DistanceTolerance << "\n";
    os << indent << "Face Orientation Tolerance: " << this->FaceOrientationTolerance << "\n";
    os << indent << "Line Orientation Tolerance: " << this->LineOrientationTolerance << "\n";
    }



double vtkCoplanarSurfaceExtractor::SafeAcos (double x){
    if (x < -1.0){
	x = -1.0;
        vtkWarningMacro(<< "   SafeAcos: x + 1.0 == " << x + 1.0 << " < 0.0");
	}
    else if (x > 1.0){
	x = 1.0;
        vtkWarningMacro(<< "   SafeAcos: x - 1.0 == " << x - 1.0 << " > 0.0");
	}
    return acos (x) ;
    }

// int vtkCoplanarSurfaceExtractor::are_coincident(double x[3], double y[3]){

//     if((x[0]==y[0])&&(x[1]==y[1])&&(x[2]==y[2]))
//         return 1;
//     else
//         return 0;
//     }

// double vtkCoplanarSurfaceExtractor::angle_in_deg(double a[3], double b[3]){
//     vtkMath::Normalize(a);
//     vtkMath::Normalize(b);
//     return(acos(vtkMath::Dot(a, b))*180/vtkMath::Pi());
//     }

// double vtkCoplanarSurfaceExtractor::angle_in_deg(double a0[3], double a1[3], double b0[3], double b1[3]){
//     double a[3];
//     double b[3];    
//     vtkMath::Subtract(a0, a1, a);
//     vtkMath::Subtract(b0, b1, b);
//     vtkMath::Normalize(a);
//     vtkMath::Normalize(b);
//     return(acos(vtkMath::Dot(a, b))*180/vtkMath::Pi());
//     }

// double vtkCoplanarSurfaceExtractor::angle_in_deg_unori(double a[3], double b[3]){
//     vtkMath::Normalize(a);
//     vtkMath::Normalize(b);
//     double t=acos(vtkMath::Dot(a, b))*180/vtkMath::Pi();
//     if(t >= 90)
//         return(180-t);
//     else
//         return(t);
//     }

double vtkCoplanarSurfaceExtractor::angle_in_deg_unori(double a0[3], double a1[3], double b0[3], double b1[3]){
    double a[3];
    double b[3];    
    vtkMath::Subtract(a0, a1, a);
    vtkMath::Subtract(b0, b1, b);
    vtkMath::Normalize(a);
    vtkMath::Normalize(b);
    double t=acos(vtkMath::Dot(a, b))*180/vtkMath::Pi();
    if(t >= 90)
        return(180-t);
    else
        return(t);
    }

int vtkCoplanarSurfaceExtractor::coplanar_check(double n0[3], double n1[3], double p0[3], double p1[3], double ft, double dt){
    //////(planar) polygons are coplanar, if the containing planes are identical, 
    //////i.e. are parallel and have the same displacement from origin

    ////check if cells are nearly parallel using face-tollerance

    // vtkMath::Normalize(n0); //see note below
    // vtkMath::Normalize(n1); //see note below
    double t= SafeAcos(vtkMath::Dot(n0, n1))*180/vtkMath::Pi();
    if(t >= 90)
        t=180-t;
    ////normalization commented out as it gives a 4x speed-up
    ////this can be done here because the inputs are only normals from vtkPolyDataNormals
    ////which employs vtkPolygon::ComputeNormal->vtkTriangle::ComputeNormal->vtkTriangle.h which normalizes with double precision
    ////HOWEVER, the normalization is different and without vtkMath::Normalize values just outside [-1.0; +1.0] can occure that are not even visible with %e after substraction of +/-1.0
    ////THERFORE, SafeAcos was introduced, http://stackoverflow.com/questions/8489792/is-it-legal-to-take-acos-of-1-0f-or-1-0f
 

    ////check if containing planes have similar displacement from origin using displacement-tolerance
    ////any point in the plane projected on the plane-normal will yield the plane-origin displacement

    double dd= vtkMath::Dot(n0, p0) - vtkMath::Dot(n1, p1);
     
    if((t < ft) && (dd < dt))
	return 1;
    else
	return 0;
    }

int vtkCoplanarSurfaceExtractor::point_in_both_cells(double p0[3], vtkCell *cell0, vtkCell *cell1, double d_tol){
    double closestPoint[3];
    int subId;
    double  	pcoords[3];
    double dist0, dist1;
    double weights[3];
  
    //check also vtkCellLocator
    int res0= cell0->EvaluatePosition(p0, closestPoint, subId, pcoords, dist0, weights);     
    int res1= cell1->EvaluatePosition(p0, closestPoint, subId, pcoords, dist1, weights);     

    if (res0 < 0){
        vtkErrorMacro(<< "computational problem encountered");
        exit(1);
        }
    if (res1 < 0){
        vtkErrorMacro(<< "computational problem encountered");
        exit(1);
        }

    if((dist0 < d_tol) && (dist1 < d_tol)){
        return 1;
        }
    else
        return 0;
    }

void vtkCoplanarSurfaceExtractor::z_normal_of_3points(double a0[3], double a1[3], double a2[3], double n[3], double &angle, double N[3]){
    double a[3];
    double b[3];
    const double z[3]={0.0,0.0,1.0}; //rotation into xy-plane

    vtkMath::Subtract(a0, a1, a);
    vtkMath::Subtract(a0, a2, b);
    vtkMath::Cross(a, b, n);
    vtkMath::Normalize(n);
    //vtkMath::Normalize(z);//not necessary for {0.0,0.0,1.0}
    vtkMath::Cross(n, z, N);
    //vtkMath::Normalize(N);//should not be necessary

    angle= acos(vtkMath::Dot(n, z))*180/vtkMath::Pi(); //RotateWXYZ(angle, N) expects angle expected in deg!!!!!!!
    
    return;
    }


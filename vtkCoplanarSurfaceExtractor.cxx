
//////boolean intersection operation on coplanar regions of two polydata meshes
//////idea and implementation by Roman Grothausmann

////vtkDistancePolyDataFilter and vtkPartialVolumeModeller taken as examples

////the filters own headers
#include "vtkCoplanarSurfaceExtractor.h"

#include <vtkSmartPointer.h>

////vtk-filter headers
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>

////headers for this specific filter
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
//#include <vtkPolygon.h>
#include <vtkGeometryFilter.h>


vtkStandardNewMacro(vtkCoplanarSurfaceExtractor);




vtkCoplanarSurfaceExtractor::vtkCoplanarSurfaceExtractor(){

    this->SetNumberOfInputPorts(2);

    this->DistanceTolerance= 0.001;
    this->FaceOrientationTolerance= 0.001;
    this->LineOrientationTolerance= this->FaceOrientationTolerance; //for now let face-tol == line-tol

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



    if(this->MeshMode== VTK_USE_DELAUNAY2D)
	std::cerr << "Using vtkDelaunay2D." << std::endl;
    else
	std::cerr << "Using vtkConvexHull2D." << std::endl; 


    int verbose;
    //verbose=1;


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
    if(verbose)
        std::cerr << "Mesh 1 contains " << N0 << " cells and " << mesh0->GetNumberOfPoints() << " points."<< std::endl; 
    vtkIdType N1= mesh1->GetNumberOfCells();
    if(verbose)
        std::cerr << "Mesh 2 contains " << N1 << " cells and " << mesh1->GetNumberOfPoints() << " points."<< std::endl; 

    double x[3];
    double y[3];

    vtkDataArray *normals0= mesh0->GetCellData()->GetNormals();
    vtkDataArray *normals1= mesh1->GetCellData()->GetNormals();

    if(normals0){
        if(verbose)
            std::cerr << "There are " << normals0->GetNumberOfTuples() << " normals in normals0." << std::endl;
        }
    else {
        vtkErrorMacro(<< "No normals in normals0! Aborting.");
        exit(1);
        }
        
    if(normals1){
        if(verbose)
            std::cerr << "There are " << normals1->GetNumberOfTuples() << " normals in normals1." << std::endl;
        }
    else {
        vtkErrorMacro(<< "No normals in normals1! Aborting.");
        exit(1);
        }
        


    vtkSmartPointer<vtkAppendPolyData> appendPD= vtkSmartPointer<vtkAppendPolyData>::New();
    ////needed for disjunct sets??? test!!!
    // vtkSmartPointer<vtkPolyData> empty_mesh= vtkSmartPointer<vtkPolyData>::New();
    // appendPD->AddInputData(empty_mesh);

    vtkIdType counter= 0;

    for (vtkIdType i= 0; i < N0; i++){

        if(verbose)
            std::cerr << "i: " << i << std::endl;
        
        vtkCell* cell0= mesh0->GetCell(i);
        if (cell0->GetCellType() != VTK_TRIANGLE){
            if(verbose)
                std::cerr << "Cell " << i << " is no triangle (" << vtkCellTypes::GetClassNameFromTypeId(cell0->GetCellType()) << ")." << std::endl;
            continue;
            }

        for (vtkIdType j= 0; j < N1; j++){

            if(verbose)
                std::cerr << "j: " << j << std::endl;
        
            vtkCell* cell1= mesh1->GetCell(j);
            if (cell1->GetCellType() != VTK_TRIANGLE){
                if(verbose)
                    std::cerr << "Cell " << j << " is no triangle (" << vtkCellTypes::GetClassNameFromTypeId(cell1->GetCellType()) << ")." << std::endl;
                continue;
                }

	    ////check if cells are nearly parallel using face-tollerance (this->FaceOrientationTolerance)
            double a=angle_in_deg_unori(normals0->GetTuple(i), normals1->GetTuple(j));
            if(verbose)
                std::cerr << "Angle: " << a << std::endl;
            if (a <= this->FaceOrientationTolerance){

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
		mergePoints->SetDivisions(10,10,10);
		//mergePoints->InitPointInsertion(point_pd->GetPoints(), point_pd->GetBounds());
		mergePoints->InitPointInsertion(point_pd->GetPoints(), mesh0->GetBounds()); //either mesh0->GetBounds() or mesh1->GetBounds() should do as no points from outside the intersection of mesh0 and mesh1 should appear!


                if(verbose)
                    std::cerr << "Checking intersections... " << std::endl;
                
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
                            if(verbose)
                                std::cerr << "Angle: " << b << " No intersection considered!" << std::endl;
                            continue;
                            }

                        ////get closest points between lines
                        double d= vtkLine::DistanceBetweenLineSegments(l00, l01, l10, l11, p0, p1, t1, t2);//return value is the shortest distance squared between the two line-segments.
                        if(verbose)
                            std::cerr << "Distance: " << d << std::endl;

			int added;
                        if(point_in_both_cells(p0, cell0, cell1, this->DistanceTolerance)){
                            vtkIdType id;
                            added= mergePoints->InsertUniquePoint(p0,id);
			    if(added)
                                if(verbose)
                                    std::cerr << "Added point: " << p0[0] << ", " << p0[1] << ", " << p0[2]  << std::endl;
                            }
			////only add one point -> therefore else if!
                        else if(point_in_both_cells(p1, cell0, cell1, this->DistanceTolerance)){
                            vtkIdType id;
                            added= mergePoints->InsertUniquePoint(p1,id);
			    if(added)
                                if(verbose)
                                    std::cerr << "Added point: " << p1[0] << ", " << p1[1] << ", " << p1[2]  << std::endl;
                            }//if(point_in_both_cells
                        }//l
                    }//k


                ////if we have at least 3 points: create polygon-cell -> add to out cell list
                if(point_pd->GetPoints()->GetNumberOfPoints() >= 3){
                    
                    //////points have to lie in xy-plane to consistently perform a correct delaunay
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
		    if(verbose)
                        std::cerr << "N: " << N[0] << ", " << N[1] << ", " << N[2] << "; angle: " << angle << std::endl;
  

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

                        if(verbose)
                            std::cerr << "Created delaunay-2d mesh with: " << out_mesh->GetNumberOfPoints() << " points and " << out_mesh->GetNumberOfCells() << " cells." << std::endl;
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

                        if(verbose)
                            std::cerr << "Created vtkConvexHull2D polygon with: " << out_mesh->GetNumberOfPoints() << " points and " << out_mesh->GetNumberOfCells() << " cells." << std::endl;
                        }//useDelaunay2D
                    }
                else if(point_pd->GetPoints()->GetNumberOfPoints() > 0){
                    if(verbose)
                        std::cerr << "Only " << point_pd->GetNumberOfPoints() << " points found. This can happen!" << std::endl;
                    }

                }//coplanar cells
	    counter++;
	    if(!(j%10000))
                printf("\r#c0: %9lld, #c1: %9lld; %6.2f%%", i, j, (counter)*100.0/(N0*N1));
            }//cell1
        }//cell0
    printf("\n");
    vtkSmartPointer<vtkGeometryFilter> usg2poly2= vtkSmartPointer<vtkGeometryFilter>::New();
    usg2poly2->SetInputConnection(appendPD->GetOutputPort());

    ////since each polygon is added individually, dublicate points and edges exist, vtkCleanPolyData can be used to clean this up as the mesh should not contain any unused points now
    vtkSmartPointer<vtkCleanPolyData> clean= vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputConnection(usg2poly2->GetOutputPort());
    clean->Update();

    output->ShallowCopy(clean->GetOutput());

    return 1;
    }


int vtkCoplanarSurfaceExtractor::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

void vtkCoplanarSurfaceExtractor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Mesh Mode: " << (this->MeshMode == VTK_USE_DELAUNAY2D ?
                                       "Using vtkDelaunay2D.\n" : "Using vtkConvexHull2D.\n");
  os << indent << "Distance Tolerance: " << this->DistanceTolerance << "\n";
  os << indent << "Face Orientation Tolerance: " << this->FaceOrientationTolerance << "\n";
  os << indent << "Line Orientation Tolerance: " << this->LineOrientationTolerance << "\n";
}



int vtkCoplanarSurfaceExtractor::are_coincident(double x[3], double y[3]){

    if((x[0]==y[0])&&(x[1]==y[1])&&(x[2]==y[2]))
        return 1;
    else
        return 0;
    }

double vtkCoplanarSurfaceExtractor::angle_in_deg(double a[3], double b[3]){
    //std::cerr << "n0: " << a[0] << ", " << a[1] << ", " << a[2] << "; n1: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
    vtkMath::Normalize(a);
    vtkMath::Normalize(b);
    return(acos(vtkMath::Dot(a, b))*180/vtkMath::Pi());
    }

double vtkCoplanarSurfaceExtractor::angle_in_deg(double a0[3], double a1[3], double b0[3], double b1[3]){
    //std::cerr << "n0: " << a[0] << ", " << a[1] << ", " << a[2] << "; n1: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
    double a[3];
    double b[3];    
    vtkMath::Subtract(a0, a1, a);
    vtkMath::Subtract(b0, b1, b);
    vtkMath::Normalize(a);
    vtkMath::Normalize(b);
    return(acos(vtkMath::Dot(a, b))*180/vtkMath::Pi());
    }

double vtkCoplanarSurfaceExtractor::angle_in_deg_unori(double a[3], double b[3]){
    //std::cerr << "n0: " << a[0] << ", " << a[1] << ", " << a[2] << "; n1: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
    vtkMath::Normalize(a);
    vtkMath::Normalize(b);
    double t=acos(vtkMath::Dot(a, b))*180/vtkMath::Pi();
    if(t >= 90)
        return(180-t);
    else
        return(t);
    }

double vtkCoplanarSurfaceExtractor::angle_in_deg_unori(double a0[3], double a1[3], double b0[3], double b1[3]){
    //std::cerr << "n0: " << a[0] << ", " << a[1] << ", " << a[2] << "; n1: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
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
    std::cerr << "computational problem encountered" << std::endl;
    exit(1);
  }
  if (res1 < 0){
    std::cerr << "computational problem encountered" << std::endl;
    exit(1);
  }

  //if((res0 > 0) && (res1 > 0)){
  ////or make return val dep on dist2 with a tolarance
  if((dist0 < d_tol) && (dist1 < d_tol)){
    return 1;
  }
  else
    return 0;
}



void vtkCoplanarSurfaceExtractor::z_normal_of_3points(double a0[3], double a1[3], double a2[3], double n[3], double &angle, double N[3]){
    double a[3];
    double b[3];
    //double n[3];
    double z[3]={0,0,1};

    vtkMath::Subtract(a0, a1, a);
    vtkMath::Subtract(a0, a2, b);
    vtkMath::Cross(a, b, n);
    vtkMath::Cross(n, z, N);
    vtkMath::Normalize(n);
    vtkMath::Normalize(N);//should not be necessary, just to be on the save side
    vtkMath::Normalize(z);

    //angle= acos(vtkMath::Dot(n, z));
    angle= acos(vtkMath::Dot(n, z))*180/vtkMath::Pi(); //RotateWXYZ(angle, N) expects angle expected in deg!!!!!!!
    
    // std::cerr << "a: " << a[0] << ", " << a[1] << ", " << a[2] << "; b: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
    //  std::cerr << "n: " << n[0] << ", " << n[1] << ", " << n[2] << std::endl;
    // std::cerr << "N: " << N[0] << ", " << N[1] << ", " << N[2] << "; angle: " << angle << std::endl;

    return;
    }


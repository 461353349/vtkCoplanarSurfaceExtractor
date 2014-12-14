/////program to do boolean intersection operation on coplanar regions of two polydata meshes using vtkCoplanarSurfaceExtractor

#include <vtkSmartPointer.h>


#include <vtkXMLPolyDataReader.h>

//for ProgressFunction(
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>


#include "vtkCoplanarSurfaceExtractor.h"
#include "testing/vtkHausdorffDistancePointSetFilter.h"
#include <vtkFieldData.h>



#define P_VERBOSE 1

#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


void ProgressFunction( vtkObject* caller, long unsigned int eventId, void* clientData, void* callData ){
  vtkAlgorithm *d= static_cast<vtkAlgorithm*>(caller);
  if(P_VERBOSE) fprintf(stderr, "\rFilter progress: %5.1f%%", 100.0 * d->GetProgress());
  if(P_VERBOSE) std::cerr.flush(); //not needed for cerr?!
}



int main(int argc, char* argv[]){


    if( argc != 8 )
	{
	std::cerr << "Usage: " << argv[0];
	std::cerr << " inputMesh0";
	std::cerr << " inputMesh1";
	std::cerr << " testMesh";
	std::cerr << " d-tolerance a-tolerance";
	std::cerr << " vtkConvexHull2D=0|vtkDelaunay2D=1";
	std::cerr << " testing-tolerance";
	std::cerr << std::endl;  
	return EXIT_FAILURE;
	}

    if(!(strcasestr(argv[1],".vtp"))) {
	std::cerr << "The input should end with .vtp" << std::endl; 
	return -1;
	}

    if(!(strcasestr(argv[2],".vtp"))) {
	std::cerr << "The input should end with .vtp" << std::endl; 
	return -1;
	}

    if(!(strcasestr(argv[3],".vtp"))) {
	std::cerr << "The output should end with .vtp" << std::endl; 
	return -1;
	}


    VTK_CREATE(vtkCallbackCommand, progressCallback);
    progressCallback->SetCallback(ProgressFunction);


    double d_tol= atof(argv[4]);
    double f_tol= atof(argv[5]);
    double l_tol= f_tol; //for now let face-tol == line-tol

    //vtkObject::SetGlobalWarningDisplay(1);

    vtkSmartPointer<vtkXMLPolyDataReader> reader0 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
 
    reader0->SetFileName(argv[1]);
    reader0->Update();

    vtkSmartPointer<vtkXMLPolyDataReader> reader1 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
 
    reader1->SetFileName(argv[2]);
    reader1->Update();

    vtkSmartPointer<vtkCoplanarSurfaceExtractor> filter = vtkSmartPointer<vtkCoplanarSurfaceExtractor>::New(); 

    filter->SetInputConnection(0, reader0->GetOutputPort());
    filter->SetInputConnection(1, reader1->GetOutputPort());


    filter->SetDistanceTolerance(atof(argv[4]));
    filter->SetFaceOrientationTolerance(atof(argv[5]));
    filter->SetLineOrientationTolerance(filter->GetFaceOrientationTolerance());

    if(atoi(argv[6]))
        filter->SetMeshMode(VTK_USE_DELAUNAY2D);
    else
        filter->SetMeshMode(VTK_USE_CONVEXHULL2D);

    filter->Update();


    ///////////// testing ///////////
    double tol= atof(argv[7]);
    bool tresult= true;

    //////load dataset for comparison
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(argv[3]);
    reader->Update();

    //////compare filter output and test dataset
    vtkSmartPointer<vtkHausdorffDistancePointSetFilter> HDdistance = vtkSmartPointer<vtkHausdorffDistancePointSetFilter>::New();
    HDdistance->SetInputConnection(0, reader->GetOutputPort());
    HDdistance->SetInputConnection(1, filter->GetOutputPort());
    HDdistance->SetTargetDistanceMethod( vtkHausdorffDistancePointSetFilter::POINT_TO_POINT );
    HDdistance->Update();
 
    printf("p2p: HausdorffDistance: %e; RelativeDistanceAtoB: %e; RelativeDistanceBtoA: %e\n", 
	static_cast<vtkPointSet*>(HDdistance->GetOutput(0))->GetFieldData()->GetArray("HausdorffDistance")->GetComponent(0,0),
	HDdistance->GetOutputDataObject(0)->GetFieldData()->GetArray("RelativeDistanceAtoB")->GetComponent(0,0),
	HDdistance->GetOutputDataObject(1)->GetFieldData()->GetArray("RelativeDistanceBtoA")->GetComponent(0,0)
	);

    // vtkSmartPointer<vtkXMLPolyDataWriter> Pwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    // Pwriter->SetFileName("HD_out_p2p.vtp");
    // Pwriter->SetInputConnection(HDdistance->GetOutputPort());
    // Pwriter->Update();

    if(tol < static_cast<vtkPointSet*>(HDdistance->GetOutput(0))->GetFieldData()->GetArray("HausdorffDistance")->GetComponent(0,0))
	tresult= false;
    if(tol < HDdistance->GetOutputDataObject(0)->GetFieldData()->GetArray("RelativeDistanceAtoB")->GetComponent(0,0))
	tresult= false;	
    if(tol < HDdistance->GetOutputDataObject(1)->GetFieldData()->GetArray("RelativeDistanceBtoA")->GetComponent(0,0))
	tresult= false;

    HDdistance->SetTargetDistanceMethod( vtkHausdorffDistancePointSetFilter::POINT_TO_CELL );
    HDdistance->Update();

    printf("p2c: HausdorffDistance: %e; RelativeDistanceAtoB: %e; RelativeDistanceBtoA: %e\n", 
	static_cast<vtkPointSet*>(HDdistance->GetOutput(0))->GetFieldData()->GetArray("HausdorffDistance")->GetComponent(0,0),
	HDdistance->GetOutputDataObject(0)->GetFieldData()->GetArray("RelativeDistanceAtoB")->GetComponent(0,0),
	HDdistance->GetOutputDataObject(1)->GetFieldData()->GetArray("RelativeDistanceBtoA")->GetComponent(0,0)
	);

    if(tol < static_cast<vtkPointSet*>(HDdistance->GetOutput(0))->GetFieldData()->GetArray("HausdorffDistance")->GetComponent(0,0))
	tresult= false;
    if(tol < HDdistance->GetOutputDataObject(0)->GetFieldData()->GetArray("RelativeDistanceAtoB")->GetComponent(0,0))
	tresult= false;
    if(tol < HDdistance->GetOutputDataObject(1)->GetFieldData()->GetArray("RelativeDistanceBtoA")->GetComponent(0,0))
	tresult= false;


    if(tresult)
	return EXIT_SUCCESS;
    else
	return EXIT_FAILURE;

    }

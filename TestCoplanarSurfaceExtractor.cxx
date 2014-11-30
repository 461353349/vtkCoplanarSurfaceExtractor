/////program to do boolean intersection operation on coplanar regions of two polydata meshes using vtkCoplanarSurfaceExtractor

#include <vtkSmartPointer.h>


#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>


#include "vtkCoplanarSurfaceExtractor.h"


// class vtkProgressCommand : public vtkCommand
// {
// public:
//   static vtkProgressCommand *New(){
//     return new vtkProgressCommand;
//   }

//   virtual void Execute(vtkObject *caller, unsigned long, void *callData){
//     double progress = *(static_cast<double*>(callData));
//     fprintf(stderr, "\rFilter progress: %5.1f\n", 100.0 * progress);
//     std::cerr.flush();
//   }
// };


int main(int argc, char* argv[]){


    if( argc != 7 )
	{
	std::cerr << "Usage: " << argv[0];
	std::cerr << " inputMesh0";
	std::cerr << " inputMesh1";
	std::cerr << " outputMesh";
	std::cerr << " d-tolarance a-tolarance";
	std::cerr << " vtkConvexHull2D=0|vtkDelaunay2D=1";
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

    double d_tol= atof(argv[4]);
    double f_tol= atof(argv[5]);
    double l_tol= f_tol; //for now let face-tol == line-tol

    int useDelaunay2D= atoi(argv[6]);


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


    filter->Print(std::cerr);
    filter->Update();
    
    vtkSmartPointer<vtkXMLPolyDataWriter> Pwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
 
    Pwriter->SetFileName(argv[3]);
    Pwriter->SetInputConnection(filter->GetOutputPort());
    Pwriter->Update();




    return EXIT_SUCCESS;
    }

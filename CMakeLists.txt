CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

#Change PROJECT_NAME to the name of your project
PROJECT(CoplanarSurfaceExtractor)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(CTest)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
#
# The packages can be specified with a version number, for example:
#
# ITK 2.8.1
# ITK 3.2.0
#
# If no version is specified, the most recent release of the package
# will be used.
SET(Required_Packages
  VTK
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
  ${VTK_LIBRARIES}
)


IF(BUILD_TESTING)

#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above.
#You can build more than one executable per project
SET(CurrentExe "TestCoplanarSurfaceExtractor")
ADD_EXECUTABLE(${CurrentExe}
  vtkCoplanarSurfaceExtractor.cxx testing/vtkHausdorffDistancePointSetFilter.cxx
  TestCoplanarSurfaceExtractor.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.

##comparison data sets created with CoplanarSurfaceExtractor + Clean:

ADD_TEST(NAME triangles_inside_01_0001-0002_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_inside_01_0001.vtp triangles_inside_01_0002.vtp triangles_inside_01_0001-0002.vtp 0.002 0.009 0 0.0001)

ADD_TEST(NAME triangles_inside_01_0001-0002_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_inside_01_0001.vtp triangles_inside_01_0002.vtp triangles_inside_01_0001-0002.vtp 0.002 0.009 1 0.0001)

ADD_TEST(NAME triangles_inside_02_0001-0002_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_inside_02_0001.vtp triangles_inside_02_0002.vtp triangles_inside_02_0001-0002.vtp 0.002 0.02 0 0.0001)

ADD_TEST(NAME triangles_inside_02_0001-0002_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_inside_02_0001.vtp triangles_inside_02_0002.vtp triangles_inside_02_0001-0002.vtp 0.002 0.02 1 0.0001)

ADD_TEST(NAME triangles_corner-on-corner_02_0001-0002_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_corner-on-corner_02_0001.vtp triangles_corner-on-corner_02_0002.vtp triangles_corner-on-corner_02_0001-0002.vtp 0.000000001 0.000000001 0 0.0001)

ADD_TEST(NAME triangles_corner-on-corner_02_0001-0002_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_corner-on-corner_02_0001.vtp triangles_corner-on-corner_02_0002.vtp triangles_corner-on-corner_02_0001-0002.vtp 0.000000001 0.000000001 1 0.0001)

ADD_TEST(NAME triangles_6intersections_02_0001-0002_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_6intersections_02_0001.vtp triangles_6intersections_02_0002.vtp triangles_6intersections_02_0001-0002.vtp 0.000001 0.007 0 0.0001)

ADD_TEST(NAME triangles_6intersections_02_0001-0002_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_6intersections_02_0001.vtp triangles_6intersections_02_0002.vtp triangles_6intersections_02_0001-0002.vtp 0.000001 0.007 1 0.0001)

ADD_TEST(NAME triangles_overlaping_0001-0002_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_overlaping_0001.vtp triangles_overlaping_0002.vtp triangles_overlaping_0001-0002.vtp 0.000000001 0.000000001 0 0.0001)

ADD_TEST(NAME triangles_overlaping_0001-0002_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_overlaping_0001.vtp triangles_overlaping_0002.vtp triangles_overlaping_0001-0002.vtp 0.000000001 0.000000001 1 0.0001)

ADD_TEST(NAME triangles_overlaping_03_0001-0002_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_overlaping_03_0001.vtp triangles_overlaping_03_0002.vtp triangles_overlaping_03_0001-0002.vtp 0.0000001 0.000000001 0 0.0001)

ADD_TEST(NAME triangles_overlaping_03_0001-0002_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_overlaping_03_0001.vtp triangles_overlaping_03_0002.vtp triangles_overlaping_03_0001-0002.vtp 0.0000001 0.000000001 1 0.0001)

ADD_TEST(NAME triangles_overlaping_3D_03_0001-0003_success_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_overlaping_3D_03_0001.vtp triangles_overlaping_3D_03_0003.vtp triangles_overlaping_3D_03_0001-0003_right.vtp 0.0001 0.0001 0 0.0001)

ADD_TEST(NAME triangles_overlaping_3D_03_0001-0003_fail_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_overlaping_3D_03_0001.vtp triangles_overlaping_3D_03_0003.vtp triangles_overlaping_3D_03_0001-0003_wrong.vtp 0.0001 0.0001 0 0.0001)
SET_TESTS_PROPERTIES(triangles_overlaping_3D_03_0001-0003_fail_CH PROPERTIES WILL_FAIL TRUE)

ADD_TEST(NAME triangles_overlaping_3D_03_0001-0003_success_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_overlaping_3D_03_0001.vtp triangles_overlaping_3D_03_0003.vtp triangles_overlaping_3D_03_0001-0003_right.vtp 0.0001 0.0001 1 0.0001)

ADD_TEST(NAME triangles_overlaping_3D_03_0001-0003_fail_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/triangle_tests/ COMMAND TestCoplanarSurfaceExtractor triangles_overlaping_3D_03_0001.vtp triangles_overlaping_3D_03_0003.vtp triangles_overlaping_3D_03_0001-0003_wrong.vtp 0.0001 0.0001 1 0.0001)
SET_TESTS_PROPERTIES(triangles_overlaping_3D_03_0001-0003_fail_DL PROPERTIES WILL_FAIL TRUE)

ADD_TEST(NAME particle_1335-1474_success_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/mc-surf_tests/ COMMAND TestCoplanarSurfaceExtractor particle_1335_raw.vtp particle_1474_raw.vtp particle_1335-1474.vtp 0.000000001 0.02 0 0.0001)

ADD_TEST(NAME particle_1335-1474_fail_CH WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/mc-surf_tests/ COMMAND TestCoplanarSurfaceExtractor particle_1335_raw.vtp particle_1474_raw.vtp particle_1335-1474.vtp 0.000000001 0.01 0 0.0001)
SET_TESTS_PROPERTIES(particle_1335-1474_fail_CH PROPERTIES WILL_FAIL TRUE)

ADD_TEST(NAME particle_1335-1474_success_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/mc-surf_tests/ COMMAND TestCoplanarSurfaceExtractor particle_1335_raw.vtp particle_1474_raw.vtp particle_1335-1474.vtp 0.000000001 0.02 1 0.0001)

ADD_TEST(NAME particle_1335-1474_fail_DL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/testing/mc-surf_tests/ COMMAND TestCoplanarSurfaceExtractor particle_1335_raw.vtp particle_1474_raw.vtp particle_1335-1474.vtp 0.000000001 0.01 1 0.0001)
SET_TESTS_PROPERTIES(particle_1335-1474_fail_DL PROPERTIES WILL_FAIL TRUE)

ENDIF(BUILD_TESTING)

OPTION( BUILD_EXAMPLE
        "Build example"
        OFF )

IF( BUILD_EXAMPLE )
# EXAMPLE ---------------------------------------------------------
   
SET(CurrentExe "CoplanarSurfaceExtractor")
ADD_EXECUTABLE(${CurrentExe}
  vtkCoplanarSurfaceExtractor.cxx
  CoplanarSurfaceExtractor.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

#------------------------------------------------------------------
ENDIF( BUILD_EXAMPLE )


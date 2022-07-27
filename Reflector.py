###############################################################################
# Reflector Script
# Author: Nathan van der Wielen
###############################################################################

# Libraries
from abaqus import *
from abaqusConstants import *
import part
import material
import section
import assembly
import step
import amplitude
import regionToolset
import __main__
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import odbAccess
import os
from operator import add

import numpy as np
from numpy import genfromtxt
import math as m
import csv

# Set Backward compatibility to control access to depricated commands
backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)

# Set the work Directory to the Reflector Folder
os.chdir(r"D:\1-Master Thesis\ABAQUS\Reflector")

# Set the Journal Options to Coordinate. Allows the user to look at the replay file and see coordinates and not Mask objects
session.journalOptions.setValues(replayGeometry=COORDINATE,
                                 recoverGeometry=COORDINATE)

# Set the Background color to white
session.graphicsOptions.setValues(backgroundStyle=SOLID, 
                                  backgroundColor='#FFFFFF')

#------------------------------------------------------------------------------
# User Input
# These Variables are used to Control different aspects of the model
#------------------------------------------------------------------------------
Model_Name = "Reflector_Thermal"
Model_Name_Displacement = 'Reflector_Displacement' # Name of the Displacement Model
Part_Name = "Reflector_Part"
Instance_Name = "Reflector_Instance"
Amp_basename = 'Amplitude-'
Load_basename = 'Load-'
Section_Name = '__Reflector_Section__'
Zipper_Section_Name = '__Zipper_Section__'
Step_Name = 'ReflectorLoad'
Viewport_Name = 'Reflector Example'

Case = 'LEO'
DirectoryName = 'D:/ESATAN_TMS/Model/Reflector_SiO/esatan/LEO_T/'
QS_Top_File = DirectoryName+'Top_QS.csv'
QA_Top_File = DirectoryName+'Top_QA.csv'
QE_Top_File = DirectoryName+'Top_QE.csv'
QS_Bottom_File = DirectoryName+'Bottom_QS.csv'
QA_Bottom_File = DirectoryName+'Bottom_QA.csv'
QE_Bottom_File = DirectoryName+'Bottom_QE.csv'
Surface_Area_File = DirectoryName+'Surface_Area.csv'

Job_Name = 'Reflector_'+Case+'_Job_Thermal'
Job_Name_Displacement = 'Reflector_'+Case+'_Job_Displacement' # Name of the Displacement Job
Report_Name_Thermal = Case+'_thermal'
Report_Name_Displacement = Case+'_deflection'

#------------------------------------------------------------------------------
# Variables
#------------------------------------------------------------------------------
P_Height = 0.34  # Height of the Parabola
P_Diameter = 3.0  # Diameter of the Parabola
P_Spline_Division = 50  # Number of divisions in spline for the Reflector
ESATAN_rad = 10  # Number of radial elements
ESATAN_ang = 40  # Number of angular elements
Mesh_Size = 0.1 # Seed Size
Section_Thickness = 0.0004
Hole_Radius = 0.2
P_Spiral_Parameter = 0.06
P_Spiral_Width = 0.01

# Computed Variables
P_Focus = pow(P_Diameter, 2)/(16*P_Height)  # Focal length of the Parabola
Num_Elem = ESATAN_ang*ESATAN_rad  # Number of elements on the reflector
Num_Face = Num_Elem*2


#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
# This function creates a Surface Set in abaqus
def get_faces(model, instance):
    face = ()
    p = mdb.models[model].rootAssembly.instances[instance]
    f = p.faces
    return f


# This function creates a Face Set in Abaqus
def Create_Set_Face(x, y, z, model, part, set_name):
    face = ()
    p = mdb.models[model].parts[part]
    f = p.faces
    myFace = f.findAt((x, y, z),)
    face = face + (f[myFace.index:myFace.index+1], )
    # print("This is a print out of face",face)
    p.Set(faces=face,
          name=set_name)
    # print("This is a print out of myFace",myFace)
    return myFace


# This function creats an Edge Set in Abaqus
def Create_Set_Edge(x, y, z, model, part, set_name):
    edge = ()
    p = mdb.models[model].parts[part]
    e = p.edges
    myEdge = e.findAt((x, y, z), )
    edge = edge + (e[myEdge.index:myEdge.index+1], )
    f = p.Set(edges=edge,
              name=set_name)
    return myEdge

# This function creates a Surface Set in Abaqus
def Create_Set_Surface(x,y,z,model,part,set_name):
    face = ()
    p = mdb.models[model].parts[part]
    s = p.faces
    myFace = s.findAt((x,y,z),)
    face = face + (s[myFace.index:myFace.index+1], )
    p.Surface(side2Faces=face,
              name=set_name)

# This function creates a set of the Edges belonging to the Spiral
def Get_Spiral_Edge(x, y, z, model, part, set_name):
    edge = ()
    p = mdb.models[model].parts[part]
    e = p.edges
    myEdge = e.findAt((x, y, z), )
    edge = edge + (e[myEdge.index:myEdge.index+1], )
    f = p.Set(edges=edge,
              name=set_name)
    return myEdge


# This function creates a Datum Plane relative to a Principal Plane
def Create_Datum_Plane_by_Principal(type_plane, part, model, offset_plane):
    p = mdb.models[model].parts[part]
    myPlane = p.DatumPlaneByPrincipalPlane(principalPlane=type_plane,
                                           offset=offset_plane)
    myID = myPlane.id
    return myID

# This function creates a Section in Abaqus
def Create_Section(model,Section_Name,material_name,thickness) :
    mdb.models[model].HomogeneousShellSection(name=Section_Name,
                                              preIntegrate=OFF,
                                              material=material_name,
                                              thicknessType=UNIFORM,
                                              thickness=Section_Thickness,
                                              thicknessField='',
                                              nodalThicknessField='',
                                              idealization=NO_IDEALIZATION,
                                              poissonDefinition=DEFAULT,
                                              thicknessModulus=None,
                                              temperature=GRADIENT,
                                              useDensity=OFF,
                                              integrationRule=SIMPSON,
                                              numIntPts=5)


# This Function assigns a section to the part
def Assign_Section(model, part, section):
    p = mdb.models[model].parts[part]
    f = p.faces[:]
    region = regionToolset.Region(faces=f)
    p.SectionAssignment(region=region,
                        sectionName=section,
                        offset=0.0,
                        offsetType=MIDDLE_SURFACE,
                        offsetField='',
                        thicknessAssignment=FROM_SECTION)


# This Function Partitions by sketcing and projecting sketch to reflector
def Create_Partition(model, part, id_plane, edge, face, radius, rad, ang, 
                     focus, height):
    rad_array = [0]*rad
    ang_array = [0]*ang
    p = mdb.models[model].parts[part]
    f1, e1, d2 = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=d2[id_plane],
                              sketchUpEdge=edge,
                              sketchPlaneSide=SIDE2,
                              origin=(0.0, 0.0, 0.0))
    s1 = mdb.models[model].ConstrainedSketch(name='__Partition__',
                                             sheetSize=9.36,
                                             gridSpacing=0.23,
                                             transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s1,
                                  filter=COPLANAR_EDGES)
    # for loop to creatsurfe the different angular cuts in the partition
    for i in range(0, ang):
        ang_array[i] = (i)*(360/ang)
        pos = (round(m.cos(m.radians(ang_array[i]))*(radius+0.1), 5),
               round(m.sin(m.radians(ang_array[i]))*(radius+0.1), 5))
        s1.Line(point1=(0.0, 0.0),
                point2=(pos))
    # for loop to create the circle in the partition
    for i in range(0, rad):
        rad_array[i] = m.sqrt((i+1)*(height/rad)*4*focus)
        s1.CircleByCenterPerimeter(center=(0.0, 0.0),
                                   point1=(rad_array[i], 0.0))
    f = p.faces[:]
    f, e, d1 = p.faces, p.edges, p.datums
    p.PartitionFaceBySketchThruAll(sketchPlane=d1[id_plane],
                                   sketchUpEdge=e[1],
                                   faces=f,
                                   sketchPlaneSide=SIDE2,
                                   sketch=s1)
    s1.unsetPrimaryObject()

# This Function creates the mesh of the reflector
def Create_Mesh(model, part, viewport, seed_size):
    p = mdb.models[model].parts[part]
    # Use elem1 and 2 for linear, elem 8 and 6 for quadratic
    elemType1 = mesh.ElemType(elemCode=DS4,
                              elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=DS3,
                              elemLibrary=STANDARD)
    # elemType1 = mesh.ElemType(elemCode=DS8,
    #                           elemLibrary=STANDARD)
    # elemType2 = mesh.ElemType(elemCode=DS6,
    #                           elemLibrary=STANDARD)
    f = p.faces[:]
    pickedRegions =(f, )
    p.setElementType(regions=pickedRegions,
                     elemTypes=(elemType1, elemType2))
    p.seedPart(size=seed_size,
               deviationFactor=0.1,
               minSizeFactor=0.1)
    pickedRegions = f
    p.setMeshControls(regions=pickedRegions,
                      elemShape=QUAD_DOMINATED)
    p.generateMesh()


# This Function creates the sketch of for the Reflector Part
def Create_Reflector_Sketch(sketch, division, focus, diameter):
    X = []  # Initialize X
    Z = []  # Initialize Z
    for i in range(0, P_Spline_Division+1):
        X.append(((P_Diameter/2)/P_Spline_Division)*i)
        Z.append(pow((((diameter / 2) / division) * i), 2)
                 / (4*focus))
    # Construction Line
    mySketch.ConstructionLine(point1=(0.0, -20.0),
                              point2=(0.0, 20.0))
    # Create the spline line to be rotated
    mySketch.Spline(zip(X, Z))


# This Function returns array of face objects in same order as the pos_array
def Find_Face(model, instance, rad, Height, pos_array, id_array):
    a1 = mdb.models[model].rootAssembly
    f1 = a1.instances[instance].faces
    Face_array = []
    for i in range (0,len(id_array)):
        face = []
        for j in range(0,len(id_array[i])):
            face.append(pos_array[id_array[i][j]])
        Face_array.append(f1.findAt(*face))
    return Face_array

# This function creates a list representing the columns that will be read in the ESATAN CSV files
def Create_Data_Selector(Number_Elements):
    data = ()
    for i in range(Number_Elements):
        data = data + (i+1,)
    return data


# Function returns the number of rows from csv file.
# This does not include number of rows of data.
# -4 is added to remove the first 3 and the last row
def Get_Number_Rows(file):
    interim_file = open(file)
    numline = len(interim_file.readlines())
    return numline-4


# Function returns the number of columns
# Excludes the blank column at the end of the csv file
def Get_Number_Columns(file):
    d = ','
    f = open(file, 'r')
    reader = csv.reader(f, delimiter=d)
    next(reader)
    next(reader)
    ncol = len(next(reader))  # Read first line and count columns
    return ncol-1


def Create_Datum_Point(model,coordinate):
    p = mdb.models[model].rootAssembly
    p.DatumPointByCoordinate(coords=coordinate)


def Create_Center_Cut(model, part, edge, hole_rad, datum):
    p = mdb.models[model].parts[part]
    e, d = p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=d[datum],
                              sketchUpEdge=edge,
                              sketchPlaneSide=SIDE1,
                              sketchOrientation=RIGHT,
                              origin=(0.0, 0.0, 0.0))
    s = mdb.models[model].ConstrainedSketch(name='__Center_Cut__',
                                            sheetSize=9.36,
                                            gridSpacing=0.23,
                                            transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s,
                                  filter=COPLANAR_EDGES)
    s.CircleByCenterPerimeter(center=(0.0, 0.0),
                              point1=(hole_rad, 0.0))
    p.CutExtrude(sketchPlane=d[datum],
                 sketchUpEdge=edge,
                 sketchPlaneSide=SIDE1,
                 sketchOrientation=RIGHT,
                 sketch=s, flipExtrudeDirection=ON)
    s.unsetPrimaryObject()
    
    All_Edges = p.edges[:]
    Cut_Edges = []
    Cut_Edges_ID = []
    # for loop that selects the edges created by Partition face-2
    for i in All_Edges:
        if i.featureName == 'Cut extrude-1':
            Cut_Edges.append(i.pointOn)
            Cut_Edges_ID.append(i.index)
    Cut_Edges_Final = p.edges.findAt(*Cut_Edges)  
    p.Set(edges=Cut_Edges_Final, name='Cut_Edge')
    
    
def cart2sph(x,y): 
    r = m.sqrt(m.pow(x, 2) + m.pow(y, 2))               # r
    theta =  abs(m.degrees(m.atan(y/x)))    # theta
    if x>0 and y>0:
        theta=theta
    if x<0 and y>0:
        theta=180-theta
    if x<0 and y<0:
        theta=180+theta
    if x>0 and y<0:
        theta = 360-theta
    return r, theta


def surf_sort(rad,ang,height,focus,surf_coord):
    rad_array = [0]*(rad+1)
    ang_array = [0]*(ang+1)
    sections = [0]*rad*ang
    surf = [0]*len(surf_coord)
    k=0
    surf_lib = []

    # Definition of the radial circle 
    for i in range(1, rad+1):
        rad_array[i] = m.sqrt((i)*(height/rad)*4*focus)
    # Definition of the angular lines
    for i in range(0, ang+1):
        ang_array[i] = (i)*(360/ang)

    for i in range(0, rad):
        for j in range(0, ang):
            sections[k] = (ang_array[j],ang_array[j+1],rad_array[i],rad_array[i+1])
            k+=1                
    for i in range(0, len(surf_coord)):
        surf[i] = (cart2sph(surf_coord[i][0][0], surf_coord[i][0][1]))
    for i in range(0,len(sections)):
        interim = []
        for j in range (0, len(surf)):
            if surf[j][0]>sections[i][2] and surf[j][0]<sections[i][3] and surf[j][1]>sections[i][0] and surf[j][1]<sections[i][1]:
                interim.append(j)
        surf_lib.append(interim)
    return surf_lib


# This funtion returns array of coordinates of center of faces on the reflector
def Get_Centroid_Coordinates(model, part, rad, ang, radius, focus, height):
    rad_array = [0]*rad
    ang_array = [0]*ang
    pos_array = []
    for i in range(0, ang):
        ang_array[i] = i*(360/ang)+((360/ang)/2)
    for i in range(0, rad):
        rad_array[i] = sqrt((i+1)*(height/rad)*4*focus)
    for i in range(0, rad):
        for j in range(0, ang):
            pos = ((round(m.cos(m.radians(ang_array[j]))*rad_array[i], 10),
                    round(pow(rad_array[i], 2)/(4*focus), 10),
                    round(m.sin(m.radians(ang_array[j]))*rad_array[i], 10)), )
            pos_array.append(pos)
            p = mdb.models[model].parts[part]
            p.DatumPointByCoordinate(coords=((round(m.cos(m.radians(ang_array[j]))*rad_array[i], 10),
                            round(pow(rad_array[i], 2)/(4*focus), 10),
                            round(m.sin(m.radians(ang_array[j]))*rad_array[i], 10))))
    return pos_array


#------------------------------------------------------------------------------
# Get Data From CSV File
#------------------------------------------------------------------------------
ncol = Get_Number_Columns(QS_Top_File)
nrow = Get_Number_Rows(QS_Top_File)

Useddata = Create_Data_Selector(ncol-1)

# Get the time from ESATAN CSV file
time = np.genfromtxt(open(QS_Top_File),
                      dtype=float,
                      delimiter=',',
                      skip_header=3,
                      usecols=(0))
# Get the Solar Fluxes for the top of the Reflector from ESATAN CSV file
QS_Top_Data = np.genfromtxt(open(QS_Top_File),
                            dtype=float,
                            delimiter=',',
                            skip_header=3,
                            usecols=Useddata)
# Get the Albedo Fluxes for the top of the Reflector from ESATAN CSV file
QA_Top_Data = np.genfromtxt(open(QA_Top_File),
                            dtype=float,
                            delimiter=',',
                            skip_header=3,
                            usecols=Useddata)
# Get the Earth IR Fluxes for the top of the Reflector from ESATAN CSV file
QE_Top_Data = np.genfromtxt(open(QE_Top_File),
                            dtype=float,
                            delimiter=',',
                            skip_header=3,
                            usecols=Useddata)
# Get the Solar Fluxes for the bottom of the Reflector from ESATAN CSV file
QS_Bottom_Data = np.genfromtxt(open(QS_Bottom_File),
                                dtype=float,
                                delimiter=',',
                                skip_header=3,
                                usecols=Useddata)
# Get the Albedo Fluxes for the bottom of the Reflector from ESATAN CSV file
QA_Bottom_Data = np.genfromtxt(open(QA_Bottom_File),
                                dtype=float,
                                delimiter=',',
                                skip_header=3,
                                usecols=Useddata)
# Get the Earth IR Fluxes for the bottom of the Reflector from ESATAN CSV file
QE_Bottom_Data = np.genfromtxt(open(QE_Bottom_File),
                                dtype=float,
                                delimiter=',',
                                skip_header=3,
                                usecols=Useddata)
# Get the Surface Areas of each face from ESATAN CSV file
Surface_Area_Data = np.genfromtxt(open(Surface_Area_File),
                                dtype=float,
                                delimiter=',',
                                skip_header=3,
                                usecols=Useddata)


Q_Top_Data = (QS_Top_Data + QA_Top_Data + QE_Top_Data)
Q_Bottom_Data = (QS_Bottom_Data + QA_Bottom_Data + QE_Bottom_Data)

CFRP_emissivity = 0.8

# Got interesting results with Q*0.2 and emissivity at 0.25


#------------------------------------------------------------------------------
# Create Model
#------------------------------------------------------------------------------
# This section create the Reflector model, assigns Constants to the model attributes, and create the Viewport
myModel = mdb.Model(name=Model_Name)

# Set the Attribute Values for the absolute zero and the stefan boltzmann constant
myModel.setValues(absoluteZero=0,
                  stefanBoltzmann=5.669E-8)


# Create viewport to display the model and the results of the analysis.
myViewport = session.Viewport(name=Viewport_Name,
                              origin=(20, 20),
                              width=150, height=120)

myViewport.makeCurrent()
myViewport.maximize()


#------------------------------------------------------------------------------
# Create Part
#------------------------------------------------------------------------------
# This Section create the Reflector Part and partitions it to match the ESATAN mesh
mySketch = myModel.ConstrainedSketch(name='ReflectorProfile',
                                     sheetSize=20)

Create_Reflector_Sketch(mySketch, P_Spline_Division, P_Focus, P_Diameter)

# Create a three-dimensional, deformable part.
myReflector = myModel.Part(name=Part_Name,
                           dimensionality=THREE_D,
                           type=DEFORMABLE_BODY)

# Create a revolve of the previously made spline
myReflector.BaseShellRevolve(sketch=mySketch,
                             angle=360.0,
                             flipRevolveDirection=OFF)

# Create a Datum Plane
myID_1 = Create_Datum_Plane_by_Principal(XZPLANE, Part_Name, Model_Name, 0.0)


# Create the Top Face
Top_Face = Create_Set_Face(0.0, 0.0, 0.0, Model_Name, Part_Name, "Ref_Face")

# Create the Top Edge
Top_Edge = Create_Set_Edge(0.0, 0.34, 1.5, Model_Name, Part_Name, "Top_Edge")

# Create the Top Surface
Top_Surface = Create_Set_Surface(0.0, 0.0, 0.0, Model_Name, Part_Name, "Top_Surface")

# Create a Partition
Create_Partition(Model_Name, Part_Name, myID_1, Top_Edge, Top_Face,
                 P_Diameter/2, ESATAN_rad, ESATAN_ang, P_Focus, P_Height)
    
Create_Center_Cut(Model_Name, Part_Name, Top_Edge, Hole_Radius, myID_1)


#------------------------------------------------------------------------------
# Create Material
#------------------------------------------------------------------------------
# Create Reflector Material
myAl = myModel.Material(name='CFRP(M55J)')
myAl.Elastic(table=((338000000000.0, 0.33), ))
myAl.Conductivity(table=((46.97, ), ))
myAl.SpecificHeat(table=((827.05, ), ))
myAl.Density(table=((1623.0, ), ))
myAl.Expansion(table=((-0.0000011, ), ))

#------------------------------------------------------------------------------
# Create Section.
#------------------------------------------------------------------------------
# Create a Section
mySection = Create_Section(Model_Name,
                           Section_Name,
                           'CFRP(M55J)',
                           Section_Thickness)

# Assign the section to the region.
Assign_Section(Model_Name,
               Part_Name,
               Section_Name)


#------------------------------------------------------------------------------
# Create Part Instance.
#------------------------------------------------------------------------------
myAssembly = myModel.rootAssembly
myInstance = myAssembly.Instance(name=Instance_Name,
                                 part=myReflector,
                                 dependent=ON)
myAssembly.rotate(instanceList=(Instance_Name, ),
                  axisPoint=(0.0, 0.0, 0.0),
                  axisDirection=(1.0, 0.0, 0.0),
                  angle=90.0)


#------------------------------------------------------------------------------
# Create TimePoints.
#------------------------------------------------------------------------------
# This section creates a timepoint in ABAQUS
time_list = []
for i in range(0, len(time)-1):
    time_list.append((time[i], ))
time_tuple = tuple(time_list)
myModel.TimePoint(name='TimePoints-1',
                  points=time_tuple)


#------------------------------------------------------------------------------
# Create Step
#------------------------------------------------------------------------------
# Create a step. The time period of the static step is 1.0,
# and the initial incrementation is 0.1; the step is created
# after the initial step.
myModel.HeatTransferStep(name='ReflectorLoad',
                         previous='Initial',
                         timePeriod=time[nrow],
                         description='Heat Transfer',
                         maxNumInc=1000,
                         deltmx=100,
                         amplitude=RAMP)


#------------------------------------------------------------------------------
# Load and Surface Definition
#------------------------------------------------------------------------------

Coordinates = []
Surf_Array_Top = []
Region_Array_Top = []
Surf_Array_Bottom = []
Region_Array_Bottom = []

faces = get_faces(Model_Name,
                  Instance_Name)

surf_Coord = []  # List of all the coordinates if all the faces on Part
for i in faces:
    j = i.pointOn
    surf_Coord.append((j[0],))
    # n= j[0]
    # Create_Datum_Point(Model_Name, (n[0], n[1], n[2]))

# Surf_ids is a 
surf_ids = surf_sort(ESATAN_rad,
                     ESATAN_ang,
                     P_Height,
                     P_Focus,
                     surf_Coord)

fa1 = Find_Face(Model_Name,
                Instance_Name,
                P_Diameter/2,
                P_Height,
                surf_Coord,
                surf_ids)

# For loop used to create the different surfaces using faces
for i in range(0, Num_Elem):
    myAssembly.Surface(side2Faces=fa1[i], name='Surf-'+str(i+1))
    Surf_Array_Top.append((fa1[i], SIDE2), )
    SurfaceRegion_Top = regionToolset.Region(side2Faces=fa1[i])
    Region_Array_Top.append(SurfaceRegion_Top)
    
for i in range(0, Num_Elem):
    myAssembly.Surface(side1Faces=fa1[i], name='Surf-'+str(i+1+Num_Elem))
    Surf_Array_Bottom.append((fa1[i], SIDE1), )
    SurfaceRegion_Bottom = regionToolset.Region(side1Faces=fa1[i])
    Region_Array_Bottom.append(SurfaceRegion_Bottom)

# Conversion of Surfaces and Region arrays to tuples
Surf_Array_Top_tuple = tuple(Surf_Array_Top)
Region_Array_Top_tuple = tuple(Region_Array_Top)

Surf_Array_Bottom_tuple = tuple(Surf_Array_Bottom)
Region_Array_Bottom_tuple = tuple(Region_Array_Bottom)


# Create Amplitudes with surface area in mind
# Create the Top Amplitudes and Loads and assign them to each other
for i in range(0, Num_Elem):
    Amp_name = Amp_basename + str(i+1)
    Load_name = Load_basename + str(i+1)
    myModel.TabularAmplitude(name=Amp_name,
                              data=zip(time[0:nrow, ], Q_Top_Data[:, i]/Surface_Area_Data[:, i]),
                              timeSpan=STEP)
    myModel.SurfaceHeatFlux(name=Load_name,
                            createStepName=Step_Name,
                            region=Region_Array_Top_tuple[i],
                            magnitude=1,
                            amplitude=Amp_name)
    
# Create the Bottom Amplitudes and Loads and assign them to each other     
for i in range(0, Num_Elem):
    Amp_name = Amp_basename + str(i+1+Num_Elem)
    Load_name = Load_basename + str(i+1+Num_Elem)
    myModel.TabularAmplitude(name=Amp_name,
                              data=zip(time[0:nrow, ], Q_Bottom_Data[:, i]/Surface_Area_Data[:, i]),
                              timeSpan=STEP)
    myModel.SurfaceHeatFlux(name=Load_name,
                            createStepName=Step_Name,
                            region=Region_Array_Bottom_tuple[i],
                            magnitude=1,
                            amplitude=Amp_name)


# # Amplitudes Using total flux
# # Create the Top Amplitudes and Loads and assign them to each other
# for i in range(0, Num_Elem):
#     Amp_name = Amp_basename + str(i+1)
#     Load_name = Load_basename + str(i+1)
#     myModel.TabularAmplitude(name=Amp_name,
#                               data=zip(time[0:nrow, ], Q_Top_Data[:, i]),
#                               timeSpan=STEP)
#     myModel.SurfaceHeatFlux(name=Load_name,
#                             createStepName=Step_Name,
#                             region=Region_Array_Top_tuple[i],
#                             magnitude=1,
#                             distributionType=TOTAL_FORCE,
#                             amplitude=Amp_name)

    
# # Create the Bottom Amplitudes and Loads and assign them to each other
# for i in range(0, Num_Elem):
#     Amp_name = Amp_basename + str(i+1+Num_Elem)
#     Load_name = Load_basename + str(i+1+Num_Elem)
#     myModel.TabularAmplitude(name=Amp_name,
#                               data=zip(time[0:nrow, ], Q_Bottom_Data[:, i]),
#                               timeSpan=STEP)
#     myModel.SurfaceHeatFlux(name=Load_name,
#                             createStepName=Step_Name,
#                             region=Region_Array_Bottom_tuple[i],
#                             magnitude=1,
#                             distributionType=TOTAL_FORCE,
#                             amplitude=Amp_name)


region = myInstance.sets['Cut_Edge']
# # myModel.PinnedBC(name='BC-1',
# #                   createStepName='Initial',
# #                   region=region,
# #                   localCsys=None)
# # a = mdb.models[Model_Name].rootAssembly
# # region = a.instances['Reflector_Instance'].sets['Cut_Edge']
# mdb.models[Model_Name].EncastreBC(name='Encastre',
#                                           createStepName=Step_Name,
#                                           region=region,
#                                           localCsys=None)

mdb.models[Model_Name].EncastreBC(name='Encastre', 
    createStepName='Initial', region=region, localCsys=None)



# addition of predifined field for initial temperature
f1 = myInstance.faces[:]
e1 = myInstance.edges[:]
v1 = myInstance.vertices[:]
region = regionToolset.Region(vertices=v1, edges=e1, faces=f1)
myModel.Temperature(name='Initial_Temperature',
                    createStepName='Initial',
                    region=region,
                    distributionType=UNIFORM,
                    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
                    magnitudes=(350, ))

# # Creating the surface radiation for top and bottom
s1 = myInstance.faces[:]
# region_top=regionToolset.Region(side2Faces=s1)
# myModel.RadiationToAmbient(name='Top Radiation',
#                             createStepName='ReflectorLoad',
#                             surface=region_top,
#                             radiationType=AMBIENT,
#                             distributionType=UNIFORM,
#                             field='',
#                             emissivity=CFRP_emissivity,
#                             ambientTemperature=3.0,
#                             ambientTemperatureAmp='')


region_bottom=regionToolset.Region(side1Faces=s1)
myModel.RadiationToAmbient(name='Bottom Radiation',
                            createStepName='ReflectorLoad',
                            surface=region_bottom,
                            radiationType=AMBIENT, 
                            distributionType=UNIFORM,
                            field='',
                            emissivity=CFRP_emissivity,
                            ambientTemperature=3.0,
                            ambientTemperatureAmp='')


s1 = myInstance.faces[:]
surf1 = regionToolset.Region(side2Faces=s1)
surfaces =(surf1, )
myModel.CavityRadiationProp(name='IntProp-1',
                            property=(( CFRP_emissivity, ), ))
props =("IntProp-1", )
myModel.CavityRadiation(name='Cavity Radiation',
                        createStepName='ReflectorLoad',
                        surfaces=surfaces, 
                        surfaceEmissivities=props,
                        blocking=NO_BLOCKING,
                        ambientTemp=3.0)


#------------------------------------------------------------------------------
# Mesh the Reflector
#------------------------------------------------------------------------------
Create_Mesh(Model_Name, Part_Name, Viewport_Name, Mesh_Size)


#------------------------------------------------------------------------------
# Outputs
#------------------------------------------------------------------------------
# Field output
myModel.FieldOutputRequest(name='F-Output-1',
                            createStepName='ReflectorLoad',
                            variables=('NT',
                                      'TEMP',
                                      'HFL',
                                      'ENER'
                                      ),
                            timePoint='TimePoints-1')

# History output
myModel.HistoryOutputRequest(name='H-Output-1',
                              createStepName='ReflectorLoad',
                              variables=('FTEMP',
                                        'HFLA',
                                        'HTL',
                                        'HTLA',
                                        'WEIGHT'),
                              timePoint='TimePoints-1')



#------------------------------------------------------------------------------
# Create Job
#------------------------------------------------------------------------------

myJob = mdb.Job(name=Job_Name, model=Model_Name,
                description='Reflector Heat Transfer Analysis')

if getInput('Do you want to submit model? (y/n)') == 'y':
    # Wait for the job to complete.
    myJob.submit()
    myJob.waitForCompletion()
    
    o3 = session.openOdb(name='D:/1-Master Thesis/ABAQUS/Reflector/'+Job_Name+'.odb')
    session.viewports['Reflector Example'].setValues(displayedObject=o3)
    session.viewports['Reflector Example'].makeCurrent()
    session.viewports['Reflector Example'].odbDisplay.setPrimaryVariable(
        variableLabel='TEMP', outputPosition=INTEGRATION_POINT, )
    session.viewports['Reflector Example'].odbDisplay.display.setValues(
        plotState=CONTOURS_ON_DEF)
    
    if os.path.exists('D:/1-Master Thesis/ABAQUS/Reports/Reflector_'+Report_Name_Thermal+'.csv'):
      os.remove('D:/1-Master Thesis/ABAQUS/Reports/Reflector_'+Report_Name_Thermal+'.csv')
    else:
      print("The file does not exist")
    
    session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
    odb = session.odbs['D:/1-Master Thesis/ABAQUS/Reflector/'+Job_Name+'.odb']
    session.writeFieldReport(fileName='D:/1-Master Thesis/ABAQUS/Reports/Reflector_'+Report_Name_Thermal+'.csv',
                             append=OFF,
                             sortItem='Node Label',
                             odb=odb,
                             step=0, 
                             frame=50,
                             outputPosition=NODAL,
                             variable=(('TEMP', INTEGRATION_POINT), ),
                             stepFrame=ALL)

if getInput('Do you want to create a displacement model? (y/n)') == 'y':
    #------------------------------------------------------------------------------
    # Create the Displacement Model
    mdb.Model(name=Model_Name_Displacement, objectToCopy=mdb.models['Reflector_Thermal'])
    #: The model "Reflector_Displacement" has been created.
    p = mdb.models[Model_Name_Displacement].parts['Reflector_Part']
    
    del mdb.models[Model_Name_Displacement].steps['ReflectorLoad']
    mdb.models[Model_Name_Displacement].StaticStep(name='Step-1',
                              previous='Initial',
                              timePeriod=time[nrow],
                              maxNumInc=200,
                              description='Deflection')
    
    
    
    a = mdb.models[Model_Name_Displacement].rootAssembly
    region = a.instances['Reflector_Instance'].sets['Ref_Face']
    mdb.models[Model_Name_Displacement].Temperature(name='Predefined Field-2',
                                                createStepName='Step-1',
                                                region=region, 
                                                distributionType=FROM_FILE, 
                                                fileName='D:\\1-Master Thesis\\ABAQUS\\Reflector\\'+str(Job_Name)+'.odb', 
                                                beginStep=0,
                                                beginIncrement=0,
                                                endStep=1,
                                                endIncrement=97,
                                                interpolate=OFF, 
                                                absoluteExteriorTolerance=0.0,
                                                exteriorTolerance=0.05)
    
        
    elemType1 = mesh.ElemType(elemCode=S4R,
                              elemLibrary=STANDARD, 
                              secondOrderAccuracy=OFF)
    elemType2 = mesh.ElemType(elemCode=S3R,
                              elemLibrary=STANDARD)
    p = mdb.models[Model_Name_Displacement].parts['Reflector_Part']
    f = p.faces[:]
    pickedRegions =(f, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    
    mdb.models[Model_Name_Displacement].FieldOutputRequest(name='F-Output-1',
                                createStepName='Step-1',
                                variables=('CDISP',
                                            'CF',
                                            'CSTRESS',
                                            'LE',
                                            'PE',
                                            'PEEQ',
                                            'PEMAG',
                                            'RF',
                                            'S',
                                            'U'),
                                timePoint='TimePoints-1')
    
    mdb.models[Model_Name_Displacement].HistoryOutputRequest(name='H-Output-1',
                                  createStepName='Step-1',
                                  variables=('ALLAE',
                                              'ALLCD',
                                              'ALLDMD',
                                              'ALLEE',
                                              'ALLFD',
                                              'ALLIE',
                                              'ALLJD',
                                              'ALLKE',
                                              'ALLKL',
                                              'ALLPD',
                                              'ALLQB',
                                              'ALLSD',
                                              'ALLSE',
                                              'ALLVD',
                                              'ALLWK',
                                              'ETOTAL'),
                                  timePoint='TimePoints-1')
    
    a.regenerate()
    
    
    myJob_Disp = mdb.Job(name=Job_Name_Displacement, model=Model_Name_Displacement,
                    description='Reflector Deflection Analysis')

    if getInput('Do you want to submit model? (y/n)') == 'y':
        myJob_Disp.submit()
        myJob_Disp.waitForCompletion()
        
        o4 = session.openOdb(name='D:/1-Master Thesis/ABAQUS/Reflector/'+Job_Name_Displacement+'.odb')
        session.viewports['Reflector Example'].setValues(displayedObject=o4)
        session.viewports['Reflector Example'].makeCurrent()
        session.viewports['Reflector Example'].odbDisplay.setPrimaryVariable(
            variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
            'Magnitude'), )
        session.viewports['Reflector Example'].odbDisplay.display.setValues(
            plotState=CONTOURS_ON_DEF)
        
        if os.path.exists('D:/1-Master Thesis/ABAQUS/Reports/Reflector_'+Report_Name_Displacement+'.csv'):
          os.remove('D:/1-Master Thesis/ABAQUS/Reports/Reflector_'+Report_Name_Displacement+'.csv')
        else:
          print("The file does not exist")
        
        session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
        odb = session.odbs['D:/1-Master Thesis/ABAQUS/Reflector/'+Job_Name_Displacement+'.odb']
        session.writeFieldReport(fileName='D:/1-Master Thesis/ABAQUS/Reports/Reflector_'+Report_Name_Displacement+'.csv',
                                 append=OFF,
                                 sortItem='Node Label',
                                 odb=odb,
                                 step=0,
                                 frame=50,
                                 outputPosition=NODAL,
                                 variable=(('U', NODAL, ((INVARIANT, 'Magnitude'), )),),
                                 stepFrame=ALL)


session.viewports['Reflector Example'].odbDisplay.commonOptions.setValues(
    deformationScaling=UNIFORM, uniformScaleFactor=1)
session.viewports['Reflector Example'].viewportAnnotationOptions.setValues(
    legendNumberFormat=FIXED)
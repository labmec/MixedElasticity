import numpy as np
import json
import gmsh
import os

def PrintJson(JsonData: dict, FileName: str)->None:
    """ 
    PrintJson(JsonData, FileName)

    Reads the information from the 'JsonData' dictionary and writes them in a Json file
    named 'FileName'.

    Return: None

    Types: 
    - 'JsonData': dict
    - 'FileName': string 
    """

    json_object = json.dumps(JsonData, indent=4)

    with open(FileName+".json", "w") as outfile:
        outfile.write(json_object)

def PrintMeshInput(FileName: str)->None:
    """
    Creates the .geo and .msh files from a gmsh model. 

    Return: None 

    Types:
    - 'FileName': string
    """

    gmsh.write(FileName + ".geo_unrolled")
    gmsh.write(FileName + ".msh")

def MoveFiles(fileName: str, JsonNewPath: str, MeshNewPath:str)->None:
    """
    Move the .json, .msh, and .geo_unrolled files 
    to an given directory

    Return: None

    Types:
        - 'fileName': string
        - 'JsonNewPath': string
        - 'MeshNewPath': string
    """
    
    oldJson = fileName + ".json"
    oldmsh = fileName + ".msh"
    oldgeo = fileName + ".geo_unrolled"

    os.rename(oldJson, JsonNewPath+oldJson)
    os.rename(oldmsh, MeshNewPath+oldmsh)
    os.rename(oldgeo, MeshNewPath+oldgeo)

def CreatePoints(PointCoordinates: list, lc: float)->list:
    """
    Creates points from a list of 'PointCoordinates' with elements spaced with lc. 
    Important: currently it only generates points with sequential numbering! 

    Return: None

    Types:
    - PointCoordinates: list of list of float
    - step: float
    """
    pointsID = []
    for coord in PointCoordinates:
        x,y,z = coord
        point = gmsh.model.occ.addPoint(x,y,z,lc)
        pointsID.append(point)

    return pointsID

def CreateLines(LineIndexes: list)->list:
    """
    Creates lines from a list of 'LineIndexes'. In this list, it must be put 
    the index of each point that belongs to each line
    Important: currently it only generates lines with sequential numbering! 

    Types: 
    - 'LineIndexes': list of list of int
    """

    lineIds = []
    for index in LineIndexes:
        init, end = index
        lineIds.append(gmsh.model.occ.addLine(init, end))
    return lineIds

def CreateCurveLoops(CurveLoopIndexes: list)->None:
    """
    Creates curve loops from the 'CurveLoopIndex' list. In this list, it must be put the
    index of each line that belongs to each curve loop.

    Important: currently it only creates curve loops with sequential numbering

    Return: None

    Types: 
    - 'CurveLoopIndex': list of int
    """

    for i, index in enumerate(CurveLoopIndexes):
        gmsh.model.occ.addCurveLoop(index, i+1)

def CreatePlanes(PlaneIndexes: list)->None:
    """
    Creates plane surfaces from the 'PlaneIndex' list. In this list, it must be put the index 
    of each curve loop that belong to each plane surface. 

    Return: None

    Types:
    - 'PlaneIndexes': list of int
    """
    planes =[]
    for i, index in enumerate(PlaneIndexes):
        plane = gmsh.model.occ.addPlaneSurface(PlaneIndexes, i+1) 
        planes.append(plane)

    return planes

def CreatePhysicalGroup(GroupDimension: list, GroupIndex: list, GroupID: list, GroupName: list)->None:
    """
    Creates the Physical Group from the 'GroupDimension', 'GroupIndex', 'GroupID', 'GroupName' information.
    This function should use the same information to write the Json file! So that, the FEM simulation will 
    be properly set.

    Return: None

    Types: 
    - 'GroupDimension': list of int
    - 'GroupIndex': list of index
    - 'GroupID': list of int
    - 'GroupName': list of string
    """

    for dimension, index, id, BCname in zip(GroupDimension, GroupIndex, GroupID, GroupName):
        gmsh.model.addPhysicalGroup(dimension, index, tag=id, name=BCname)

def CreateCircles(Xcenter: float, Ycenter: float, Zcenter: float, Radius: float) -> int:
    """
    Creates surface circles from lists of 'Xcenter', 'Ycenter', 'Zcenter', and 'Radius'. 
    It first generates the circle as a drawong element in gmsh. Then, it transforms these circles 
    into curve loops to finally converts them into surfaces. Returns the surface circle TAG

    Return: surfaceCircle

    Types:
    - 'Xcenter': float
    - 'Ycenter': float
    - 'Zcenter': float
    - 'Radius': float
    - 'surfaceCircle': int
    """

    circle = gmsh.model.occ.addCircle(Xcenter,Ycenter,Zcenter,Radius)
    curveLoopCircle = gmsh.model.occ.addCurveLoop([circle])
    surfaceCircle = gmsh.model.occ.addPlaneSurface([curveLoopCircle])
    
    return surfaceCircle

def MakeHoles(domain: int, holesList: list, meshDim: int)->None:
    """
    Makes holes in a surface domain. Given a 'domain' tag, the 'holeList' tags, and the 'meshDim', it uses the 
    gmsh module cut to calculate the boolean difference the object domain and the object to be cut from it.

    Return: None

    Types:
    - 'domain': int
    - 'holesList': list of int
    - 'meshDim': int
    """
    
    holesTuple = [(meshDim, hole) for hole in holesList]
    
    gmsh.model.occ.cut([(meshDim,domain)], holesTuple)

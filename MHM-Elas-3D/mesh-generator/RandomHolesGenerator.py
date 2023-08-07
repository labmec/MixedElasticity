#%%
# ======================================
#           MODULES IMPORTED
# ======================================
import sys
import gmsh 
import random
import numpy as np
import TPZStokesInputCreator as Input 

#%%
# ======================================
#               FUNCTIONS 
# ======================================
def EuclideanDistance(xa: float, ya: float, xb: float, yb: float)->float:
    """
    Calculates the euclidean distance between the points (xa, ya) and (xb, yb)

    Return: Euclidean Distance
    """
    Xa = np.array([xa,ya])
    Xb = np.array([xb,yb])

    return np.linalg.norm(Xa-Xb)

def NoOverlappingCircles(domainX: tuple, domainY: tuple, radius: float, n_samples: int):
    """
    Constructs a list containing the (x,y,z) coordinates of 'n_samples' random circles, with no 
    overlapping between them. The range of the coordinates is within the ['domainX' x 'domainY'] values
    and the circles have radius = 'radius'

    Return: circleList
    """
    circleList = []

    xmin, xmax = domainX
    ymin, ymax = domainY

    while len(circleList) < n_samples:
        x = random.uniform(xmin, xmax)
        y = random.uniform(ymin, ymax)

        if not any((Xcenter, Ycenter, Zcenter) for Xcenter, Ycenter, Zcenter in circleList if EuclideanDistance(x, y,Xcenter, Ycenter)< 2.5*radius):
            circleList.append((x, y, 0))  

    return circleList

#%%
# ======================================
#               MAIN CODE
# ======================================
def main():
    # json file data
    json_Data = {
        "MeshName": "../Meshes/HolesEverywhere",
        "CreateMsh":False,
        "HdivType": 1,
        "VelpOrder": 1,
        "TracpOrder": 0,
        "Dim": 2,
        "Resolution": 0,
        "StaticCondensation": True,
        "isAxisymmetric": 0,
        "Domain": [
            {
                "name": "Stokes_Domain",
                "matID": 1,
                "viscosity": 1
            }
        ],
        "NormalBoundary": [
            {
                "name": "NoPenetration",
                "type": 0,
                "value": 0,
                "matID": 2
            },
            {
                "name": "PressureIn",
                "type": 2,
                "value": -10,
                "matID": 3
            },
            {
                "name": "PressureOut",
                "type": 2,
                "value": 0,
                "matID": 4
            }
        ],
        "TangentialBoundary": [
            {
                "name": "NoSlip",
                "type": 1,
                "value": 0,
                "matID": 5
            }
        ],
        "AxisymmetryDomain": [
        ],
        "AxisymmetryBoundary":[
        ],
        "LambdaID": 10,
        "InterfaceID": 20,
        "AxiLambdaID": 30,
        "AxiInterfaceID": 40,
        "FluxInterfaceID": 50
    }

    # input to create the mesh
    fileName = "RandomTest"

    random.seed(10)

    gmsh.initialize()

    n_samples = 20
    radius = .5
    recx = 10
    cmin  = .1*recx
    cmax = .9*recx
    
    lc = 1
    meshDim = 2
    points = [[0,0,0],[recx, 0, 0], [recx, recx, 0], [0,recx,0]]
    lineIndex = [[1,2],[2,3],[3,4],[4,1]]
    curveLoopIndex = [[1,2,3,4]]
    planeIndex = [1]
    physicalDim = [2,1,1,1,1]

    # creating the mesh using gmsh
    Input.CreatePoints(points, lc)
    Input.CreateLines(lineIndex)
    Input.CreateCurveLoops(curveLoopIndex)
    domain = Input.CreatePlanes(planeIndex)

    circleCoordinates = NoOverlappingCircles((cmin, cmax), (cmin, cmax), radius, n_samples)

    physicalName: list = [domain["name"] for domain in json_Data["Domain"]]+[condition["name"] for condition in json_Data["NormalBoundary"]] + \
    [condition["name"] for condition in json_Data["TangentialBoundary"]]

    physicalID: list = [domain["matID"] for domain in json_Data["Domain"]]+[condition["matID"] for condition in json_Data["NormalBoundary"]] + \
    [condition["matID"] for condition in json_Data["TangentialBoundary"]]

    circles = []
    for coordinate in circleCoordinates:
        Xcenter,Ycenter,Zcenter = coordinate 
        circles.append(Input.CreateCircles(Xcenter,Ycenter,Zcenter,radius))

    holesIndexes = [gmsh.model.occ.getCurveLoops(circle)[1][0][0] for circle in circles]

    Input.MakeHoles(domain[0], circles, meshDim)
    physicalIndex = [[1], [2,3]+holesIndexes, [4], [1], [1,2,3,4]+holesIndexes]

    gmsh.model.occ.synchronize()

    Input.CreatePhysicalGroup(physicalDim, physicalIndex, physicalID, physicalName)

    gmsh.model.mesh.generate(meshDim)

    # writing the .json, .msh, and .geo_unrolled files
    Input.PrintMeshInput(fileName)
    Input.PrintJson(json_Data, fileName)

    # moving the files to their repctive directories
    #Input.MoveFiles(fileName, "../DataInput/", "../Meshes/")

    # displaying the mesh created
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()

# calling the main function
if __name__ == "__main__":
    main()
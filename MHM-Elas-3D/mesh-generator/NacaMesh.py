
import sys

import gmsh

nacamat = -2
boundmat = -3
domainmat = 1

def Naca(par):
    yx = 5 * 0.12 * (0.2969 * par - 0.1260 * par*par - 0.3516 * par**4 + 0.2843 * par**6 - 0.1036 * par**8)
    return [par*par,yx,0]

gmsh.initialize()
meshdim = 2
PointCoordinates = [[-10,-10,0],[10,-10,0],[10,10,0],[-10,10,0]]
pointsID = []
lc = 1
for coord in PointCoordinates:
    x,y,z = coord
    point = gmsh.model.occ.addPoint(x,y,z,lc)
    pointsID.append(point)
quadlineID = [gmsh.model.occ.addLine(pointsID[t], pointsID[(t+1)%4]) for t in range(4)]
npt = 4
boundloop = gmsh.model.occ.add_curve_loop(quadlineID,boundmat)
gmsh.model.occ.synchronize()

NacaPointsId = []
lc = 0.1
for pt in range(10):
    par = pt/10.
    coord = Naca(par)
    x,y,z = coord
    pointid = gmsh.model.occ.addPoint(x,y,z,lc)
    NacaPointsId.append(pointid)
for pt in range(10):
    par = 1.-pt/10.
    coord = Naca(par)
    coord[1] = -coord[1]
    x,y,z = coord
    pointid = gmsh.model.occ.addPoint(x,y,z,lc)
    NacaPointsId.append(pointid)

gmsh.model.occ.synchronize()
  
npt = 20
NacaLineId = [gmsh.model.occ.addLine(NacaPointsId[t], NacaPointsId[(t+1)%npt]) for t in range(npt)]
gmsh.model.occ.synchronize()
nacaloop = gmsh.model.occ.add_curve_loop(NacaLineId,nacamat)
domain = gmsh.model.occ.addPlaneSurface([boundloop,nacaloop]) 

gmsh.model.occ.synchronize()
gmsh.model.addPhysicalGroup(1,quadlineID,boundmat,"bound")
gmsh.model.addPhysicalGroup(1,NacaLineId,nacamat,"naca")
gmsh.model.addPhysicalGroup(2,[domain],domainmat,"domain")

#gmsh.model.occ.cut([(meshdim,domain)], [(meshdim,nacaloop)])


gmsh.model.occ.synchronize()
# Generate mesh:
gmsh.model.mesh.generate(meshdim)
# displaying the mesh created
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()


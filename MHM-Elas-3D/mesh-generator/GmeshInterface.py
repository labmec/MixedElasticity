
import sys

from MHMesh1 import points,polylist
import TPZStokesInputCreator as Input
import gmsh

gmsh.initialize()
pointsIDs = Input.CreatePoints(points,1)
box = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
# box = gmsh.model.occ.addBox(-0.001, -0.001, -0.001, 1.001, 1.001, 1.001)
gmsh.model.occ.synchronize()


pairlist  = dict()
for poly in polylist:
    np = len(poly)
    for ip in range(np):
        p1 = poly[ip]
        ip2 = (ip+1)%np
        p2 = poly[ip2]
        p1id = pointsIDs[p1-1]
        p2id = pointsIDs[p2-1]
        if p1 < p2 : 
            p1p2 = tuple([p1id,p2id])
            if p1p2 not in pairlist:
                lineid = gmsh.model.occ.addLine(p1id, p2id)
                pairlist[p1p2] = lineid
gmsh.model.occ.synchronize()
faces = []
facebound = 4
faceinternal = 6
facename = dict()
facename[facebound] = "bound"
facename[faceinternal] = "internal"
for poly in polylist:
    np = len(poly)
    taglist = []

    delx=[0,0,0]
    for ip in range(np):
        p1 = poly[ip]
        ip2 = (ip+1)%np
        p2 = poly[ip2]
        delx=[max((points[p2-1][t]-points[p1-1][t]),delx[t]) for t in range(3)]
        p1id = pointsIDs[p1-1]
        p2id = pointsIDs[p2-1]
        if p1 < p2 : 
            p1p2 = tuple([p1id,p2id])
            lineid = pairlist[p1p2]
            taglist.append(lineid)
        else :
            p2p1 = tuple([p2id,p1id])
            lineid = -pairlist[p2p1]
            taglist.append(lineid)
    aligned=-1
    for i in range(3): 
        if delx[i] < 1.e-6: aligned = i
    if aligned != -1 :
        p1= poly[0]
        xco = points[p1-1][aligned]
        if (xco < 1.e-6 or xco > 1.-1.e-6) :
            facemat = facebound
        else:
            facemat = faceinternal
            faces.append([gmsh.model.occ.add_curve_loop(taglist),facemat])
    else:
        facemat = faceinternal
        faces.append([gmsh.model.occ.add_curve_loop(taglist),facemat])

    
facesurfacetag = []
facebytype = dict()
facebytype[facebound] = []
facebytype[faceinternal] = []          
for face in faces:
    facetag = gmsh.model.occ.add_plane_surface([face[0]])
    facesurfacetag.append(facetag)
    facebytype[face[1]].append(facetag)

gmsh.model.occ.synchronize()
# gmsh.model.addPhysicalGroup(2,facebytype[facebound],facebound,facename[facebound])
gmsh.model.addPhysicalGroup(2,facebytype[faceinternal],faceinternal,facename[faceinternal])



obj = [(3,box)]
tool = [(2,t) for t in facesurfacetag]
_,frag  = gmsh.model.occ.fragment(obj,tool,removeObject=True,removeTool=False)
gmsh.model.occ.synchronize()

boxlist = [tag for _,tag in frag[0] ]

boundaries = gmsh.model.getBoundary(frag[0],True,False,False)
regionsNF = [tag for _,tag in boundaries ]
tagnoflux = -1
namenoflux = "bound"
gmsh.model.add_physical_group(2, regionsNF, tagnoflux, namenoflux)


gmsh.model.occ.synchronize()
gmsh.model.addPhysicalGroup(3,boxlist,1,"volume")
# Generate mesh:
gmsh.model.mesh.generate()
# displaying the mesh created
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

FileName = "MHMesh1"
gmsh.write(FileName + ".geo_unrolled")
gmsh.write(FileName + ".msh")

gmsh.finalize()


# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 13:11:56 2015

@author: HannahWolfe
"""

vertices = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]]
sides = [[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],[4,7],[5,7],[6,7]]#verticesIndex
faces = [[1,2,6,8],[3,4,9,10],[0,2,4,7],[5,6,9,11],[0,1,3,5],[7,8,10,11]]#sidesIndex
iterations = 3

def checkAndAdd (item, array):
    try:
        array.index(item)
    except ValueError:
        array.append(item)
    return array.index(item)

def centroid (face, vertices, sides):
    x = []
    y = []
    z = []
    for side in face:
        for vertex in sides[side]:
            x.append(vertices[vertex][0])
            y.append(vertices[vertex][1])
            z.append(vertices[vertex][2])    
    return [sum(x) / (2.0*len(face)), sum(y) / (2.0*len(face)), sum(z) / (2.0*len(face))]

def midPoint (side, newVertices):
    v1= newVertices[side[0]]
    v2= newVertices[side[1]]
    return [(v1[0] + v2[0]) / 2.0, (v1[1] + v2[1]) / 2.0, (v1[2] + v2[2]) / 2.0]

def force (vertex,vertices):
    x = 0
    y = 0
    z = 0
    for v in vertices:
        x += v[0]
        y += v[1]
        z += v[2]
    return [vertex[0] + (.2* (vertex[0] - x / (1.0*len(vertices)))),
            vertex[1] + (.2* (vertex[1] - y / (1.0*len(vertices)))),
            vertex[2] + (.2* (vertex[2] - z / (1.0*len(vertices))))]
    
def iteration(vertices,sides,faces):
    newVertices = vertices[:]
    newSides = []
    newFaces = []
    
    for face in faces:
        print "face " + str(face)
        center = centroid(face, vertices, sides)
        centroidIndex = checkAndAdd(center, newVertices)   
        
        sideV = [sides[i] for i in face]
        points = list(set([item for sublist in sideV for item in sublist])) 
        for p in points:
            newFace = []
            for s in sideV:
                if p in s:
                    midpoint = midPoint(s, newVertices)
                    midPointIndex = checkAndAdd(midpoint, newVertices)
                    newFace.append(checkAndAdd(sorted([p,midPointIndex]), newSides))
                    newFace.append(checkAndAdd(sorted([centroidIndex,midPointIndex]), newSides))
            checkAndAdd(sorted(newFace), newFaces)
        
    for vertex in vertices:
        v = newVertices.index(vertex)
        newVertices[v] = force(vertex, vertices)
        
    print newVertices
    print newSides
    print newFaces
    return newVertices, newSides, newFaces
 
   
for i in range(0,iterations):
    vertices, sides, faces = iteration(vertices,sides,faces)
    
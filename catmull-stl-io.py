# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 13:11:56 2015

@author: HannahWolfe
"""

import numpy
from stl import mesh
import sys

vertices = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]]
sides = [[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],[4,7],[5,7],[6,7]]#verticesIndex
faces = [[1,2,6,8],[3,4,9,10],[0,2,4,7],[5,6,9,11],[0,1,3,5],[7,8,10,11]]#sidesIndex
iterations = 1
outputFile = "default_output.stl"

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
        
 #   print newVertices
 #   print newSides
 #   print newFaces
    return newVertices, newSides, newFaces

if len(sys.argv) > 1:
    inputFile = sys.argv[1]
    input_mesh = mesh.Mesh.from_file(inputFile)
    outputFile = "modified-"+inputFile
    vertices = []
    sides = []
    faces = []
    for v in range(0, len(input_mesh.v0)):
        v0 = checkAndAdd(list(input_mesh.v0[v]),vertices)
        v1 = checkAndAdd(list(input_mesh.v1[v]),vertices)
        v2 = checkAndAdd(list(input_mesh.v2[v]),vertices)
        s0 = checkAndAdd(sorted([v0,v1]),sides)
        s1 = checkAndAdd(sorted([v2,v1]),sides)
        s2 = checkAndAdd(sorted([v0,v2]),sides)
        checkAndAdd(sorted([s0,s1,s2]),faces)

print vertices
print sides
print faces
   
for i in range(0,iterations):
    print "iteration: " + str(i)
    vertices, sides, faces = iteration(vertices,sides,faces)
  
data = numpy.zeros(2*len(faces), dtype=mesh.Mesh.dtype)
final_mesh = mesh.Mesh(data, remove_empty_areas=False)

for f in range(0,len(faces)):
    sideV = [sides[i] for i in faces[f]]
    points = list(set([item for sublist in sideV for item in sublist])) 
    p = [points[0]]
    for s in sideV:
        if p[0] in s:
           p.append([item for item in s if item not in p][0])
    final_mesh.v0[2*f] = vertices[p[0]]
    final_mesh.v1[2*f] = vertices[p[1]]
    final_mesh.v2[2*f] = vertices[p[2]]
    p1 = [item for item in points if item not in p]
    for s in sideV:
	if p1[0] in s:
            p1.append([item for item in s if item not in p1][0])
    final_mesh.v0[2*f+1] = vertices[p1[0]] 
    final_mesh.v1[2*f+1] = vertices[p1[1]]
    final_mesh.v2[2*f+1] = vertices[p1[2]]
print final_mesh.v0
print final_mesh.v1
print final_mesh.v2
print final_mesh.points[0]
final_mesh.save(outputFile)

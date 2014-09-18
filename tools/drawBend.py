from pgmagick import *
import math, sys, re
from os import path

argv = Argv()
#argv = sys.argv[1:]

DefineScalarExpression("x_comp", "array_decompose(<E-Field>,0)")
DefineScalarExpression("y_comp", "array_decompose(<E-Field>,1)")

def drawStep(stepNum, s, g1, g2, g3):
    baseName = s.fileName
    serieNumber = "%04d" % stepNum
    fname = "c_" + serieNumber + ".pvtr"
    s.fileName = baseName + serieNumber
    print("processing %s" % fname)

    OpenDatabase(fname)
    SetSaveWindowAttributes(s)

    # print x component
    AddPlot("Pseudocolor", "x_comp")
    DrawPlots()

    SaveWindow()
    DeleteAllPlots()

    im = Image(s.fileName + ".png")
    im.fillColor(Color("black"))
    im.draw(g1)
    im.draw(g2)
    im.draw(g3)

    imfnameOut = baseName + "_xcomp" + serieNumber + ".png"
    im.write(imfnameOut)

    # print y component
    AddPlot("Pseudocolor", "y_comp")
    DrawPlots()

    SetSaveWindowAttributes(s)
    SaveWindow()
    DeleteAllPlots()

    im = Image(s.fileName + ".png")
    im.fillColor(Color("black"))
    im.draw(g1)
    im.draw(g2)
    im.draw(g3)

    imfnameOut = baseName + "_ycomp" + serieNumber + ".png"
    im.write(imfnameOut)

    # print magnitude
    AddPlot("Pseudocolor", "E-Field_magnitude")
    DrawPlots()

    SetSaveWindowAttributes(s)
    SaveWindow()
    DeleteAllPlots()

    im = Image(s.fileName + ".png")
    im.fillColor(Color("black"))
    im.draw(g1)
    im.draw(g2)
    im.draw(g3)

    imfnameOut = baseName + "_magnitude" + serieNumber + ".png"
    im.write(imfnameOut)

    CloseDatabase(fname)
    s.fileName = baseName

if len(argv) == 0:
    print("no arguments given")
    sys.exit()

l, lb, w, ekin, bz, p = -0.1, -0.1, -0.1, -0.1, 0.0, 0.0
infile = argv[0]
fh = open(infile,'r')
while fh:
    line = fh.readline()
    if not line:
        break
    a = re.search('--length (\d*\.\d*) ', line)
    if a:
        l = float(a.group(1))
    a = re.search('--lengthbefore (\d*\.\d*) ', line)
    if a:
        lb = float(a.group(1))
    a = re.search('--width (\d*\.\d*) ', line)
    if a:
        w = float(a.group(1))
    a = re.search('--Bz (\d*\.\d*) ', line)
    if a:
        bz = float(a.group(1))
    a = re.search('--Ekin (\d*\.?\d*) ', line)
    if a:
        ekin = float(a.group(1))
    a = re.search('--phi (\d*\.?\d*) ', line)
    if a:
        p = float(a.group(1)) / 180.0 * math.pi

if l < 0.0 or lb < 0.0 or w < 0.0 or ekin < 0.0:
    print("could not find length or width")
    sys.exit()

if math.fabs(bz) > 1e-8:
    R = math.sqrt(ekin**2 + 2.0 * ekin * 0.511) * 1e6 / (299792458 * math.fabs(bz))
else:
    R = 99999999

outfileBase = "overview"
ii = 0
while path.exists(outfileBase + "%04d.png" % (ii)):
    os.remove(outfileBase + "%04d.png" % (ii))
    ii += 1

fname = 'c_0000.pvtr'
if not path.exists(fname):
    print("can't find files c_*.pvtr")
    sys.exit()

OpenDatabase("c_*.pvtr database")
# OpenDatabase(fname)
# AddPlot("Pseudocolor", "x_comp") #E-Field_magnitude")
# DrawPlots()
# Query("MinMax")
# lmin,lmax = GetQueryOutputValue()
# DeleteAllPlots()
nSteps = TimeSliderGetNStates()
CloseDatabase(fname)

domainLength = lb + (R + w / 2) * math.sin(p) + l * math.cos(p)
domainWidth = (R - w / 2) * (1 - math.cos(p)) + w + l * math.sin(p)

s = SaveWindowAttributes()
s.format = s.PNG
s.fileName = outfileBase
s.family = 0
s.resConstraint = s.NoConstraint
s.width = 4096
s.height = int(math.floor(s.width / domainLength * domainWidth))
s.screenCapture = 0
SetSaveWindowAttributes(s)

dx = domainLength / (s.width / 4 * 3 - 2)

x0 = s.width / 5 + 1 #820
x1 = x0 + int(math.floor(lb/dx + 0.5))
x2 = x1 + int(math.floor(math.sin(p) * (R - w / 2) / dx + 0.5))
x3 = x2 + int(math.floor(math.cos(p) * l / dx + 0.5))
x4 = x3 + int(math.floor(math.sin(p) * w / dx + 0.5))
x5 = x4 - int(math.floor(math.cos(p) * l / dx + 0.5))

y0 = s.height * 85 / 100 - w / dx + 1 #523
y1 = y0 - int(math.floor((1 - math.cos(p)) * (R - w / 2) / dx + 0.5))
y2 = y1 - int(math.floor(math.sin(p) * l / dx + 0.5))
y3 = y2 + int(math.floor(math.cos(p) * w / dx + 0.5))
y4 = y3 + int(math.floor(math.sin(p) * l / dx + 0.5))
y5 = y4 + int(math.floor((1 - math.cos(p)) * (R + w / 2) / dx + 0.5))

rx0 = x1 - int(math.floor((R - w / 2) / dx + 0.5))
rx1 = x1 + int(math.floor((R - w / 2) / dx + 0.5))
rx2 = x1 - int(math.floor((R + w / 2) / dx + 0.5))
rx3 = x1 + int(math.floor((R + w / 2) / dx + 0.5))

ry0 = y0
ry1 = ry0 - int(math.floor(2 * (R - w / 2) / dx + 0.5))
ry2 = y5
ry3 = ry2 - int(math.floor(2 * (R + w / 2) / dx + 0.5))

#g0 = DrawableArc(rx0,ry0-4,rx1,ry1-4,85,90)

l1 = CoordinateList()
l1.push_back(Coordinate(x0,y0-4))
l1.push_back(Coordinate(x1,y0-4))
for i in range(20):
    t = (5.0 * (i + 1) / 22.0) * math.pi / 180.0
    deltaX = int(math.floor((R - w / 2) * math.sin(t) / dx + 0.5))
    deltaY = int(math.floor((R - w / 2) * (1 - math.cos(t)) / dx + 0.5))
    l1.push_back(Coordinate(x1 + deltaX, y0 - deltaY - 4))
l1.push_back(Coordinate(x2,y1-4))
l1.push_back(Coordinate(x3,y2-4))
l1.push_back(Coordinate(x0,y2-4))
g1 = DrawablePolygon(l1)

l2 = CoordinateList()
l2.push_back(Coordinate(x3,y2-4))
l2.push_back(Coordinate(x4,y3-4))
l2.push_back(Coordinate(x4,y2-3))
g2 = DrawablePolygon(l2)

l3 = CoordinateList()
l3.push_back(Coordinate(x4,y3-4))
l3.push_back(Coordinate(x5,y4-2))

for i in range(20):
    t = (5.0 - 5.0 * (i + 1) / 22.0) * math.pi / 180.0
    deltaX = int(math.floor((R + w / 2) * math.sin(t) / dx + 0.5))
    deltaY = int(math.floor((R + w / 2) * (1 - math.cos(t)) / dx + 0.5))
    l3.push_back(Coordinate(x1 + deltaX, y5 - deltaY - 2))

l3.push_back(Coordinate(x1,y5-2))
l3.push_back(Coordinate(x4,y5-2))
g3 = DrawablePolygon(l3)

for ii in range(nSteps):
    drawStep(ii, s, g1, g2, g3)

sys.exit()

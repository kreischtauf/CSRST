#!/usr/bin/python2

from pgmagick import *
import sys, os, math

baseName = "overview_xcomp"
for ii in range(52,55):
    serieNumber = "%04d" % ii
    fname = baseName + serieNumber + ".png"
    imTotal = Image("EndFrame.png")
    imTotal.write(fname)

ib = 55
for ii in range(55,62):
    serieNumber = "%04d" % ii
    fname = baseName + serieNumber + ".png"
    imTotal = Image("EndFrame.png")
    imBlack = Image("BlackFrame%d.png" % (ii - ib))
    imTotal.composite(imBlack, 0, 0, CompositeOperator.AtopCompositeOp)

    imTotal.write(fname)

ib = 62
for ii in range(62, 82):
    serieNumber = "%04d" % ii
    fname = baseName + serieNumber + ".png"
    imTotal = Image("EndCredits.png")
    imTotal.crop(Geometry(1505, 1205, 0, int(math.floor((ii - ib) * 115.32 + 0.5))))
    imTotal.write(fname)

for ii in range(82,85):
    serieNumber = "%04d" % ii
    fname = baseName + serieNumber + ".png"
    imTotal = Image("EndCredits.png")
    imTotal.crop(Geometry(1505, 1205, 0, 2191))
    imTotal.write(fname)

ib = 85
for ii in range(85,92):
    serieNumber = "%04d" % ii
    fname = baseName + serieNumber + ".png"
    imTotal = Image("EndCredits.png")
    imTotal.crop(Geometry(1505, 1205, 0, 2191))
    imBlack = Image("BlackFrame%d.png" % (ii - ib))
    imTotal.composite(imBlack, 0, 0, CompositeOperator.AtopCompositeOp)

    imTotal.write(fname)

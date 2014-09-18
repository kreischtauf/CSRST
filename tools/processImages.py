#!/usr/bin/python2

from pgmagick import *
import sys, os

baseName = "overview_ycomp"
for ii in range(52):
    serieNumber = "%04d" % ii
    fname = baseName + serieNumber + ".png"
    imTotal = Image(fname)

    imLegend = Image(imTotal)
    imLegend.crop(Geometry(444, 563, 297, 19))

    imData = Image(imTotal)
    imData.fillColor(Color("white"))
    imData.draw(DrawableRectangle(4583, 1312, 5860, 1492))
    imData.draw(DrawableRectangle(5111, 1309, 5120, 1312))
    imData.draw(DrawableRectangle(4583, 1288, 5099, 1312))
    imData.draw(DrawableRectangle(3469, 1320, 3601, 1347))
    imData.fillColor(Color("black"))
    imData.draw(DrawableRectangle(5837, 149, 6144, 1278))
    imData.crop(Geometry(754, 1165, 1230 + ii * 84, 148))

    imAxis = Image(imTotal)
    imAxis.crop(Geometry(1505, 1205, 576, 131))
    imAxis.fillColor(Color("white"))
    imAxis.draw(DrawableRectangle(0, 0, 165, 451))
    imAxis.draw(DrawableRectangle(625, 1161, 661, 1183))
    imAxis.draw(DrawableRectangle(654, 0, 1505, 1205))
    imAxis.fillColor(Color("black"))
    imAxis.draw(DrawableRectangle(654, 18, 1407, 1147))
    imAxis.composite(imLegend, 67, 0)
    imAxis.composite(imData, 654, 17)
    imAxis.write("processed/" + fname)

sys.exit()

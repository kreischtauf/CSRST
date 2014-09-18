#!/usr/bin/python2

import sys

component = "xcomp" # xcomp ycomp, magnitude
componentString = "X Component"
argv = Argv()
numFiles = len(argv)
print argv

DefineScalarExpression("E-Field_xcomp", "array_decompose(<E-Field>,0)")
DefineScalarExpression("E-Field_ycomp", "array_decompose(<E-Field>,1)")
# DefineScalarExpression("E-Field_xcomp", "array_decompose(<Evel-Field>,0)+array_decompose(<Edrift-Field>,0)")
# DefineScalarExpression("E-Field_ycomp", "array_decompose(<Evel-Field>,1)+array_decompose(<Edrift-Field>,1)")
# DefineScalarExpression("E-Field_xcomp", "array_decompose(<Eacc-Field>,0)")
# DefineScalarExpression("E-Field_ycomp", "array_decompose(<Eacc-Field>,1)")
outfname = "frame_" + component

if numFiles < 1:
    print("no arguments given")
    sys.exit()

ii = 0
while os.path.exists(outfname + "%04d.png" % (ii)):
    os.remove(outfname + "%04d.png" % (ii))
    ii += 1

ii = 1
for v in argv:
    OpenDatabase(v + " database")
    nsteps = TimeSliderGetNStates()
    AddPlot("Pseudocolor", "E-Field_" + component)
    # p = PseudocolorAttributes()
    # p.colorTableName = "hot_and_cold"
    # SetPlotOptions(p)
    #SetTimeSliderState(10)
    #if ii != numFiles:
    #    p = PseudocolorAttributes()
    #    p.legendFlag = 0
    #    SetPlotOptions(p)
    #ii += 1

#DrawPlots()
#Query("MinMax")
#lmin,lmax = GetQueryOutputValue()
def animate():
   print "animating"
   fname = string.replace(argv[0], '*', "0000", 1)
   OpenDatabase(argv[0] + " database")
   s = SaveWindowAttributes()
   s.format = s.PNG
   s.resConstraint = s.NoConstraint
   s.fileName = outfname
   s.width,s.height = 2126,2126 #3072,3072
   s.screenCapture = 0
   SetSaveWindowAttributes(s)
   nsteps = TimeSliderGetNStates()
   componentTitle = CreateAnnotationObject("Text2D")
   componentTitle.text = componentString
   componentTitle.position = (0.5, 0.94)
   componentTitle.fontBold = 0
   componentTitle.height = 0.025
   DeleteAllPlots()

#   fh = open('MinMaxValues.dat','w')
   fh = open('MinMaxValues.dat','r')
   for i in range(nsteps):
       serieNumber = "%04d" % (i)
       fname = string.replace(argv[0], '*', serieNumber, 1)
       # fname = "mm_%04d.pvtr" % (i)
       OpenDatabase(fname)
       AddPlot("Pseudocolor", "E-Field_" + component) #"E-Field_magnitude")
       AddOperator("Transform")
       DrawPlots()
       #
       Query("SpatialExtents")
       cvals = GetQueryOutputValue()
       # Query("MinMax")
       # lmin,lmax = GetQueryOutputValue()
       # fh.write("%12.4f %12.4f\n" % (lmin, lmax))

       lmin,lmax = [float(x) for x in fh.readline().split()]
       p = PseudocolorAttributes()
       p.min,p.minFlag = lmin,1
       p.max,p.maxFlag = lmax,1
       # p.colorTableName = "hot_and_cold"
       SetPlotOptions(p)

       #
       v = GetView2D()
       # cX = 0.5 * (cvals[0] + cvals[1])
       # cY = 0.5 * (cvals[2] + cvals[3])
       # v.windowCoords = (-4.1056e-3, 4.1056e-3, -5.51181e-3, 5.51181e-3)
       # v.viewportCoords = (0.38, 0.99, 0.1, 0.93)
       v.windowCoords = (-4.1056e-3, 4.1056e-3, -5.11811e-3, 5.11811e-3)
       v.viewportCoords = (0.2955, 0.981, 0.08, 0.93)
       SetView2D(v)
       #
       if (cvals[0] > 0.0):
           tr_atts = TransformAttributes()
           tr_atts.doTranslate = 1
           tr_atts.translateX = -0.5 * (cvals[0] + cvals[1])
           tr_atts.translateY = -0.5 * (cvals[2] + cvals[3])
           SetOperatorOptions(tr_atts)
       #
       a = AnnotationAttributes()
       a.userInfoFlag = 0
       a.timeInfoFlag = 0
       a.axes2D.xAxis.title.font.scale = 1.5
       a.axes2D.xAxis.label.font.scale = 1.5
       a.axes2D.yAxis.title.font.scale = 1.5
       a.axes2D.yAxis.label.font.scale = 1.5
       a.axes2D.autoSetTicks = 0
       a.axes2D.xAxis.tickMarks.majorMinimum = -0.006
       a.axes2D.xAxis.tickMarks.majorMaximum = 0.006
       a.axes2D.xAxis.tickMarks.majorSpacing = 0.002
       a.axes2D.xAxis.tickMarks.minorSpacing = 0.0002
       a.axes2D.yAxis.tickMarks.majorMinimum = -0.006
       a.axes2D.yAxis.tickMarks.majorMaximum = 0.006
       a.axes2D.yAxis.tickMarks.majorSpacing = 0.002
       a.axes2D.yAxis.tickMarks.minorSpacing = 0.0002
       SetAnnotationAttributes(a)
       #
       refNames = GetAnnotationObjectNames()
       ref = GetAnnotationObject(refNames[-1])
       ref.fontHeight = 0.03
       ref.yScale = 2
       ref.drawTitle = 0
       #
       DrawPlots()
#       SetTimeSliderState(i)
#       p = PseudocolorAttributes()
#       p.legendFlag = 0
#       p.min,p.minFlag = lmin,1
#       p.max,p.maxFlag = lmax,1
#       SetPlotOptions(p)
#       for j in range(numFiles - 1):
#           SetActivePlots(j)
#       SetActivePlots(tuple(range(numFiles)))
#       RedrawWindow()
       SaveWindow()
       DeleteAllPlots()
       CloseDatabase(fname)
   fh.close()
   print "done"

animate()

sys.exit()

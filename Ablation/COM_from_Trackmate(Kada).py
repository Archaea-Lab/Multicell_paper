from ij import IJ
from ij import WindowManager 
from ij.plugin.frame import RoiManager 
from ij.gui import Roi 
from ij.process import ImageStatistics as IS 
from ij.measure import ResultsTable 

imp = WindowManager.getCurrentImage() 
# Add points to ROI manager 

#print(imp.getStack()) 
rm = RoiManager.getInstance() 
if not rm: 
    rm = RoiManager() 
options = IS.MEAN | IS.AREA | IS.CENTROID | IS.CENTER_OF_MASS 
count=rm.getCount() 
rois=rm.getRoisAsArray() 
rt=ResultsTable() 

for roi in rois: 

	imp.setRoi(roi) 
	stats=imp.getStatistics(options) 
	rt.incrementCounter() 
	rt.addValue("Label",roi.getName()) 
	name=roi.getName()
	Slice=rm.getSliceNumber(name)
	IJ.setSlice(Slice)
	rt.addValue("Slice",rm.getSliceNumber(name)) 
	rt.addValue("xCOM",stats.xCenterOfMass) 
	rt.addValue("yCOM",stats.yCenterOfMass) 
	rt.addValue("xCentroid",stats.xCentroid) 
	rt.addValue("yCentroid",stats.yCentroid) 
	rt.addValue("Area",stats.area) 
	rt.sort("Slice") 


length=rt.size() 

print(r'Total number of frames in the roi are ',length) 


rt.show("Single_Track_X_Y_COM_Centroid") 
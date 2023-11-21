dir_in=getDirectory("Choose_input_directory"); 
dir_out=getDirectory("Choose_output_directory"); 
list= getFileList (dir_in); 
path = dir_in + list[0];
run("Image Sequence...", "open=&path file=Index sort");
run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=150 frames=20 display=Color");
makeRectangle(388, 420, 1200, 1200);
run("Crop");
run("Bleach Correction", "correction=[Simple Ratio]");

for (i=0; i<list.length; i++)
{
	Stack.setFrame(i+1) 
	run("Reduce Dimensionality...", "slices keep");
	saveAs("Tiff", dir_out+list[i]);
	close();
	}
run("Close All");	




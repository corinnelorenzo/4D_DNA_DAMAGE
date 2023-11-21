path_tiff=getDirectory("Choose_input_directory");
path_M=getDirectory("Choose_output_directory")
list=getFileList (path_tiff); 
path=dir_in + list[0];
setBatchMode(true);
// OPEN STACK SERIES// OPEN STACK SERIES
for (i=0; i<list.length; i++) {
     showProgress(i+1, list.length);
     filename = dir_in + list[i];
     if (endsWith(filename, "tif")) {
     open(filename);for (i=0; i<list.length; i++) {
     showProgress(i+1, list.length);
     filename = dir_in + list[i];
     if (endsWith(filename, "tif")) {
     open(filename);
 // 3D WATERSHED & 3D MANAGR
run("3D Watershed", "seeds_threshold=1 image_threshold=1 image=170906BleoR1 seeds=Automatic radius=2");
run("3D Manager");
Ext.Manager3D_AddImage();
// do some measurements, save measurements and close window
Ext.Manager3D_Measure();
Ext.Manager3D_SaveResult("M","/Users/lorenzo/Desktop/Output/");
Ext.Manager3D_CloseResult("M");


dataFileName = "fileName.csv";

path = "/path/to/image/folder/";

imgLst = getFileList(path);


for (i = 0; i < imgLst.length; i++)
{
	
	open( path+imgLst[i] );
	run("Flip Horizontally");
	makeRectangle(0, 30, 1280, 660);
	run("Crop");
	setAutoThreshold("Default");
	setThreshold(0, 90);
	setOption("BlackBackground", false);
	run("Convert to Mask");

	makeRectangle(0, 164, 184, 126);
	changeValues(0, 255, 0);

	makeRectangle(0, 225, 28, 152);
	run("Set Measurements...", "area centroid display redirect=None decimal=3");
	run("Analyze Particles...", "size=10-40 display");

	makeRectangle(28, 0, 1280, 660);
	run("Analyze Particles...", "size=10-40 display");

	close();



}


selectWindow("Results");
saveAs("Text", "/results/saving/folder/"+dataFileName);

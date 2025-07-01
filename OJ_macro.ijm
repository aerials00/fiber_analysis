close("*");
List.clear();

dir = "maindirectory"; // Directory that has everything (reflection + GFP)
coldir  =  dir + "collagen_test/"; // Directory containing the test collagen images (reflection only)

ojdir   = dir + "Orientation_Masks/"; // Directory where orientation masks are stored
File.makeDirectory(ojdir);

num_files = 10; //Insert number of  files that are present in the collagen_test folder 

files = getFileList(coldir); // Get the files inside the collagen directory

// Start analysis 
for(f=0; f<num_files; f++){
	
	id = coldir+files[f]; // the path to file
	open(id); // open file
	fullname = getTitle(); // Read the title of open image in ImageJ
	name = substring(fullname, 0,lastIndexOf(fullname, ".")); // Print the name of the file
	print(name);
	
	// Run the Orientation analysis
	run("OrientationJ Vector Field", "tensor=1 gradient=4 orientation=on radian=off"); // tensor = 1 is the sigma and gradient 0 is cubic spline, 4 is gaussian blur radian = 0ff--> degrees
	IJ.renameResults("Results");
		
	col = getTitle();
	print(col);

	selectImage("OJ-Orientation-1");
	OJ_name = "Orientation_Mask_" + name;
	saveAs("tiff", ojdir + OJ_name);

		
	run("Clear Results");
	close("*");
}

run("Close All");
showMessage("Done");
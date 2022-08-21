#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix
// See also Process_Folder.py for a version of this code
// in the Python scripting language.
processFolder(input);
// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
	print("Finished!");
}
function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	setBatchMode(true); // true prevents image windows from opening while the script is running
	open(input + File.separator + file); //open file

	print("Processing: " + file);
	//run("Bandpass Filter...", "filter_large=4 filter_small=2 suppress=None tolerance=5 autoscale saturate process");
    //run("Enhance Contrast", "saturated=0.35");
    rpx = "5";
	run("Subtract Background...", "rolling=" + rpx +" stack");

	print("Saving");
	//saveAs("tiff", output + File.separator + "bandpass4-2px_" + file);
	saveAs("tiff", output + File.separator + "rollingball" + rpx + "px_" + file);
}

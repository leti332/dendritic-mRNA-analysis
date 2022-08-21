// Save Dendritic Outlines
// WORKS - 18MAR2021
// Editor: Leti Nunez
//////////////////////////////////////////////////////////////////

// 1. Load Average Z-projection (see Kymograph Analysis MATLAB script) and raw movie
// 2. Use Segmented Line on Avg Z-project to segment dendrites
// Note: Start from soma (if no soma in movie, base on branching direction. Bifurcating
// branching usually are in the anterograde direction. In other words, the "Y" should
// be to the right when segmenting.
// 3. Run this script to save dendritic outline


// Specify the dendrite number
	dNum = "-" + "2";

// Apply dendritic outline to duplicated image (for saving purposes)
	orig = getImageID(); //get imageid of selected image
	//print(title);//sanity check
	run("Duplicate...", " "); // make duplicate image it can be saved with outline

	new = getImageID(); // save imageid of duplicate image

	selectImage(orig); //select original image to get the outline selection
	selectImage(new); //select duplicate image
	run("Restore Selection"); // apply outline to duplicate image

// Change Title to reflect dendrite number
	selectImage(orig); //select original image for title
	title = getTitle();
	idx = indexOf(title, "."); // get index of period
	//print(idx);//sanity check

	Title = substring(title, 4, idx); //get title before the period
	print(Title);//sanity check


	Title = "Outline_" + Title + dNum + ".tif"// save outline
	//print(Title);//sanity check

// Get Output Directory for Image
//out = getDirectory("Output Path"); // GUI selects directory
	out = "C:/Users/letinunez/Dropbox (EinsteinMed)/Projects/RNA-Protein_Dynamic_Modeling/ACTBmRNAtracking_ZBP1KO/PAPER_FINAL/Data/2_Dendrite Outlines/WT/exp3_20150717"
	//print(out);//sanity check

// Save Image
	selectImage(new);
	saveAs("Tiff", out + File.separator + Title);

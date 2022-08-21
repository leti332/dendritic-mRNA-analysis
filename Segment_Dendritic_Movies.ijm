#@ File(label="Raw Movie directory", style = "directory") movies
#@ File(label="Dendritic Outline directory", style = "directory") outlines
#@ File(label="Segmented Dendrite Output directory", style = "directory") output
#@ String(label = "File Suffix", value = ".tif") suffix


processFolder(movies, outlines);

function processFolder(movies, outlines) {
	print("Processing...");
	mlist = getFileList(movies);

	olist = getFileList(outlines);

	for (i = 0; i < mlist.length; i++) {
		if(File.isDirectory(movies + mlist[i]))
			processFolder("" + movies + mlist[i]);
			//print("Movie: " + mlist[i]);
			idx = indexOf(mlist[i], suffix); //get index of suffix
			Movie = substring(mlist[i],0, idx);// get name without suffix
			//print("Movie: " + Movie); //print name without suffix
		if(endsWith(mlist[i], suffix)){
			//Loop through outline folder
			for (k = 0; k < olist.length; k++){
				if(File.isDirectory(outlines + olist[k]))
					processFolder("" + outlines + olist[k]);
					//print("Outline: " + olist[k]);

				// Conditional if the Movie name is a substring of the outline name
				//Execute following code
				if(indexOf(olist[k], Movie) >= 0){
					strDend(outlines, movies, olist[k], mlist[i], output);
				}else{
					print(olist[k] + " & " + mlist[i] + " are not Matches");
				}
			}
		}
	}
	print("Finished!");

}

function strDend(outlines, movies, outline, movie, output){
	setBatchMode(true); // true prevents image windows from opening while the script is running

	print("Matches: " + outline + " & " + movie);
	open(outlines + File.separator + outline);
	print("Processing: " + outlines + File.separator + outline);
	oId = getImageID();
	otitle = getTitle();


	open(movies + File.separator + movie);
	mId = getImageID();
	mtitle = getTitle();
	selectImage(oId);
	selectImage(mId);
	run("Restore Selection");
	run("Straighten...", "title line=50 process");

	// Remove string "Outline_" from title name
	rm = "Outline_";
	index = lengthOf(rm);

	ix = lastIndexOf(otitle, "Outline_");
	index = index + ix;

	file = substring(otitle,index);

	// Save File
	print("Saving Dendritic: " + file);
	saveAs("tiff", output + File.separator + file);
}

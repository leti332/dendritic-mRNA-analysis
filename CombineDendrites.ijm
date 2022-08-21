#@ File(label="Dendrite directory", style = "directory") movies
#@ File(label="Combined Dendrite Output directory", style = "directory") output
#@ String(label = "File Suffix", value = ".tif") suffix


processFolder(movies);

function processFolder(movies) {
	print("Processing...");
	mlist = getFileList(movies);

	for (i = 0; i < mlist.length; i++) {

		file = "KO_Comb_20140717_KO23.tif";

		if(endsWith(mlist[i], suffix)){
			if(i == 0){
				print(i);
				comb1stDend(movies, mlist[i], mlist[1], output, file);
			}
			if(i > 1){
				print("Progress: " + i + "/" + mlist.length);
				combDend(movies, file, mlist[i], output);
			}

		}
	}
	print("Finished!");
}

function comb1stDend(movies, movie1, movie2, output, file){
	setBatchMode(true); // true prevents image windows from opening while the script is running

	//Movie1
	open(movies + File.separator + movie1); //open file
	print("Movie1: " + movie1);
	idx = indexOf(movie1, suffix); //get index of suffix
	Movie1 = substring(movie1,0, idx);// get name without suffix
	//M1title = getTitle();
	//print("Movie1: " + M1title);

	//Movie2
	open(movies + File.separator + movie2); //open file
	print("Movie2: " + movie2);
	idx = indexOf(movie2, suffix); //get index of suffix
	Movie2 = substring(movie2,0, idx);// get name without suffix

	//Combine Dendrites
	run("Combine...", "stack1=" + movie1 + " stack2=" + movie2 + " combine");

	print(output);
	saveAs("tiff", output + File.separator + file);

	close();
	print("1st Combination Completed");
}

function combDend(movies, movie1, movie2, output){
	setBatchMode(true); // true prevents image windows from opening while the script is running

	//Movie1 - Combined Movie
	open(output + File.separator + movie1); //open file
	print("Movie1: " + movie1);
	idx = indexOf(movie1, suffix); //get index of suffix
	Movie1 = substring(movie1,0, idx);// get name without suffix
	//M1title = getTitle();
	//print("Movie1: " + M1title);

	//Movie2
	open(movies + File.separator + movie2); //open file
	print("Movie2: " + movie2);
	idx = indexOf(movie2, suffix); //get index of suffix
	Movie2 = substring(movie2,0, idx);// get name without suffix

	//Combine Dendrites
	print("Combining: " + movie1 + " & " + movie2);
	run("Combine...", "stack1=[" + movie1 +"] stack2="+ movie2 +" combine");

	saveAs("tiff", output + File.separator + movie1);
	close();
	print("Completed");
}

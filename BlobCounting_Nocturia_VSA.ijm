macro "Automated-Particle-Analysis" {
waitForUser("Choose INPUT folder");	
inputFolder=getDirectory("Choose input folder");
waitForUser("Choose OUTPUT folder");
outputFolder=getDirectory("Choose output folder");
list=getFileList(inputFolder);

//setBatchMode(true);

for(i=0; i<list.length; i++) {
 path=inputFolder+list[i];
 name=list[i];
 if(endsWith(path,".tiff")) open(path);
 showProgress(i, list.length);
 if(nImages>=1) {
  if(i==0) {
	
  }
	setTool("line");
	waitForUser("Draw the scale bar");
	run("Set Scale...", "known=10 unit=cm");
	setTool("polygon");
	waitForUser("Crop the arena");
	run("Crop");
	run("Clear Outside");
	run("8-bit");
	run("Set Measurements...", "area redirect=None decimal=3");	
	run("Duplicate...", "title=duplicate");
	//run("Brightness/Contrast...");
	run("Enhance Contrast", "saturated=0.1");
	run("Apply LUT");
	run("Gaussian Blur...", "sigma=1.0");
	//run("Subtract Background...", "rolling=15");
	run("Auto Threshold", "method=Triangle white");
	run("Analyze Particles...", "size=16-Infinity pixel circularity=0.1-1.00 show=Outlines display clear");
	
	//Record zeros. Found this fix online.
	currentNResults = nResults;
	if (nResults == currentNResults) {  // if Analyze Particles did not find any particles
    setResult("Label", nResults, getTitle());
	}
	
  //output path
  outputPath=outputFolder+list[i];
  	//save the mask
	selectWindow("duplicate");
	//saveAs("TIFF", outputPath + list[i] + ".tif");
	saveAs("TIFF", outputPath  + ".tif");
	//run("Close");
	//save the outlines
	selectWindow("Drawing of duplicate");
	//saveAs("TIFF", outputPath + list[i] + ".tif");
	saveAs("TIFF", outputPath  + "_outlines.tif");
	//run("Close");
  //save the results
  selectWindow("Results");
  //outputPath=outputFolder+list[i];
  //The following two lines removes the file extension
  fileExtension=lastIndexOf(outputPath,"."); 
  if(fileExtension!=-1) outputPath=substring(outputPath,0,fileExtension);
  saveAs("Results", outputPath+".csv");
  //run("Close"); //closes Results window
  run("Clear Results");  //testing now
  close(); //closes the current image
  }
 }
//setBatchMode(false);
}